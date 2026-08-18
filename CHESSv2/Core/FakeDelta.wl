(* ::Package:: *)

(*
  FakeDelta.wl

  This module handles the special equation

      dY(t)/dt = B0(t) Y(t)

  by introducing an auxiliary parameter delta,

      dY(t,delta)/dt = delta B0(t) Y(t,delta),
      Y(t,delta) = Sum[k>=0] delta^k C_k(t).

  The coefficient equations are exactly canonical:

      C_0' = 0,
      C_k' = B0 C_(k-1).

  Consequently the original canonical CHESS solver can be reused unchanged.
  Only the construction of the auxiliary boundary columns, the evaluation at
  delta=1, and the adaptive order selection are new.

  Important numerical status:
  the ratio-based tail estimate below is an empirical stopping rule.  It is
  useful and was validated in the AMFlow transport experiments, but it is not
  a proof that every future coefficient obeys the observed geometric envelope.
*)

ClearAll[
  CHESSFakeDeltaRunAtOrder,
  CHESSDeltaCoefficientMatrices,
  CHESSDeltaCoefficientMetrics,
  CHESSDeltaSafeRatio,
  CHESSDeltaTailEstimates,
  CHESSFindAcceptedDeltaOrder,
  CHESSFakeDeltaAssembleResult,
  CHESSFakeDeltaPropagateSingle,
  CHESSFakeDeltaPropagate
];

CHESSFakeDeltaPropagate::boundary =
  "The B0-only boundary must be a non-empty numerical matrix; received dimensions `1`.";
CHESSFakeDeltaPropagate::order =
  "DeltaOrder must be Automatic or a non-negative integer; received `1`.";
CHESSFakeDeltaPropagate::maxorder =
  "MaxDeltaOrder must be an integer not smaller than four for automatic order selection; received `1`.";
CHESSFakeDeltaPropagate::tolerance =
  "DeltaTolerance must be Automatic or a positive number; received `1`.";
CHESSFakeDeltaPropagate::segments =
  "Segments must be a positive integer; received `1`.";
CHESSFakeDeltaPropagate::solve =
  "The canonical auxiliary-delta solve failed at expansion order `1`.";
CHESSFakeDeltaPropagate::noconvergence =
  "No acceptable fake-delta tail was found through order `1`; the last empirical tail estimate was `2` with target `3`.";

Options[CHESSFakeDeltaPropagate] = {
  "DeltaOrder" -> Automatic,
  "DeltaTolerance" -> Automatic,
  "MaxDeltaOrder" -> 128,
  "Segments" -> 1
};

(* Run the unchanged canonical propagator through one requested delta order.
   If the physical boundary has q columns, each column is transported in a
   separate auxiliary series.  B0 does not mix physical epsilon columns, so the
   results can be recombined column by column at delta=1.  This simple approach
   avoids constructing an artificial Kronecker-product system. *)
CHESSFakeDeltaRunAtOrder[
  b0_, boundary_?MatrixQ, interval_, order_Integer?NonNegative,
  canonicalRules_List
] := Module[
  {n, physicalColumns, evaluator, fakeBoundary, results},
  n = Dimensions[boundary][[1]];
  physicalColumns = Dimensions[boundary][[2]];
  (* Validate B0 against the physical system size before it reaches the copied
     canonical coefficient loop, whose historical failure propagation is not
     reliable for a dimension-mismatched operator. *)
  evaluator = CHESSCanonicalNodeEvaluator[b0, {n, n}];

  results = Table[
    fakeBoundary = ConstantArray[0, {n, order + 1}];
    fakeBoundary[[All, 1]] = boundary[[All, column]];
    CHESSCheckedCanonicalPropagate[
      evaluator,
      fakeBoundary,
      interval,
      canonicalRules
    ],
    {column, physicalColumns}
  ];

  If[MemberQ[results, $Failed], $Failed, results]
];

(* Convert q canonical endpoint matrices of shape n x (K+1) into a list

      {C_0, C_1, ..., C_K},

   where every C_k is n x q and therefore has exactly the shape of the physical
   boundary.  The same coefficient list drives both the tail estimator and the
   final delta=1 sum. *)
CHESSDeltaCoefficientMatrices[results_List] := Module[
  {orderCount},
  orderCount = Dimensions[results[[1, 3]]][[2]];
  Table[
    Transpose[results[[All, 3, All, orderIndex]]],
    {orderIndex, orderCount}
  ]
];

(* The scale-invariant coefficient metric used in the AMFlow order scans is

      a_k = max_i |C_k(i)| / (1 + |S_k(i)|),
      S_k = C_0 + ... + C_k.

   Dividing by 1+|S_k| avoids declaring a large component converged merely
   because all coefficients share a large overall normalization. *)
CHESSDeltaCoefficientMetrics[coefficientMatrices_List] := Module[
  {partial, metrics, coefficient, scaled},
  partial = ConstantArray[0, Dimensions[First[coefficientMatrices]]];
  metrics = Table[
    coefficient = coefficientMatrices[[index]];
    partial = partial + coefficient;
    scaled = Abs[Flatten[coefficient]]/(1 + Abs[Flatten[partial]]);
    If[scaled === {}, 0, Max[scaled]],
    {index, Length[coefficientMatrices]}
  ];
  metrics
];

(* Define 0/0 as zero because two consecutive exactly vanishing coefficients do
   not signal growth.  A nonzero coefficient following an exact zero produces
   Infinity and therefore cannot pass the convergence gate. *)
CHESSDeltaSafeRatio[new_, old_] := Which[
  TrueQ[PossibleZeroQ[old]] && TrueQ[PossibleZeroQ[new]], 0,
  TrueQ[PossibleZeroQ[old]], Infinity,
  True, new/old
];

(* For k>=3, envelope the last three observed ratios and estimate the omitted
   tail by

      E_k = 2 a_k q_k/(1-q_k),  0 <= q_k < 1.

   Earlier orders cannot form three ratios and are assigned Infinity.  The
   factor two is the same conservative empirical margin used in the verified
   AMFlow scan. *)
CHESSDeltaTailEstimates[metrics_List] := Module[
  {estimates, ratios, q},
  estimates = ConstantArray[Infinity, Length[metrics]];
  Do[
    ratios = {
      CHESSDeltaSafeRatio[metrics[[index]], metrics[[index - 1]]],
      CHESSDeltaSafeRatio[metrics[[index - 1]], metrics[[index - 2]]],
      CHESSDeltaSafeRatio[metrics[[index - 2]], metrics[[index - 3]]]
    };
    q = Max[ratios];
    estimates[[index]] = If[
      NumberQ[q] && FreeQ[q, _DirectedInfinity] && 0 <= q < 1,
      2 metrics[[index]] q/(1 - q),
      Infinity
    ],
    {index, 4, Length[metrics]}
  ];
  estimates
];

(* Find the first order whose estimate is below tolerance and whose immediately
   following guard layer is also below tolerance.  List index i corresponds to
   delta order i-1; accepting crossing index i and guard index i+1 therefore
   returns physical truncation order i. *)
CHESSFindAcceptedDeltaOrder[coefficientMatrices_List, tolerance_] := Module[
  {metrics, estimates, crossing},
  metrics = CHESSDeltaCoefficientMetrics[coefficientMatrices];
  estimates = CHESSDeltaTailEstimates[metrics];
  crossing = SelectFirst[
    Range[4, Length[estimates] - 1],
    FreeQ[estimates[[#]], _DirectedInfinity] &&
      FreeQ[estimates[[# + 1]], _DirectedInfinity] &&
      estimates[[#]] < tolerance && estimates[[# + 1]] < tolerance &,
    Missing["NotFound"]
  ];
  If[
    MissingQ[crossing],
    {Missing["NotFound"], Last[estimates]},
    {crossing, estimates[[crossing + 1]]}
  ]
];

(* Repackage the accepted auxiliary series into the ordinary CHESS six-part
   result.  Dummy coefficients are summed internally at delta=1; the public
   final state therefore has the same n x q dimensions as the physical input.

   Node values use the historical coefficient-first flattening convention:
   for q physical columns, each node stores column 1's n components, then
   column 2's n components, and so on. *)
CHESSFakeDeltaAssembleResult[
  results_List, selectedOrder_Integer?NonNegative, tailEstimate_, tolerance_,
  selectionMode_, computedOrder_Integer?NonNegative
] := Module[
  {n, physicalColumns, coefficientMatrices, finalState, nodeCount,
   nodeMatrices, coefficientRows, columnVectors, nodeValues},
  n = Dimensions[results[[1, 3]]][[1]];
  physicalColumns = Length[results];
  coefficientMatrices = CHESSDeltaCoefficientMatrices[results];
  finalState = Total[Take[coefficientMatrices, selectedOrder + 1]];
  nodeCount = Length[results[[1, 1]]];

  nodeMatrices = Table[
    columnVectors = Table[
      coefficientRows = Partition[results[[column, 2, node]], n];
      Total[Take[coefficientRows, selectedOrder + 1]],
      {column, physicalColumns}
    ];
    Transpose[columnVectors],
    {node, nodeCount}
  ];
  nodeValues = Flatten[Transpose[#], 1] & /@ nodeMatrices;

  {
    results[[1, 1]],
    nodeValues,
    finalState,
    {},
    {},
    {
      "FakeDelta",
      selectionMode,
      "SelectedOrder", selectedOrder,
      "ComputedOrder", computedOrder,
      "TailEstimate", tailEstimate,
      "Tolerance", tolerance,
      "GuardLayers", If[selectionMode === "Automatic", 1, 0]
    }
  }
];

(* Solve one interval.  Automatic mode grows the available coefficient range by
   doubling 8,16,32,... until the first crossing plus guard layer is present.
   Recomputing at logarithmically many orders keeps this wrapper small and leaves
   the verified canonical solver untouched; the largest solve dominates cost. *)
CHESSFakeDeltaPropagateSingle[
  b0_, boundary_?MatrixQ, interval_, canonicalRules_List,
  deltaOrder_, tolerance_, maxOrder_Integer
] := Module[
  {currentOrder, results, coefficientMatrices, accepted, selectedOrder,
   tailEstimate, metrics, estimates},

  If[IntegerQ[deltaOrder],
    results = CHESSFakeDeltaRunAtOrder[
      b0, boundary, interval, deltaOrder, canonicalRules
    ];
    If[results === $Failed,
      Message[CHESSFakeDeltaPropagate::solve, deltaOrder];
      Return[$Failed]
    ];
    coefficientMatrices = CHESSDeltaCoefficientMatrices[results];
    metrics = CHESSDeltaCoefficientMetrics[coefficientMatrices];
    estimates = CHESSDeltaTailEstimates[metrics];
    tailEstimate = If[deltaOrder + 1 <= Length[estimates],
      estimates[[deltaOrder + 1]], Infinity
    ];
    Return @ CHESSFakeDeltaAssembleResult[
      results, deltaOrder, tailEstimate, tolerance, "Fixed", deltaOrder
    ]
  ];

  currentOrder = Min[8, maxOrder];
  While[True,
    results = CHESSFakeDeltaRunAtOrder[
      b0, boundary, interval, currentOrder, canonicalRules
    ];
    If[results === $Failed,
      Message[CHESSFakeDeltaPropagate::solve, currentOrder];
      Return[$Failed]
    ];
    coefficientMatrices = CHESSDeltaCoefficientMatrices[results];
    accepted = CHESSFindAcceptedDeltaOrder[coefficientMatrices, tolerance];
    selectedOrder = accepted[[1]];
    tailEstimate = accepted[[2]];
    If[!MissingQ[selectedOrder],
      Return @ CHESSFakeDeltaAssembleResult[
        results, selectedOrder, tailEstimate, tolerance,
        "Automatic", currentOrder
      ]
    ];
    If[currentOrder >= maxOrder,
      Message[CHESSFakeDeltaPropagate::noconvergence,
        maxOrder, tailEstimate, tolerance];
      Return[$Failed]
    ];
    currentOrder = Min[2 currentOrder, maxOrder]
  ]
];

(* Public internal entry for B0-only dispatch.  Segmentation is implemented here
   rather than modifying either historical core.  Every segment receives the
   physical delta=1 state from the previous segment and chooses its own order.
   This is the correct composition law for the actual B0 transport. *)
CHESSFakeDeltaPropagate[
  b0_, boundary_?MatrixQ, {x0_, x1_}, canonicalRules_List,
  OptionsPattern[]
] := Module[
  {deltaOrder, tolerance, maxOrder, segments, precision, segmentPoints,
   state, result, allNodes, allValues, segmentInfo, segmentFailure},

  If[Dimensions[boundary][[1]] <= 0 || Dimensions[boundary][[2]] <= 0 ||
      !MatrixQ[boundary, NumericQ],
    Message[CHESSFakeDeltaPropagate::boundary, Dimensions[boundary]];
    Return[$Failed]
  ];

  deltaOrder = OptionValue["DeltaOrder"];
  If[!(deltaOrder === Automatic || IntegerQ[deltaOrder] && deltaOrder >= 0),
    Message[CHESSFakeDeltaPropagate::order, deltaOrder];
    Return[$Failed]
  ];

  maxOrder = OptionValue["MaxDeltaOrder"];
  If[deltaOrder === Automatic && (!IntegerQ[maxOrder] || maxOrder < 4),
    Message[CHESSFakeDeltaPropagate::maxorder, maxOrder];
    Return[$Failed]
  ];

  precision = Replace["Precision" /. canonicalRules, "Precision" -> 160];
  tolerance = Replace[
    OptionValue["DeltaTolerance"],
    Automatic :> 10^-Max[10, Floor[7 N[precision]/10]]
  ];
  If[!NumberQ[tolerance] || !TrueQ[tolerance > 0],
    Message[CHESSFakeDeltaPropagate::tolerance,
      OptionValue["DeltaTolerance"]];
    Return[$Failed]
  ];

  segments = OptionValue["Segments"];
  If[!IntegerQ[segments] || segments <= 0,
    Message[CHESSFakeDeltaPropagate::segments, segments];
    Return[$Failed]
  ];

  (* Dispatch has already filtered canonicalRules, because the historical
     canonical core must never receive the fake-delta Segments option. *)
  segmentPoints = N[Subdivide[x0, x1, segments], precision];
  state = boundary;
  allNodes = {};
  allValues = {};
  segmentInfo = {};
  segmentFailure = False;

  Do[
    result = CHESSFakeDeltaPropagateSingle[
      b0,
      state,
      {segmentPoints[[segment]], segmentPoints[[segment + 1]]},
      canonicalRules,
      deltaOrder,
      tolerance,
      maxOrder
    ];
    If[result === $Failed,
      segmentFailure = True;
      Break[]
    ];
    If[segment === 1,
      allNodes = result[[1]];
      allValues = result[[2]],
      allNodes = Join[allNodes, Rest[result[[1]]]];
      allValues = Join[allValues, Rest[result[[2]]]]
    ];
    state = result[[3]];
    AppendTo[segmentInfo, result[[6]]],
    {segment, 1, segments}
  ];
  If[TrueQ[segmentFailure], Return[$Failed]];

  {
    allNodes,
    allValues,
    state,
    {},
    {},
    If[segments === 1, First[segmentInfo],
      {"FakeDeltaSegments", "Segments", segments, "SegmentInfo", segmentInfo}
    ]
  }
];
