(* ::Package:: *)

(*
  Dispatch.wl

  This is the only layer which changes the user-facing propagation API.  The
  historical solvers remain available internally as

      SpectralPropagateSequential              canonical core
      SpectralPropagatePolynomialEpsilon       general non-canonical core

  while the public SpectralPropagate chooses between them according to the
  explicit epsilon-coefficient list:

      B1                         -> canonical (legacy syntax)
      {B0}                       -> B0-only fake delta
      {0,B1}                     -> canonical
      {B0,B1,...} or any other   -> original non-canonical algorithm

  The classification is structural.  We never sample a user function and call
  it zero merely because it vanishes at selected numerical points.
*)

ClearAll[
  CHESSCanonicalSegments,
  CHESSRunNativeCanonicalUnified,
  CHESSRunCanonicalUnified,
  CHESSRunNonCanonicalUnified,
  CHESSRunFakeDeltaUnified,
  CHESSNonzeroCoefficientPositions
];

(* Preserve the two historical option contracts before replacing the public
   SpectralPropagate definition and extending its option list. *)
$CHESSCanonicalOptions = Options[SpectralPropagate];
$CHESSNonCanonicalOptions = Options[SpectralPropagatePolynomialEpsilon];

ClearAll[SpectralPropagate];

SpectralPropagate::usage =
  "SpectralPropagate[matrixSpec,boundary,{x0,x1},opts] automatically selects canonical, B0-only fake-delta, or general non-canonical CHESS propagation.";
SpectralPropagate::boundary =
  "Boundary data must be a non-empty numerical n x k matrix; received dimensions `1`.";
SpectralPropagate::spec =
  "The explicit coefficient list must contain at least one B_p entry.";
SpectralPropagate::segments =
  "Segments must be a positive integer; received `1`.";
SpectralPropagate::canonicalsegments =
  "Canonical segmentation with singular endpoint regularization is not supported by the unchanged canonical core. Use one segment or split the regular path explicitly.";
SpectralPropagate::endpointdata =
  "EndpointData belongs to the general non-canonical solver and cannot be used on a route automatically classified as canonical.";
SpectralPropagate::orphanendpointdata =
  "EndpointData was supplied without an active RegularizedEndpoints specification on the general non-canonical route.";
SpectralPropagate::noncanonicalsegments =
  "Segmented general non-canonical propagation cannot repeat singular endpoint regularization independently on every segment. Use one segment or split only regular subpaths explicitly.";
SpectralPropagate::b0endpoint =
  "Automatic B0-only fake-delta propagation currently supports regular endpoints only. Use the general non-canonical representation when explicit endpoint data are required.";
SpectralPropagate::normalize =
  "The polynomial coefficient list could not be normalized to the boundary dimensions.";
SpectralPropagate::nativebackend =
  "NumericalBackend must be Mathematica, Native, or Automatic; received `1`.";
SpectralPropagate::nativeunsupported =
  "The native propagation backend currently supports one regular canonical interval only. Use NumericalBackend->Mathematica for segmented or endpoint-regularized propagation.";
SpectralPropagate::nativenoncanonical =
  "The native propagation backend does not yet implement a genuinely mixed or higher-degree non-canonical operator. This route remains on the Mathematica non-canonical core.";
SpectralPropagate::nativeresult =
  "ResultData must be Full or Endpoint; received `1`.";
SpectralPropagate::nativehandle =
  "The supplied NativeHandle is incompatible with this dimension, interval, node count, or precision.";

Options[SpectralPropagate] = DeleteDuplicatesBy[
  Join[
    $CHESSCanonicalOptions,
    $CHESSNonCanonicalOptions,
    {
      (* Automatic fake-delta selection.  An integer requests a fixed order. *)
      "DeltaOrder" -> Automatic,

      (* Automatic reserves roughly thirty percent of Precision as guard
         digits; callers may provide an explicit positive tolerance. *)
      "DeltaTolerance" -> Automatic,

      (* Automatic mode fails rather than returning an unconverged sum after
         this order.  The verified AMFlow example required order 56. *)
      "MaxDeltaOrder" -> 128,

      (* Native remains opt-in while the experiment is generalized.  Passing a
         prepared handle amortizes matrix evaluation and LU factorization. *)
      "NumericalBackend" -> "Mathematica",
      "NativeHandle" -> Automatic,

      (* Endpoint skips the large coefficient-by-node WSTP payload and is the
         high-throughput mode for repeated boundary propagation. *)
      "ResultData" -> "Full"
    }
  ],
  First
];

(* Return 1-based positions of coefficients which are not structurally zero.
   Position one is B0, position two is B1, and so on. *)
CHESSNonzeroCoefficientPositions[specs_List] := Flatten @ Position[
  CHESSStructuralZeroCoefficientQ /@ specs,
  False
];

(* The original canonical solver has no Segments option.  For regular paths we
   compose several unchanged canonical solves, passing the complete physical
   epsilon coefficient matrix from one segment to the next.  Endpoint
   regularization is deliberately rejected here because interpreting "Left"
   or "Right" independently on every subinterval would be mathematically
   wrong. *)
CHESSCanonicalSegments[
  evaluator_, boundary_?MatrixQ, {x0_, x1_}, canonicalRules_List,
  segments_Integer?Positive
] := Module[
  {precision, points, state, result, allNodes, allValues, segment,
   segmentFailure},
  precision = Replace["Precision" /. canonicalRules, "Precision" -> 160];
  points = N[Subdivide[x0, x1, segments], precision];
  state = boundary;
  allNodes = {};
  allValues = {};
  segmentFailure = False;

  Do[
    result = CHESSCheckedCanonicalPropagate[
      evaluator,
      state,
      {points[[segment]], points[[segment + 1]]},
      canonicalRules
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
    state = result[[3]],
    {segment, segments}
  ];
  If[TrueQ[segmentFailure], Return[$Failed]];
  {allNodes, allValues, state, {}, {}, {"CanonicalSegments", segments}}
];

(* Native canonical execution is kept in this thin wrapper.  It neither changes
   the copied canonical core nor leaks WSTP protocol details into public route
   classification. *)
CHESSRunNativeCanonicalUnified[
  evaluator_, systemSpec_, boundary_?MatrixQ, interval_, canonicalRules_List,
  segments_Integer, regularizedEndpoints_, nativeHandle_, resultData_
] := Module[{dimension, handle, tensor, run},
  If[
    segments =!= 1 ||
    !MemberQ[{None, False, {}}, regularizedEndpoints],
    Message[SpectralPropagate::nativeunsupported];
    Return[$Failed]
  ];
  If[!MemberQ[{"Full", "Endpoint"}, resultData],
    Message[SpectralPropagate::nativeresult, resultData];
    Return[$Failed]
  ];
  dimension = Dimensions[boundary][[1]];
  handle = If[
    nativeHandle === Automatic,
    CHESSNativePrepare[systemSpec, dimension, interval, canonicalRules],
    nativeHandle
  ];
  If[
    handle === $Failed ||
    !CHESSNativeHandleCompatibleQ[
      handle, systemSpec, dimension, interval, canonicalRules
    ],
    Message[SpectralPropagate::nativehandle];
    Return[$Failed]
  ];
  tensor = {Transpose[boundary]};
  run = CHESSNativeRun[handle, tensor, resultData === "Full"];
  If[run === $Failed, Return[$Failed]];
  CHESSNativeCanonicalResult[handle, run, resultData]
];

(* Execute a canonical route after filtering the unified option list.  Both
   explicit B1 coefficients and legacy bare functions have already been wrapped
   by the common adapter, which preserves their call conventions while checking
   the boundary-derived operator dimensions. *)
CHESSRunCanonicalUnified[
  evaluator_, systemSpec_, boundary_?MatrixQ, interval_, rules_List,
  segments_Integer,
  regularizedEndpoints_, endpointData_, numericalBackend_, nativeHandle_,
  resultData_
] := Module[{canonicalRules, resolvedBackend},
  If[endpointData =!= {},
    Message[SpectralPropagate::endpointdata];
    Return[$Failed]
  ];
  If[segments > 1 && !MemberQ[{None, False, {}}, regularizedEndpoints],
    Message[SpectralPropagate::canonicalsegments];
    Return[$Failed]
  ];
  canonicalRules = CHESSFilterOptionRules[rules, $CHESSCanonicalOptions];
  resolvedBackend = Replace[
    numericalBackend,
    Automatic :> If[AssociationQ[nativeHandle], "Native", "Mathematica"]
  ];
  If[!MemberQ[{"Mathematica", "Native"}, resolvedBackend],
    Message[SpectralPropagate::nativebackend, numericalBackend];
    Return[$Failed]
  ];
  If[
    resolvedBackend === "Native",
    CHESSRunNativeCanonicalUnified[
      evaluator, systemSpec, boundary, interval, canonicalRules, segments,
      regularizedEndpoints, nativeHandle, resultData
    ],
    If[segments === 1,
      CHESSCheckedCanonicalPropagate[
        evaluator, boundary, interval, canonicalRules
      ],
      CHESSCanonicalSegments[
        evaluator, boundary, interval, canonicalRules, segments
      ]
    ]
  ]
];

(* Execute the original non-canonical numerical algorithm.  The normalization
   layer wraps zero/constant entries and genuine coefficient functions alike so
   all three call conventions and the boundary-derived dimensions are enforced
   before values enter the historical solver. *)
CHESSRunNonCanonicalUnified[
  specs_, boundary_?MatrixQ, interval_, rules_List, segments_Integer,
  regularizedEndpoints_, endpointData_, numericalBackend_, nativeHandle_
] := Module[{normalized, nonCanonicalRules, resolvedBackend},
  resolvedBackend = Replace[
    numericalBackend,
    Automatic :> If[AssociationQ[nativeHandle], "Native", "Mathematica"]
  ];
  If[!MemberQ[{"Mathematica", "Native"}, resolvedBackend],
    Message[SpectralPropagate::nativebackend, numericalBackend];
    Return[$Failed]
  ];
  If[resolvedBackend === "Native",
    Message[SpectralPropagate::nativenoncanonical];
    Return[$Failed]
  ];
  If[endpointData =!= {} &&
      MemberQ[{None, False, {}}, regularizedEndpoints],
    Message[SpectralPropagate::orphanendpointdata];
    Return[$Failed]
  ];
  If[segments > 1 &&
      !MemberQ[{None, False, {}}, regularizedEndpoints],
    Message[SpectralPropagate::noncanonicalsegments];
    Return[$Failed]
  ];
  (* A CHESSNodeEvaluator with an explicit positive EpsilonDegree already
     returns the complete {B0,B1,...} list at each node.  Keep that batch
     evaluator intact so a FORM/FLINT multipoint call is made only once. *)
  normalized = If[
    CHESSNodeEvaluatorQ[specs],
    specs,
    CHESSNormalizePolynomialSpec[specs, boundary]
  ];
  If[normalized === $Failed,
    Message[SpectralPropagate::normalize];
    Return[$Failed]
  ];
  nonCanonicalRules = CHESSFilterOptionRules[
    rules, $CHESSNonCanonicalOptions
  ];
  SpectralPropagatePolynomialEpsilon[
    normalized,
    boundary,
    interval,
    Sequence @@ nonCanonicalRules
  ]
];

(* Execute B0-only fake-delta propagation.  Explicit endpoint data cannot be
   inserted into the unchanged canonical endpoint machinery, so this new route
   is intentionally regular-point only. *)
CHESSRunFakeDeltaUnified[
  b0_, boundary_?MatrixQ, interval_, rules_List, segments_Integer,
  regularizedEndpoints_, endpointData_, deltaOrder_, deltaTolerance_,
  maxDeltaOrder_, numericalBackend_, nativeHandle_, resultData_
] := Module[
  {canonicalRules, resolvedBackend, precision, resolvedTolerance},
  (* EndpointData is meaningful only together with the non-canonical endpoint
     construction.  Silently dropping it here would be especially dangerous:
     the numerical answer could look plausible while ignoring user-supplied
     boundary information.  Until fake-delta endpoint matching is derived and
     tested, both endpoint controls therefore fail closed. *)
  If[endpointData =!= {} ||
      !MemberQ[{None, False, {}}, regularizedEndpoints],
    Message[SpectralPropagate::b0endpoint];
    Return[$Failed]
  ];
  canonicalRules = CHESSFilterOptionRules[rules, $CHESSCanonicalOptions];
  resolvedBackend = Replace[
    numericalBackend,
    Automatic :> If[AssociationQ[nativeHandle], "Native", "Mathematica"]
  ];
  If[!MemberQ[{"Mathematica", "Native"}, resolvedBackend],
    Message[SpectralPropagate::nativebackend, numericalBackend];
    Return[$Failed]
  ];
  If[
    resolvedBackend === "Native",
    If[!(deltaOrder === Automatic || IntegerQ[deltaOrder] && deltaOrder >= 0),
      Message[CHESSFakeDeltaPropagate::order, deltaOrder];
      Return[$Failed]
    ];
    If[
      !IntegerQ[maxDeltaOrder] || maxDeltaOrder < 0 ||
      deltaOrder === Automatic && maxDeltaOrder < 4,
      Message[CHESSFakeDeltaPropagate::maxorder, maxDeltaOrder];
      Return[$Failed]
    ];
    precision = CHESSNativeOptionValue[
      canonicalRules, "Precision", 160
    ];
    resolvedTolerance = Replace[
      deltaTolerance,
      Automatic :> 10^-Max[10, Floor[7 N[precision]/10]]
    ];
    If[!NumberQ[resolvedTolerance] || !TrueQ[resolvedTolerance > 0],
      Message[CHESSFakeDeltaPropagate::tolerance, deltaTolerance];
      Return[$Failed]
    ];
    CHESSNativeFakeDeltaPropagate[
      b0, boundary, interval, canonicalRules, deltaOrder,
      resolvedTolerance, maxDeltaOrder, segments, nativeHandle, resultData
    ],
    CHESSFakeDeltaPropagate[
      b0,
      boundary,
      interval,
      canonicalRules,
      "DeltaOrder" -> deltaOrder,
      "DeltaTolerance" -> deltaTolerance,
      "MaxDeltaOrder" -> maxDeltaOrder,
      "Segments" -> segments
    ]
  ]
];

(* Legacy syntax: a bare evaluator is B1 and therefore canonical.  A dense
   constant matrix also reaches this definition unless Mathematica first
   matches the List definition below; that definition detects MatrixQ and sends
   it back to this canonical route explicitly. *)
SpectralPropagate[
  evaluator_, boundary_?MatrixQ, interval : {_, _},
  opts : OptionsPattern[]
] := Module[
  {
    rules, segments, regularizedEndpoints, endpointData, canonicalEvaluator,
    numericalBackend, nativeHandle, resultData, epsilonDegree
  },
  If[Dimensions[boundary][[1]] <= 0 || Dimensions[boundary][[2]] <= 0 ||
      !MatrixQ[boundary, NumericQ],
    Message[SpectralPropagate::boundary, Dimensions[boundary]];
    Return[$Failed]
  ];
  rules = {opts};
  segments = OptionValue["Segments"];
  If[!IntegerQ[segments] || segments <= 0,
    Message[SpectralPropagate::segments, segments];
    Return[$Failed]
  ];
  regularizedEndpoints = OptionValue["RegularizedEndpoints"];
  endpointData = OptionValue["EndpointData"];
  numericalBackend = OptionValue["NumericalBackend"];
  nativeHandle = OptionValue["NativeHandle"];
  resultData = OptionValue["ResultData"];
  epsilonDegree = OptionValue["EpsilonDegree"];
  (* Batch evaluators which return polynomial coefficient matrices cannot be
     classified by their expression head.  An explicit positive degree is the
     unambiguous signal to retain the evaluator and use the non-canonical core.
     Without it, legacy CHESSNodeEvaluator objects remain canonical. *)
  If[
    CHESSNodeEvaluatorQ[evaluator] &&
    IntegerQ[epsilonDegree] && epsilonDegree >= 1,
    Return @ CHESSRunNonCanonicalUnified[
      evaluator, boundary, interval, rules, segments,
      regularizedEndpoints, endpointData, numericalBackend, nativeHandle
    ]
  ];
  (* SparseArray has head SparseArray rather than List, so a sparse constant
     reaches this generic definition.  Wrap every numerical matrix explicitly;
     otherwise the canonical core would try to call it as a function and can
     silently obtain zero source terms. *)
  canonicalEvaluator = If[
    CHESSNodeEvaluatorQ[evaluator],
    CHESSDimensionCheckedNodeEvaluator[
      evaluator, ConstantArray[Dimensions[boundary][[1]], 2]
    ],
    CHESSCanonicalNodeEvaluator[
      evaluator, ConstantArray[Dimensions[boundary][[1]], 2]
    ]
  ];
  CHESSRunCanonicalUnified[
    canonicalEvaluator, evaluator, boundary, interval, rules, segments,
    regularizedEndpoints, endpointData, numericalBackend, nativeHandle,
    resultData
  ]
];

(* Explicit polynomial syntax.  The only new decision logic in CHESSv2 lives in
   this definition; the actual numerical work remains in the three solver
   functions selected below. *)
SpectralPropagate[
  specs_List, boundary_?MatrixQ, interval : {_, _},
  opts : OptionsPattern[]
] := Module[
  {rules, segments, regularizedEndpoints, endpointData, positions,
   deltaOrder, deltaTolerance, maxDeltaOrder, zeroMatrix,
   numericalBackend, nativeHandle, resultData},

  (* A dense constant n x n matrix is legacy canonical input, not a list of
     epsilon coefficients. *)
  If[MatrixQ[specs] && ArrayDepth[specs] === 2,
    Return @ SpectralPropagate[
      (* SparseArray changes the expression head so the recursive call reaches
         the generic constant-matrix path exactly once. *)
      SparseArray[specs],
      boundary,
      interval,
      opts
    ]
  ];
  If[specs === {},
    Message[SpectralPropagate::spec];
    Return[$Failed]
  ];
  If[Dimensions[boundary][[1]] <= 0 || Dimensions[boundary][[2]] <= 0 ||
      !MatrixQ[boundary, NumericQ],
    Message[SpectralPropagate::boundary, Dimensions[boundary]];
    Return[$Failed]
  ];

  rules = {opts};
  segments = OptionValue["Segments"];
  If[!IntegerQ[segments] || segments <= 0,
    Message[SpectralPropagate::segments, segments];
    Return[$Failed]
  ];
  regularizedEndpoints = OptionValue["RegularizedEndpoints"];
  endpointData = OptionValue["EndpointData"];
  deltaOrder = OptionValue["DeltaOrder"];
  deltaTolerance = OptionValue["DeltaTolerance"];
  maxDeltaOrder = OptionValue["MaxDeltaOrder"];
  numericalBackend = OptionValue["NumericalBackend"];
  nativeHandle = OptionValue["NativeHandle"];
  resultData = OptionValue["ResultData"];
  positions = CHESSNonzeroCoefficientPositions[specs];

  Which[
    (* The identically zero equation is safely handled as B0=0. *)
    positions === {},
      zeroMatrix = SparseArray[{}, ConstantArray[Dimensions[boundary][[1]], 2]];
      CHESSRunFakeDeltaUnified[
        zeroMatrix, boundary, interval, rules, segments,
        regularizedEndpoints, endpointData, deltaOrder, deltaTolerance,
        maxDeltaOrder, numericalBackend, nativeHandle, resultData
      ],

    (* B0 is the only nonzero coefficient. *)
    positions === {1},
      CHESSRunFakeDeltaUnified[
        specs[[1]], boundary, interval, rules, segments,
        regularizedEndpoints, endpointData, deltaOrder, deltaTolerance,
        maxDeltaOrder, numericalBackend, nativeHandle, resultData
      ],

    (* B1 is the only nonzero coefficient; use the original canonical method. *)
    positions === {2},
      CHESSRunCanonicalUnified[
        CHESSCanonicalNodeEvaluator[
          specs[[2]], ConstantArray[Dimensions[boundary][[1]], 2]
        ],
        specs[[2]],
        boundary,
        interval,
        rules,
        segments,
        regularizedEndpoints,
        endpointData,
        numericalBackend,
        nativeHandle,
        resultData
      ],

    (* Every genuinely mixed or higher-degree case stays on the unchanged
       non-canonical active-support algorithm. *)
    True,
      CHESSRunNonCanonicalUnified[
        specs, boundary, interval, rules, segments,
        regularizedEndpoints, endpointData, numericalBackend, nativeHandle
      ]
  ]
];
