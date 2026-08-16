(* ::Package:: *)

(* Non-canonical epsilon-polynomial extension for CHESS.
   The implementation starts from the canonical sequential-epsilon solver:
   Chebyshev nodes and scalar collocation solves are reused directly from the
   canonical Chess.wl package distributed beside this file.

   It solves systems of the form
     J'(t,eps) = Sum[p>=0] eps^p B_p(t) J(t,eps)
   for the epsilon coefficients of J.  If B_0 is zero, the code reduces to the
   canonical CHESS solve.  If B_0 is sparse, only the components touched by B_0
   are solved as a coupled Chebyshev system; all other components keep the fast
   scalar collocation path. *)

(* CHESSv2 keeps this numerical algorithm intact.  Its only behavioral repair
   below is explicit failure propagation out of Do loops: a failed later path
   segment or epsilon layer must return $Failed from the enclosing solver,
   rather than accidentally returning a partially filled success result. *)

ChessNonCanonical::base =
  "The canonical CHESS core was not found at `1`. Load Core/Canonical.wl before Core/NonCanonical.wl.";

(* Load the canonical implementation from the extension's own directory. *)
If[!ValueQ[$CHESSVersion],
  With[
    {
      packageDirectory = If[
        StringQ[$InputFileName] && $InputFileName =!= "",
        DirectoryName[ExpandFileName[$InputFileName]],
        Directory[]
      ]
    },
    With[{basePackage = FileNameJoin[{packageDirectory, "Canonical.wl"}]},
      If[FileExistsQ[basePackage],
        Get[basePackage],
        Message[ChessNonCanonical::base, FileNameJoin[{packageDirectory, "Canonical.wl"}]];
        Abort[]
      ]
    ]
  ]
];

(* Stop immediately if an incomplete or incompatible base package was loaded. *)
If[!ValueQ[$CHESSVersion] || DownValues[ChebyshevLobattoData] === {},
  Message[
    ChessNonCanonical::base,
    FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "Canonical.wl"}]
  ];
  Abort[]
];

ClearAll[
  SpectralPropagatePolynomialEpsilon,
  CHESSNonCanonicalNumericSparseQ,
  CHESSNonCanonicalNumericRowsQ,
  CHESSNonCanonicalEvalMatrix,
  CHESSNonCanonicalEvalMatrixList,
  CHESSNonCanonicalPrepareParallel,
  CHESSNonCanonicalActiveSupport,
  CHESSNonCanonicalEndpointLookup,
  CHESSNonCanonicalEndpointSupport,
  CHESSNonCanonicalBlockRules,
  CHESSNonCanonicalCoupledMatrix,
  CHESSNonCanonicalCoupledSolve,
  CHESSNonCanonicalBuildRHS,
  CHESSNonCanonicalEndpointSource,
  CHESSNonCanonicalBuildOperatorMatrix
];

SpectralPropagatePolynomialEpsilon::badfun =
  "The coefficient matrix evaluator failed at node `1`.";
SpectralPropagatePolynomialEpsilon::badbatch =
  "The batch coefficient-matrix evaluator returned `1`; expected data for `2` ordinary nodes.";
SpectralPropagatePolynomialEpsilon::baddegree =
  "The matrix evaluator returned only `1` epsilon-power matrices, but degree `2` was requested.";
SpectralPropagatePolynomialEpsilon::solvefail =
  "Linear solve failed at epsilon coefficient `1`.";
SpectralPropagatePolynomialEpsilon::endpoint =
  "Endpoint regularization is not implemented in the clean non-canonical solver.";
SpectralPropagatePolynomialEpsilon::badnodes =
  "Nodes must be a positive integer; received `1`.";

Options[SpectralPropagatePolynomialEpsilon] = {
  "Nodes" -> 48,
  "Precision" -> 160,
  "WorkingPrecisionA" -> 220,
  "LinearSolverMethod" -> Automatic,
  "Segments" -> 1,
  "EpsilonDegree" -> Automatic,
  "RegularizedEndpoints" -> {},
  "EndpointData" -> {},
  "ParallelEvaluation" -> Automatic,
  "ParallelKernels" -> 1,
  "ParallelDefinitionsReady" -> False
};

CHESSNonCanonicalNumericSparseQ[m_] := Module[{rules},
  rules = Most[ArrayRules[SparseArray[m]]];
  rules === {} || VectorQ[rules[[All, 2]], NumericQ]
];

CHESSNonCanonicalNumericRowsQ[rows_] := VectorQ[Flatten[rows], NumericQ];

CHESSNonCanonicalEvalMatrix[f_, tvalue_, prec_] := Module[{value},
  value = Quiet @ Check[f[tvalue, prec], $Failed];
  If[value === $Failed, value = Quiet @ Check[f[tvalue], $Failed]];
  If[value === $Failed, Return[$Failed]];
  value = SparseArray[N[value, prec]];
  If[CHESSNonCanonicalNumericSparseQ[value], value, $Failed]
];

CHESSNonCanonicalEvalMatrixList[spec_List, tvalue_, prec_] := Module[{mats},
  mats = CHESSNonCanonicalEvalMatrix[#, tvalue, prec] & /@ spec;
  If[MemberQ[mats, $Failed], $Failed, mats]
];

CHESSNonCanonicalEvalMatrixList[spec_CHESSNodeEvaluator, tvalue_, prec_] :=
  CHESSNonCanonicalEvalMatrixList[CHESSNodeEvaluatorScalar[spec], tvalue, prec];

CHESSNonCanonicalEvalMatrixList[spec_, tvalue_, prec_] := Module[{value},
  value = Quiet @ Check[spec[tvalue, prec], $Failed];
  If[value === $Failed, value = Quiet @ Check[spec[tvalue], $Failed]];
  If[value === $Failed || !ListQ[value], Return[$Failed]];
  value = SparseArray[N[#, prec]] & /@ value;
  If[AllTrue[value, CHESSNonCanonicalNumericSparseQ], value, $Failed]
];

CHESSNonCanonicalPrepareParallel[matrixSpec_] := Module[{dir},
  dir = Directory[];
  Quiet @ Check[ParallelEvaluate[SetDirectory[dir]], Null];
  Quiet @ Check[
    DistributeDefinitions[
      CHESSNonCanonicalEvalMatrix,
      CHESSNonCanonicalEvalMatrixList
    ],
    Null
  ];
  Quiet @ Check[DistributeDefinitions[matrixSpec], Null]
];

CHESSNonCanonicalActiveSupport[b0Nodes_] := Module[
  {rules, positions},
  rules = Flatten[Most[ArrayRules[SparseArray[#]]] & /@ b0Nodes, 1];
  If[rules === {}, Return[{}]];
  positions = rules[[All, 1]];
  Sort @ DeleteDuplicates @ Flatten[positions]
];

CHESSNonCanonicalEndpointLookup[point_, endpointData_, degree_Integer, n_Integer, prec_] := Module[
  {entry, data},
  entry = FirstCase[
    Flatten[{endpointData}, 1],
    (Rule[p_, d_] | {p_, d_}) /; Quiet[PossibleZeroQ[N[p - point, prec]]] :> d,
    $Failed
  ];
  If[entry === $Failed || Length[entry] < degree + 1, Return[$Failed]];
  data = ({SparseArray[N[#[[1]], prec]], SparseArray[N[#[[2]], prec]]} &) /@ Take[entry, degree + 1];
  If[
    !And @@ (Dimensions[#[[1]]] == {n, n} && Dimensions[#[[2]]] == {n, n} & /@ data),
    Return[$Failed]
  ];
  data
];

CHESSNonCanonicalEndpointSupport[endpointIndexData_] := Module[
  {rules},
  rules = Flatten[
    Table[
      Join[
        Most[ArrayRules[endpointIndexData[[e, 3, 1, 1]]]],
        Most[ArrayRules[endpointIndexData[[e, 3, 1, 2]]]]
      ],
      {e, Length[endpointIndexData]}
    ],
    1
  ];
  If[rules === {}, {}, Sort @ DeleteDuplicates @ Flatten[rules[[All, 1]]]]
];

CHESSNonCanonicalBlockRules[rowNode_, colNode_, block_, active_] := Module[
  {rules, s},
  s = Length[active];
  rules = Most[ArrayRules[SparseArray[block]]];
  ({(rowNode - 1) s + #[[1, 1]], (colNode - 1) s + #[[1, 2]]} -> #[[2]]) & /@ rules
];

CHESSNonCanonicalBuildOperatorMatrix[
  derivMat_, nodes_, x0_, b0Nodes_, active_, endpointIndexData_, leftQ_, prec_
] := Module[
  {m, s, rules, id, activeBlock, endpointIndices, endpointQ, endpointAt, data, r0, m0, sj},
  m = Length[derivMat] - 1;
  s = Length[active];
  id = SparseArray[Band[{1, 1}] -> 1, {s, s}];
  endpointIndices = endpointIndexData[[All, 1]];
  endpointQ[idx_] := MemberQ[endpointIndices, idx];
  endpointAt[idx_] := endpointIndexData[[FirstPosition[endpointIndices, idx][[1]], 3]];
  rules = Flatten[
    Flatten[
      Table[
        activeBlock = Which[
          TrueQ[leftQ] && j === 1,
            data = endpointAt[j];
            r0 = data[[1, 1]];
            If[l === 1, id - r0[[active, active]], SparseArray[{}, {s, s}]],
          TrueQ[leftQ],
            sj = nodes[[j]] - x0;
            derivMat[[j, l]] sj id + If[j === l, id - sj b0Nodes[[j]][[active, active]], SparseArray[{}, {s, s}]],
          j === 1,
            If[l === 1, id, SparseArray[{}, {s, s}]],
          endpointQ[j],
            data = endpointAt[j];
            r0 = data[[1, 1]];
            m0 = data[[1, 2]];
            derivMat[[j, l]] (id - r0[[active, active]]) - If[j === l, m0[[active, active]], SparseArray[{}, {s, s}]],
          True,
            derivMat[[j, l]] id - If[j === l, b0Nodes[[j]][[active, active]], SparseArray[{}, {s, s}]]
        ];
        CHESSNonCanonicalBlockRules[j, l, activeBlock, active],
        {j, 1, m + 1}, {l, 1, m + 1}
      ],
      2
    ],
    1
  ];
  N[SparseArray[rules, {(m + 1) s, (m + 1) s}], prec]
];

CHESSNonCanonicalCoupledMatrix[derivMat_, b0Nodes_, active_, prec_] := Module[
  {m, s, rules, id, activeBlock},
  m = Length[derivMat] - 1;
  s = Length[active];
  id = SparseArray[Band[{1, 1}] -> 1, {s, s}];
  rules = Join[
    Table[{a, a} -> 1, {a, 1, s}],
    Flatten[
      Table[
        activeBlock = If[
          j === l,
          derivMat[[j, l]] id - b0Nodes[[j]][[active, active]],
          derivMat[[j, l]] id
        ];
        CHESSNonCanonicalBlockRules[j, l, activeBlock, active],
        {j, 2, m + 1}, {l, 1, m + 1}
      ],
      2
    ]
  ];
  N[SparseArray[rules, {(m + 1) s, (m + 1) s}], prec]
];

CHESSNonCanonicalCoupledSolve[solver_, rhsRows_, active_, prec_] := Module[
  {rhs, sol},
  rhs = SetPrecision[rhsRows[[All, active]], prec];
  If[!CHESSNonCanonicalNumericRowsQ[rhs], Return[$Failed]];
  sol = Quiet @ Check[solver[Flatten[rhs]], $Failed];
  If[sol === $Failed, Return[$Failed]];
  Partition[SetPrecision[sol, prec], Length[active]]
];

CHESSNonCanonicalBuildRHS[bByPower_, coeffSolutions_, coeff_Integer, degree_Integer, n_Integer] := Module[
  {m, rhsRows, j, p, source},
  m = Length[bByPower[[1]]] - 1;
  rhsRows = ConstantArray[0, {m + 1, n}];
  Do[
    source = ConstantArray[0, n];
    Do[
      If[coeff - p >= 1,
        source = source + bByPower[[p + 1, j]] . coeffSolutions[[coeff - p, j]]
      ],
      {p, 1, Min[degree, coeff - 1]}
    ];
    rhsRows[[j]] = source,
    {j, 1, m + 1}
  ];
  rhsRows
];

CHESSNonCanonicalEndpointSource[
  endpointData_, coeffSolutions_, coeffDerivatives_, coeff_Integer, nodeIndex_Integer, degree_Integer, n_Integer
] := Module[
  {source, p},
  source = ConstantArray[0, n];
  Do[
    If[coeff - p >= 1,
      source = source +
        endpointData[[p + 1, 1]] . coeffDerivatives[[coeff - p, nodeIndex]] +
        endpointData[[p + 1, 2]] . coeffSolutions[[coeff - p, nodeIndex]]
    ],
    {p, 1, Min[degree, coeff - 1]}
  ];
  source
];

SpectralPropagatePolynomialEpsilon[matrixSpec_, y0_, {x0_, x1_}, opts : OptionsPattern[]] := Module[
  {
    m, prec, precA, solverMethod, segments, regularizedEndpoints, endpointDataSpec, degreeOpt,
    kernels, definitionsReady, collocation, nodes, derivMat, n, k, parallelQ,
    probe, degree, matrixAtNode, matricesByNode, bByPower, b0Nodes,
    scalarBaseBlock, scalarSolver, active, inactive, coupledMatrix, coupledSolver,
    coeffSolutions, coeff, rhsRows, activeRows, inactiveSolution, solutionRows,
    values, finalValue, segmentPoints, state, allNodes, allValues, segResult,
    segmentNodes, segmentValues, singleOpts, s, activeRegularizedEndpoints,
    endpointIndexData, endpointDatum, endpointIndices, leftEndpointRegularized,
    coeffDerivatives, zRows, sj, endpointDataAt, jpos, batchEvaluatorQ,
    batchThreadCount, ordinaryIndices, batchValues, zeroMatrixList, index,
    segmentFailure, coefficientFailure, batchFailure
  },
  m = OptionValue["Nodes"];
  m = Which[
    IntegerQ[m] && m > 0,
      m,
    NumericQ[m] && TrueQ[PossibleZeroQ[m - Round[m]]] && Round[m] > 0,
      Round[m],
    True,
      Message[SpectralPropagatePolynomialEpsilon::badnodes, m];
      Return[$Failed]
  ];
  prec = OptionValue["Precision"];
  precA = OptionValue["WorkingPrecisionA"];
  solverMethod = OptionValue["LinearSolverMethod"];
  segments = OptionValue["Segments"];
  degreeOpt = OptionValue["EpsilonDegree"];
  kernels = OptionValue["ParallelKernels"];
  definitionsReady = TrueQ[OptionValue["ParallelDefinitionsReady"]];
  regularizedEndpoints = NormalizeRegularizedEndpoints[OptionValue["RegularizedEndpoints"], {x0, x1}];
  endpointDataSpec = OptionValue["EndpointData"];

  If[IntegerQ[segments] && segments > 1,
    singleOpts = DeleteCases[
      {opts},
      Rule["Segments", _] |
      RuleDelayed["Segments", _] |
      Rule["ParallelDefinitionsReady", _] |
      RuleDelayed["ParallelDefinitionsReady", _]
    ];
    segmentPoints = N[Subdivide[x0, x1, segments], prec];
    state = y0;
    allNodes = {};
    allValues = {};
    segmentFailure = False;
    Do[
      segResult = SpectralPropagatePolynomialEpsilon[
        matrixSpec,
        state,
        {segmentPoints[[s]], segmentPoints[[s + 1]]},
        Sequence @@ Join[singleOpts, {"Segments" -> 1, "ParallelDefinitionsReady" -> definitionsReady}]
      ];
      (* Return inside Do does not reliably propagate through the enclosing
         Module.  Record the failure and test it after the loop so a failed
         later segment can never masquerade as a successful shorter run. *)
      If[segResult === $Failed,
        segmentFailure = True;
        Break[]
      ];
      segmentNodes = segResult[[1]];
      segmentValues = segResult[[2]];
      If[s === 1,
        allNodes = segmentNodes;
        allValues = segmentValues,
        allNodes = Join[allNodes, Rest[segmentNodes]];
        allValues = Join[allValues, Rest[segmentValues]]
      ];
      state = segResult[[3]],
      {s, 1, segments}
    ];
    If[TrueQ[segmentFailure], Return[$Failed]];
    Return[{allNodes, allValues, state, {}, {}, "PolynomialEpsilonSegments"}]
  ];

  collocation = ChebyshevLobattoData[m, {x0, x1}, prec];
  nodes = N[collocation[[1]], prec];
  derivMat = N[collocation[[2]], prec];
  n = Dimensions[y0][[1]];
  k = Dimensions[y0][[2]];
  leftEndpointRegularized = LeftEndpointRegularizedQ[regularizedEndpoints, x0];
  activeRegularizedEndpoints = If[
    TrueQ[leftEndpointRegularized],
    Select[regularizedEndpoints, EndpointPointIndex[#, nodes] > 0 &],
    ActiveRegularizedEndpoints[regularizedEndpoints, nodes]
  ];
  (* A requested point outside the actual collocation grid used to disappear
     from activeRegularizedEndpoints silently, which could also discard its
     EndpointData.  Reject every such request before constructing operators. *)
  If[Length[activeRegularizedEndpoints] =!= Length[regularizedEndpoints],
    Message[SpectralPropagatePolynomialEpsilon::endpoint];
    Return[$Failed]
  ];
  endpointIndices = EndpointPointIndex[#, nodes] & /@ activeRegularizedEndpoints;
  batchEvaluatorQ = CHESSNodeEvaluatorQ[matrixSpec];
  batchThreadCount = If[
    batchEvaluatorQ,
    CHESSBatchThreadCount[OptionValue["ParallelEvaluation"], m, kernels],
    1
  ];

  degree = Which[
    ListQ[matrixSpec],
      Length[matrixSpec] - 1,
    IntegerQ[degreeOpt] && degreeOpt >= 0,
      degreeOpt,
    True,
      probe = CHESSNonCanonicalEvalMatrixList[matrixSpec, nodes[[2]], precA];
      If[probe === $Failed, Return[$Failed]];
      Length[probe] - 1
  ];

  endpointIndexData = Table[
    endpointDatum = CHESSNonCanonicalEndpointLookup[point, endpointDataSpec, degree, n, precA];
    If[endpointDatum === $Failed,
      Message[SpectralPropagatePolynomialEpsilon::endpoint];
      $Failed,
      {EndpointPointIndex[point, nodes], point, endpointDatum}
    ],
    {point, activeRegularizedEndpoints}
  ];
  If[MemberQ[endpointIndexData, $Failed], Return[$Failed]];

  parallelQ = If[
    batchEvaluatorQ,
    False,
    If[
      Length[DownValues[CHESSResolveParallelQ]] > 0,
      CHESSResolveParallelQ[OptionValue["ParallelEvaluation"], m, kernels],
      TrueQ[OptionValue["ParallelEvaluation"]] && kernels =!= 1
    ]
  ];
  If[parallelQ && kernels =!= 1 && !definitionsReady,
    CHESSNonCanonicalPrepareParallel[matrixSpec]
  ];

  matrixAtNode[j_] := Module[{mats},
    If[j === 1 || MemberQ[endpointIndices, j],
      Return[ConstantArray[SparseArray[{}, {n, n}], degree + 1]]
    ];
    mats = CHESSNonCanonicalEvalMatrixList[matrixSpec, nodes[[j]], precA];
    If[mats === $Failed,
      Message[SpectralPropagatePolynomialEpsilon::badfun, j];
      Return[$Failed]
    ];
    If[Length[mats] < degree + 1,
      Message[SpectralPropagatePolynomialEpsilon::baddegree, Length[mats], degree];
      Return[$Failed]
    ];
    Take[mats, degree + 1]
  ];

  matricesByNode = If[
    batchEvaluatorQ,
    ordinaryIndices = Select[
      Range[2, m + 1],
      !MemberQ[endpointIndices, #] &
    ];
    batchValues = If[
      ordinaryIndices === {},
      {},
      Quiet @ Check[
        CHESSNodeEvaluatorBatch[matrixSpec][
          nodes[[ordinaryIndices]],
          precA,
          batchThreadCount
        ],
        $Failed
      ]
    ];
    If[batchValues === $Failed || !ListQ[batchValues] ||
        Length[batchValues] =!= Length[ordinaryIndices],
      Message[
        SpectralPropagatePolynomialEpsilon::badbatch,
        If[batchValues === $Failed, $Failed, Head[batchValues]],
        Length[ordinaryIndices]
      ];
      Return[$Failed]
    ];
    batchFailure = False;
    batchValues = Table[
      If[!ListQ[batchValues[[index]]] || Length[batchValues[[index]]] < degree + 1,
        Message[
          SpectralPropagatePolynomialEpsilon::baddegree,
          If[ListQ[batchValues[[index]]], Length[batchValues[[index]]], 0],
          degree
        ];
        batchFailure = True;
        $Failed,
        SparseArray[N[#, precA]] & /@ Take[batchValues[[index]], degree + 1]
      ],
      {index, Length[batchValues]}
    ];
    If[TrueQ[batchFailure], Return[$Failed]];
    zeroMatrixList = ConstantArray[SparseArray[{}, {n, n}], degree + 1];
    matricesByNode = ConstantArray[zeroMatrixList, m + 1];
    Do[
      matricesByNode[[ordinaryIndices[[index]]]] = batchValues[[index]],
      {index, Length[ordinaryIndices]}
    ];
    matricesByNode,
    If[parallelQ && kernels =!= 1,
      ParallelTable[matrixAtNode[j], {j, 1, m + 1}, Method -> "CoarsestGrained"],
      Table[matrixAtNode[j], {j, 1, m + 1}]
    ]
  ];
  If[!FreeQ[matricesByNode, $Failed], Return[$Failed]];
  bByPower = Table[matricesByNode[[j, p + 1]], {p, 0, degree}, {j, 1, m + 1}];
  b0Nodes = bByPower[[1]];

  scalarBaseBlock = N[
    If[
      TrueQ[leftEndpointRegularized],
      LiftedScalarCollocationMatrix[derivMat, nodes, x0],
      BaseScalarCollocationMatrix[derivMat]
    ],
    prec
  ];
  scalarSolver = CHESSBuildLinearSolver[scalarBaseBlock, solverMethod];
  If[Head[scalarSolver] =!= LinearSolveFunction, Return[$Failed]];

  active = Sort @ DeleteDuplicates @ Join[
    CHESSNonCanonicalActiveSupport[b0Nodes],
    CHESSNonCanonicalEndpointSupport[endpointIndexData]
  ];
  inactive = Complement[Range[n], active];
  coupledSolver = If[
    active === {},
    None,
    coupledMatrix = If[
      endpointIndexData === {},
      CHESSNonCanonicalCoupledMatrix[derivMat, b0Nodes, active, prec],
      CHESSNonCanonicalBuildOperatorMatrix[
        derivMat, nodes, x0, b0Nodes, active, endpointIndexData,
        leftEndpointRegularized, prec
      ]
    ];
    Quiet @ Check[LinearSolve[coupledMatrix], $Failed]
  ];
  If[active =!= {} && Head[coupledSolver] =!= LinearSolveFunction, Return[$Failed]];

  coeffSolutions = ConstantArray[0, {k, m + 1, n}];
  coeffDerivatives = ConstantArray[0, {k, m + 1, n}];
  coefficientFailure = False;
  Do[
    rhsRows = CHESSNonCanonicalBuildRHS[bByPower, coeffSolutions, coeff, degree, n];
    solutionRows = ConstantArray[0, {m + 1, n}];
    If[endpointIndexData =!= {},
      If[
        TrueQ[leftEndpointRegularized],
        endpointDataAt = endpointIndexData[[FirstPosition[endpointIndexData[[All, 1]], 1][[1]], 3]];
        rhsRows[[1]] =
          endpointDataAt[[1, 2]] . y0[[All, coeff]] +
          CHESSNonCanonicalEndpointSource[
            endpointDataAt, coeffSolutions, coeffDerivatives, coeff, 1, degree, n
          ];
        Do[
          rhsRows[[j]] = rhsRows[[j]] + b0Nodes[[j]] . y0[[All, coeff]],
          {j, 2, m + 1}
        ],
        Do[
          endpointDataAt = endpointIndexData[[jpos, 3]];
          rhsRows[[endpointIndexData[[jpos, 1]]]] =
            CHESSNonCanonicalEndpointSource[
              endpointDataAt,
              coeffSolutions,
              coeffDerivatives,
              coeff,
              endpointIndexData[[jpos, 1]],
              degree,
              n
            ],
          {jpos, Length[endpointIndexData]}
        ]
      ]
    ];

    If[!TrueQ[leftEndpointRegularized],
      rhsRows[[1]] = SetPrecision[y0[[All, coeff]], prec]
    ];
    If[!CHESSNonCanonicalNumericRowsQ[rhsRows],
      Message[SpectralPropagatePolynomialEpsilon::solvefail, coeff];
      coefficientFailure = True;
      Break[]
    ];

    If[active === {},
      solutionRows = CHESSScalarCollocationSolve[scalarSolver, rhsRows, n, solverMethod, prec];
      If[solutionRows === $Failed,
        Message[SpectralPropagatePolynomialEpsilon::solvefail, coeff];
        coefficientFailure = True;
        Break[]
      ],
      activeRows = CHESSNonCanonicalCoupledSolve[coupledSolver, rhsRows, active, prec];
      If[activeRows === $Failed,
        Message[SpectralPropagatePolynomialEpsilon::solvefail, coeff];
        coefficientFailure = True;
        Break[]
      ];
      solutionRows[[All, active]] = activeRows;
      If[inactive =!= {},
        If[!TrueQ[leftEndpointRegularized],
          rhsRows[[1, inactive]] = SetPrecision[y0[[inactive, coeff]], prec]
        ];
        inactiveSolution = CHESSScalarCollocationSolve[
          scalarSolver,
          rhsRows[[All, inactive]],
          Length[inactive],
          solverMethod,
          prec
        ];
        If[inactiveSolution === $Failed,
          Message[SpectralPropagatePolynomialEpsilon::solvefail, coeff];
          coefficientFailure = True;
          Break[]
        ];
        solutionRows[[All, inactive]] = inactiveSolution
      ]
    ];
    If[
      TrueQ[leftEndpointRegularized],
      zRows = solutionRows;
      coeffSolutions[[coeff]] = Table[
        SetPrecision[y0[[All, coeff]] + (nodes[[j]] - x0) zRows[[j]], prec],
        {j, 1, m + 1}
      ];
      coeffDerivatives[[coeff]] = Table[
        SetPrecision[zRows[[j]] + (nodes[[j]] - x0) (derivMat . zRows)[[j]], prec],
        {j, 1, m + 1}
      ],
      coeffSolutions[[coeff]] = SetPrecision[solutionRows, prec];
      coeffDerivatives[[coeff]] = SetPrecision[derivMat . solutionRows, prec]
    ],
    {coeff, 1, k}
  ];
  If[TrueQ[coefficientFailure], Return[$Failed]];

  values = Table[Flatten[coeffSolutions[[All, j]], 1], {j, 1, m + 1}];
  finalValue = Transpose[coeffSolutions[[All, -1]]];
  {
    nodes,
    values,
    finalValue,
    {},
    {},
    If[
      active === {}, "PolynomialEpsilonCanonicalCompatible", {"PolynomialEpsilonActiveSupport", Length[active]}
    ]
  }
];
