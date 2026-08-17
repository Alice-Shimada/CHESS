#!/usr/bin/env wolframscript

(* Assemble a fresh set of all-Mathematica F3 node matrices and propagate them
   through CHESSv2's unified entry.  Node generation is timed by the external
   worker orchestration; this script times only matrix import and propagation. *)

ClearAll["Global`*"];
$HistoryLength = 0;

If[!MemberQ[{6, 8}, Length[$ScriptCommandLine]],
  Print[
    "usage: f3_zoia_mathematica_assemble.wls EXAMPLE_ROOT NODES PRECISION ",
    "MATRIX_DIRECTORY OUTPUT_FILE ",
    "[OBSERVED_WORKER_WALL_SECONDS WORKER_TIME_LOG]"
  ];
  Quit[2]
];

packageRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
packageFile = FileNameJoin[{packageRoot, "CHESSv2.wl"}];
exampleRoot = ExpandFileName[$ScriptCommandLine[[2]]];
nodeCount = ToExpression[$ScriptCommandLine[[3]]];
workingPrecision = ToExpression[$ScriptCommandLine[[4]]];
matrixDirectory = ExpandFileName[$ScriptCommandLine[[5]]];
outputFile = ExpandFileName[$ScriptCommandLine[[6]]];
observedWorkerWallSeconds = If[
  Length[$ScriptCommandLine] === 8,
  ToExpression[$ScriptCommandLine[[7]]],
  Missing["NotRecorded"]
];
workerTimeLog = If[
  Length[$ScriptCommandLine] === 8,
  ExpandFileName[$ScriptCommandLine[[8]]],
  Missing["NotRecorded"]
];
If[
  (!MissingQ[observedWorkerWallSeconds] &&
    !(NumericQ[observedWorkerWallSeconds] && observedWorkerWallSeconds > 0)) ||
  (!MissingQ[workerTimeLog] && !FileExistsQ[workerTimeLog]),
  Print["CHESSV2_F3_MMA_WALL_TIME_INPUT_FAILED=True"];
  Quit[2]
];
If[!MissingQ[workerTimeLog],
  workerTimeLogText = Quiet @ Check[Import[workerTimeLog, "Text"], $Failed];
  workerTimeLogValues = If[
    StringQ[workerTimeLogText],
    StringCases[
      workerTimeLogText,
      "CHESSV2_F3_MMA_WORKER_WALL_SECONDS=" ~~ value : NumberString :>
        ToExpression[value]
    ],
    {}
  ];
  If[
    Length[workerTimeLogValues] =!= 1 ||
    !TrueQ[Abs[First[workerTimeLogValues] - observedWorkerWallSeconds] < 10^-6],
    Print["CHESSV2_F3_MMA_WALL_TIME_EVIDENCE_FAILED=True"];
    Quit[2]
  ]
];
projectRoot = DirectoryName[exampleRoot];
boundaryFile = FileNameJoin[{
  projectRoot, "ttW_F3", "ancillary_files", "F3", "values", "F3_X1.m"
}];
referenceFile = FileNameJoin[{
  projectRoot, "ttW_F3", "ancillary_files", "F3", "values", "F3_X2.m"
}];

requiredFiles = {packageFile, boundaryFile, referenceFile};
If[
  !And @@ (FileExistsQ /@ requiredFiles) || !DirectoryQ[matrixDirectory],
  Print["CHESSV2_F3_MMA_REQUIRED_INPUT_FAILED=True"];
  Quit[2]
];
If[FileExistsQ[outputFile],
  Print["CHESSV2_F3_MMA_REFUSE_OVERWRITE=", outputFile];
  Quit[2]
];

setupResult = Quiet @ Check[
  Get[packageFile];
  {initialPoint, initialValues} = Get[boundaryFile];
  {finalPoint, finalValues} = Get[referenceFile];
  True,
  $Failed
];
If[setupResult === $Failed,
  Print["CHESSV2_F3_MMA_SETUP_FAILED=True"];
  Quit[2]
];

importSeconds = AbsoluteTiming[
  matrixFiles = FileNames["matrix_*.m", matrixDirectory];
  entries = Quiet @ Check[Get /@ matrixFiles, $Failed];
  If[entries =!= $Failed, entries = SortBy[entries, First]];
][[1]];
If[entries === $Failed,
  Print["CHESSV2_F3_MMA_MATRIX_IMPORT_FAILED=True"];
  Quit[3]
];
indices = If[entries === {}, {}, Quiet @ Check[entries[[All, 1]], $Failed]];
If[indices =!= Range[2, nodeCount + 1],
  Print["CHESSV2_F3_MMA_MATRIX_SET_FAILED=", InputForm[indices]];
  Quit[3]
];

boundaryValues = N[
  Table[
    Coefficient[initialValues[[index]], eps, order],
    {index, Length[initialValues]},
    {order, 0, 4}
  ],
  workingPrecision
];
referenceValues = N[
  Table[
    Coefficient[finalValues[[index]], eps, order],
    {index, Length[finalValues]},
    {order, 0, 4}
  ],
  workingPrecision
];
If[
  !And @@ Map[
    Function[entry,
      ListQ[entry] && Length[entry] === 3 && IntegerQ[entry[[1]]] &&
      NumericQ[entry[[2]]] && entry[[2]] >= 0 &&
      ListQ[entry[[3]]] && Length[entry[[3]]] === 3 &&
      And @@ (
        MatrixQ[#, NumericQ] &&
        Dimensions[#] === ConstantArray[Length[boundaryValues], 2] & /@
          entry[[3]]
      )
    ],
    entries
  ],
  Print["CHESSV2_F3_MMA_MATRIX_SHAPE_FAILED=True"];
  Quit[3]
];
nodes = N[
  ChebyshevLobattoData[nodeCount, {0, 1}, workingPrecision][[1]],
  workingPrecision
];
zeroMatrices = ConstantArray[
  SparseArray[{}, {Length[boundaryValues], Length[boundaryValues]}],
  3
];
precomputedMatrices = Join[{zeroMatrices}, entries[[All, 3]]];
nodeSeconds = entries[[All, 2]];

ClearAll[F3PrecomputedMatrixList, F3PrecomputedMatrixListBatch];
F3PrecomputedMatrixList[tValue_, requestedPrecision_] :=
  precomputedMatrices[[
    First @ Ordering[
      Abs[nodes - SetPrecision[tValue, workingPrecision]], 1
    ]
  ]];
F3PrecomputedMatrixListBatch[
  requestedNodes_, requestedPrecision_, requestedThreads_
] := F3PrecomputedMatrixList[#, requestedPrecision] & /@ requestedNodes;

precomputedEvaluator = CHESSNodeEvaluator[
  F3PrecomputedMatrixList,
  F3PrecomputedMatrixListBatch
];
propagationSeconds = AbsoluteTiming[
  solution = SpectralPropagate[
    precomputedEvaluator,
    boundaryValues,
    {0, 1},
    "Nodes" -> nodeCount,
    "Precision" -> workingPrecision,
    "WorkingPrecisionA" -> workingPrecision,
    "EpsilonDegree" -> 2,
    "ParallelEvaluation" -> False,
    "ParallelKernels" -> 1
  ];
][[1]];
If[solution === $Failed,
  Print["CHESSV2_F3_MMA_PROPAGATION_FAILED=True"];
  Quit[4]
];

finalState = CHESSFinalState[solution];
perOrderErrors = Max /@ Transpose[Abs[finalState - referenceValues]];
maxError = Max[perOrderErrors];
HashHex[file_String] := IntegerString[FileHash[file, "SHA256"], 16, 64];
workerScript = FileNameJoin[{
  DirectoryName[ExpandFileName[$InputFileName]],
  "f3_zoia_mathematica_worker.wls"
}];
workerLoads = Table[
  Total @ Pick[
    nodeSeconds,
    Mod[indices - 2, 8] + 1,
    worker
  ],
  {worker, 1, 8}
];

record = <|
  "Example" -> "Simone Zoia F3 all-Mathematica node baseline",
  "Dimension" -> Length[boundaryValues],
  "Nodes" -> nodeCount,
  "WorkingPrecision" -> workingPrecision,
  "Environment" -> <|
    "WolframVersion" -> $Version,
    "SystemID" -> $SystemID,
    "ProcessorType" -> $ProcessorType,
    "MachineName" -> $MachineName,
    "CPU" -> "Intel Core i9-13950HX"
  |>,
  "ImportSeconds" -> importSeconds,
  "NodeEvaluationSecondsSum" -> Total[nodeSeconds],
  "WorkerEvaluationSeconds" -> workerLoads,
  "EstimatedParallelNodeSeconds" -> Max[workerLoads],
  "ObservedWorkerWallSeconds" -> observedWorkerWallSeconds,
  "WorkerTimeLog" -> workerTimeLog,
  "WorkerTimeLogSHA256" -> If[
    MissingQ[workerTimeLog], Missing["NotRecorded"], HashHex[workerTimeLog]
  ],
  "PropagationSeconds" -> propagationSeconds,
  "ObservedEndToEndSeconds" -> If[
    MissingQ[observedWorkerWallSeconds],
    Missing["NotRecorded"],
    observedWorkerWallSeconds + importSeconds + propagationSeconds
  ],
  "ProvenanceSHA256" -> <|
    "CHESSv2Loader" -> HashHex[packageFile],
    "CanonicalCore" -> HashHex[
      FileNameJoin[{packageRoot, "Core", "Canonical.wl"}]
    ],
    "NonCanonicalCore" -> HashHex[
      FileNameJoin[{packageRoot, "Core", "NonCanonical.wl"}]
    ],
    "Dispatch" -> HashHex[
      FileNameJoin[{packageRoot, "Core", "Dispatch.wl"}]
    ],
    "WorkerScript" -> HashHex[workerScript],
    "AssemblerScript" -> HashHex[ExpandFileName[$InputFileName]],
    "LegacyNonCanonicalPackage" -> HashHex[
      FileNameJoin[{projectRoot, "ChessNonCanonical.wl"}]
    ],
    "PreparedDAG" -> HashHex[
      FileNameJoin[{exampleRoot, "prepared_F3_physical_DAG_X1_X2.mx"}]
    ],
    "DAGEvaluator" -> HashHex[
      FileNameJoin[{exampleRoot, "scripts", "f3_physical_dag_evaluator.wl"}]
    ],
    "BoundaryX1" -> HashHex[boundaryFile],
    "ReferenceX2" -> HashHex[referenceFile]
  |>,
  "MaxAbsError" -> maxError,
  "PerOrderAbsError" -> perOrderErrors,
  "FinalState" -> finalState,
  "ReferenceValues" -> referenceValues
|>;
recordWrite = Quiet @ Check[
  Put[record, outputFile];
  outputFile,
  $Failed
];
If[recordWrite === $Failed || !FileExistsQ[outputFile],
  Print["CHESSV2_F3_MMA_RECORD_WRITE_FAILED=", outputFile];
  Quit[5]
];

Print["CHESSV2_F3_MMA_IMPORT_SECONDS=", N[importSeconds, 12]];
Print["CHESSV2_F3_MMA_WORKER_SECONDS=", N[workerLoads, 10]];
Print[
  "CHESSV2_F3_MMA_ESTIMATED_PARALLEL_NODE_SECONDS=",
  N[Max[workerLoads], 12]
];
Print[
  "CHESSV2_F3_MMA_OBSERVED_WORKER_WALL_SECONDS=",
  observedWorkerWallSeconds
];
Print["CHESSV2_F3_MMA_PROPAGATION_SECONDS=", N[propagationSeconds, 12]];
Print["CHESSV2_F3_MMA_MAX_ABS_ERROR=", N[maxError, 12]];
Print["CHESSV2_F3_MMA_OUTPUT=", outputFile];
Print["CHESSV2_F3_MMA_BENCHMARK_PASS=", TrueQ[maxError < 10^-30]];
If[TrueQ[maxError < 10^-30], Exit[0], Exit[5]];
