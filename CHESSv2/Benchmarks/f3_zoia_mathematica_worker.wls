#!/usr/bin/env wolframscript

(* One isolated worker for the historical all-Mathematica F3 node evaluator.
   The output directory is supplied explicitly so a fresh benchmark never
   overwrites the preserved CHESS_non_canonical matrices. *)

ClearAll["Global`*"];
$HistoryLength = 0;

If[Length[$ScriptCommandLine] =!= 7,
  Print[
    "usage: f3_zoia_mathematica_worker.wls EXAMPLE_ROOT NODES PRECISION ",
    "WORKERS WORKER_INDEX OUTPUT_DIRECTORY"
  ];
  Quit[2]
];

exampleRoot = ExpandFileName[$ScriptCommandLine[[2]]];
nodeCount = ToExpression[$ScriptCommandLine[[3]]];
workingPrecision = ToExpression[$ScriptCommandLine[[4]]];
workerCount = ToExpression[$ScriptCommandLine[[5]]];
workerIndex = ToExpression[$ScriptCommandLine[[6]]];
outputDirectory = ExpandFileName[$ScriptCommandLine[[7]]];
projectRoot = DirectoryName[exampleRoot];
packageFile = FileNameJoin[{projectRoot, "ChessNonCanonical.wl"}];
cacheFile = FileNameJoin[{
  exampleRoot, "prepared_F3_physical_DAG_X1_X2.mx"
}];
evaluatorFile = FileNameJoin[{
  exampleRoot, "scripts", "f3_physical_dag_evaluator.wl"
}];

If[
  !And @@ (FileExistsQ /@ {packageFile, cacheFile, evaluatorFile}) ||
  !IntegerQ[nodeCount] || nodeCount <= 0 ||
  !IntegerQ[workerCount] || workerCount <= 0 ||
  !IntegerQ[workerIndex] || !Between[workerIndex, {1, workerCount}] ||
  !DirectoryQ[outputDirectory],
  Print["CHESSV2_F3_MMA_WORKER_INPUT_FAILED=True"];
  Quit[2]
];

loadResult = Quiet @ Check[
  Get[packageFile];
  Get[cacheFile];
  Get[evaluatorFile];
  True,
  $Failed
];
If[loadResult === $Failed,
  Print["CHESSV2_F3_MMA_WORKER_LOAD_FAILED=True"];
  Quit[2]
];

nodes = N[
  ChebyshevLobattoData[nodeCount, {0, 1}, workingPrecision][[1]],
  workingPrecision
];
indices = Select[
  Range[2, nodeCount + 1],
  Mod[# - 2, workerCount] + 1 === workerIndex &
];

Do[
  nodeSeconds = AbsoluteTiming[
    nodeMatrices = F3PhysicalDAGMatrixList[
      nodes[[nodeIndex]], workingPrecision
    ];
  ][[1]];
  If[
    !ListQ[nodeMatrices] || Length[nodeMatrices] =!= 3 ||
    !And @@ (
      MatrixQ[#, NumericQ] && Dimensions[#] === {131, 131} & /@
        nodeMatrices
    ),
    Print[
      "CHESSV2_F3_MMA_WORKER_MATRIX_FAILED=", workerIndex,
      " NODE=", nodeIndex
    ];
    Quit[3]
  ];
  matrixFile = FileNameJoin[{
    outputDirectory,
    "matrix_" <> IntegerString[nodeIndex, 10, 4] <> ".m"
  }];
  If[FileExistsQ[matrixFile],
    Print["CHESSV2_F3_MMA_WORKER_REFUSE_OVERWRITE=", matrixFile];
    Quit[3]
  ];
  writeResult = Quiet @ Check[
    Put[{nodeIndex, nodeSeconds, nodeMatrices}, matrixFile];
    matrixFile,
    $Failed
  ];
  If[writeResult === $Failed || !FileExistsQ[matrixFile],
    Print["CHESSV2_F3_MMA_WORKER_WRITE_FAILED=", matrixFile];
    Quit[4]
  ];
  Print[
    "CHESSV2_F3_MMA_WORKER=", workerIndex,
    " NODE=", nodeIndex,
    " SECONDS=", N[nodeSeconds, 8]
  ];
  Clear[nodeMatrices];
  ClearSystemCache[],
  {nodeIndex, indices}
];

Print["CHESSV2_F3_MMA_WORKER_PASS=", workerIndex];
