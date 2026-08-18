#!/usr/bin/env wolframscript

(*
  Reproduce the large symbolic F3 differential-equation benchmark supplied by
  Simone Zoia.  The 639 expanded polynomials are evaluated at all 96 Lobatto
  nodes by the standalone FORM/FLINT batch module.  The example-specific path,
  roots, one-forms, and sparse matrix assembly are loaded from the preserved
  CHESS_non_canonical experiment; spectral transport goes through CHESSv2's
  unified SpectralPropagate entry.

  Timing boundaries are intentionally explicit:

    FLINT process       memory-map and evaluate the prepared 790 MiB SLP;
    FLINT interface     point serialization + process + binary decoding;
    matrix assembly     compact roots/one-forms + sparse {B0,B1,B2};
    node backend        path coordinates + FLINT interface + assembly;
    propagation         node backend + unchanged non-canonical solve.

  The one-time FORM optimisation and prepared-cache construction are excluded,
  just as in the historical 462.1 s versus 38.8--39.0 s comparison.
*)

ClearAll["Global`*"];
$HistoryLength = 0;
$MaxExtraPrecision = 5000;

packageRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
packageFile = FileNameJoin[{packageRoot, "CHESSv2.wl"}];
exampleRoot = Replace[
  Environment["CHESSV2_F3_EXAMPLE_ROOT"],
  value_ /; !StringQ[value] || value === "" :>
    "/home/liuyuanche/Program/Chebyshev_nonUT/ttW_physical_variables_2608_03746"
];
projectRoot = DirectoryName[exampleRoot];
externalBackendFile = FileNameJoin[{
  exampleRoot, "backend", "F3PhysicalBackend.wl"
}];
evaluatorSourceFile = FileNameJoin[{
  exampleRoot, "backend", "chess_slp_evaluator.cpp"
}];
flintLibraryFile = "/home/liuyuanche/MMA_Package/flint-3.4.0/libflint.so.22";
boundaryFile = FileNameJoin[{
  projectRoot, "ttW_F3", "ancillary_files", "F3", "values", "F3_X1.m"
}];
referenceFile = FileNameJoin[{
  projectRoot, "ttW_F3", "ancillary_files", "F3", "values", "F3_X2.m"
}];
outputFile = Replace[
  Environment["CHESSV2_F3_BENCHMARK_OUTPUT"],
  value_ /; !StringQ[value] || value === "" :>
    FileNameJoin[{$TemporaryDirectory, "chessv2-f3-zoia-flint.wl"}]
];

requiredFiles = {
  packageFile,
  externalBackendFile,
  evaluatorSourceFile,
  flintLibraryFile,
  FileNameJoin[{exampleRoot, "backend", "chess_slp_evaluator"}],
  FileNameJoin[{exampleRoot, "backend", "F3_polynomials.chslp"}],
  FileNameJoin[{exampleRoot, "backend", "prepared_F3_backend.mx"}],
  boundaryFile,
  referenceFile
};
If[!And @@ (FileExistsQ /@ requiredFiles),
  Print["CHESSV2_F3_REQUIRED_FILE_MISSING=", Pick[
    requiredFiles, FileExistsQ /@ requiredFiles, False
  ]];
  Quit[2]
];
If[FileExistsQ[outputFile],
  Print["CHESSV2_F3_REFUSE_OVERWRITE=", outputFile];
  Quit[2]
];

nodeCount = 96;
workingPrecision = 100;
threadCount = 8;

(* Loading F3PhysicalBackend.wl supplies only the example-specific prepared
   path and matrix-assembly functions.  Its duplicated process/decoder path is
   deliberately not called; CHESSFlintRunSLP below owns the FLINT interface. *)
{setupSeconds, setupResult} = AbsoluteTiming[
  Quiet @ Check[
    Get[packageFile];
    Get[externalBackendFile];
    initialValues = Last[Get[boundaryFile]];
    finalValues = Last[Get[referenceFile]];
    True,
    $Failed
  ]
];
If[setupResult === $Failed,
  Print["CHESSV2_F3_SETUP_FAILED=True"];
  Quit[2]
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

ClearAll[
  CHESSV2F3MatrixList,
  CHESSV2F3MatrixListBatch,
  CHESSV2F3LastTiming,
  $CHESSV2F3LastTiming
];
$CHESSV2F3LastTiming = <||>;

(* This is the thin problem adapter: seven physical path coordinates go into
   the generic FLINT batch protocol; the returned 639 values feed the already
   prepared sparse F3 matrix assembly. *)
CHESSV2F3MatrixListBatch[
  nodes_List, requestedPrecision_, requestedThreads_
] := Module[
  {
    totalStart, digits, threads, pathRows, pathSeconds, flintResult,
    flintInterfaceSeconds, tasks, matrices, assemblySeconds, totalSeconds,
    matricesValidQ
  },
  If[nodes === {}, Return[{}]];
  totalStart = AbsoluteTime[];
  digits = Max[20, Ceiling[N[requestedPrecision]]];
  threads = Max[1, Min[
    Length[nodes],
    Replace[requestedThreads, {value_Integer?Positive :> value, _ :> 1}]
  ]];

  pathSeconds = AbsoluteTiming[
    pathRows = Table[
      N[#[SetPrecision[node, digits]], digits] & /@ F3PhysicalPathFunctions,
      {node, nodes}
    ];
  ][[1]];

  flintInterfaceSeconds = AbsoluteTiming[
    flintResult = CHESSFlintRunSLP[
      $F3PhysicalBackendExecutable,
      $F3PhysicalBackendSLP,
      pathRows,
      digits,
      threads
    ];
  ][[1]];
  If[
    !AssociationQ[flintResult] || flintResult["ValueCount"] =!= 639,
    Return[$Failed]
  ];

  tasks = MapThread[
    {#1, #2, #3, digits} &,
    {nodes, pathRows, flintResult["Values"]}
  ];
  assemblySeconds = AbsoluteTiming[
    matrices = If[threads > 1 && Length[nodes] >= 8,
      F3PhysicalBackendPrepareParallel[threads];
      ParallelMap[
        F3PhysicalBackendAssembleTask,
        tasks,
        Method -> "CoarsestGrained"
      ],
      F3PhysicalBackendAssembleTask /@ tasks
    ];
  ][[1]];
  matricesValidQ = ListQ[matrices] && Length[matrices] === Length[nodes] &&
    And @@ Map[
      Function[nodeMatrices,
        ListQ[nodeMatrices] && Length[nodeMatrices] === 3 &&
        And @@ (
          MatrixQ[#, NumericQ] &&
          Dimensions[#] === ConstantArray[Length[boundaryValues], 2] & /@
            nodeMatrices
        )
      ],
      matrices
    ];
  If[!TrueQ[matricesValidQ], Return[$Failed]];

  totalSeconds = N[AbsoluteTime[] - totalStart];
  $CHESSV2F3LastTiming = <|
    "PathSeconds" -> pathSeconds,
    "FLINTInterfaceSeconds" -> flintInterfaceSeconds,
    "FLINTProcessSeconds" -> flintResult["FLINTProcessSeconds"],
    "AssemblySeconds" -> assemblySeconds,
    "TotalNodeBackendSeconds" -> totalSeconds,
    "Bits" -> flintResult["Bits"],
    "Threads" -> flintResult["Threads"],
    "NodeCount" -> Length[nodes]
  |>;
  matrices
];

CHESSV2F3MatrixList[tValue_, requestedPrecision_?NumericQ] := First[
  CHESSV2F3MatrixListBatch[{tValue}, requestedPrecision, 1]
];
CHESSV2F3LastTiming[] := $CHESSV2F3LastTiming;

f3Evaluator = CHESSNodeEvaluator[
  CHESSV2F3MatrixList,
  CHESSV2F3MatrixListBatch
];

propagationSeconds = AbsoluteTiming[
  solution = SpectralPropagate[
    f3Evaluator,
    boundaryValues,
    {0, 1},
    "Nodes" -> nodeCount,
    "Precision" -> workingPrecision,
    "WorkingPrecisionA" -> workingPrecision,
    "EpsilonDegree" -> 2,
    "ParallelEvaluation" -> True,
    "ParallelKernels" -> threadCount
  ];
][[1]];
If[solution === $Failed,
  Print["CHESSV2_F3_PROPAGATION_FAILED=True"];
  Quit[3]
];

nodeTiming = CHESSV2F3LastTiming[];
finalState = CHESSFinalState[solution];
perOrderErrors = Max /@ Transpose[Abs[finalState - referenceValues]];
maxError = Max[perOrderErrors];
HashHex[file_String] := IntegerString[FileHash[file, "SHA256"], 16, 64];

record = <|
  "Example" -> "Simone Zoia F3 X1 to X2 physical-coordinate path",
  "Dimension" -> Length[boundaryValues],
  "Nodes" -> nodeCount,
  "WorkingPrecision" -> workingPrecision,
  "Threads" -> threadCount,
  "EpsilonMatrixDegree" -> 2,
  "SolutionOrders" -> Range[0, 4],
  "Environment" -> <|
    "WolframVersion" -> $Version,
    "SystemID" -> $SystemID,
    "ProcessorType" -> $ProcessorType,
    "MachineName" -> $MachineName,
    "CPU" -> "Intel Core i9-13950HX",
    "Compiler" -> "g++ 11.4.0",
    "FLINT" -> "3.4.0"
  |>,
  "SetupSeconds" -> setupSeconds,
  "NodeTiming" -> nodeTiming,
  "PropagationSeconds" -> propagationSeconds,
  "PostNodeSolveSeconds" ->
    propagationSeconds - nodeTiming["TotalNodeBackendSeconds"],
  "HistoricalMathematicaNodeSeconds" -> 462.1,
  "HistoricalFLINTNodeSecondsRange" -> {38.769, 39.010},
  "NodeSpeedupVersusHistoricalMathematica" ->
    462.1/nodeTiming["TotalNodeBackendSeconds"],
  "MaxAbsError" -> maxError,
  "PerOrderAbsError" -> perOrderErrors,
  "BoundarySource" -> boundaryFile,
  "ReferenceSource" -> referenceFile,
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
    "FLINTBatchEvaluation" -> HashHex[
      FileNameJoin[{packageRoot, "Core", "FLINTBatchEvaluation.wl"}]
    ],
    "BenchmarkDriver" -> HashHex[ExpandFileName[$InputFileName]],
    "F3AssemblyAdapter" -> HashHex[externalBackendFile],
    "EvaluatorSource" -> HashHex[evaluatorSourceFile],
    "FLINTLibrary" -> HashHex[flintLibraryFile],
    "BoundaryX1" -> HashHex[boundaryFile],
    "ReferenceX2" -> HashHex[referenceFile]
  |>,
  "EvaluatorExecutableSHA256" ->
    HashHex[$F3PhysicalBackendExecutable],
  "PreparedSLPSHA256" ->
    HashHex[$F3PhysicalBackendSLP],
  "PreparedCacheSHA256" ->
    HashHex[$F3PhysicalBackendCache],
  "FinalState" -> finalState,
  "ReferenceValues" -> referenceValues
|>;
recordWrite = Quiet @ Check[
  Put[record, outputFile];
  outputFile,
  $Failed
];
If[recordWrite === $Failed || !FileExistsQ[outputFile],
  Print["CHESSV2_F3_RECORD_WRITE_FAILED=", outputFile];
  Quit[5]
];

Print["CHESSV2_F3_NODE_TIMING=", InputForm[nodeTiming]];
Print["CHESSV2_F3_PROPAGATION_SECONDS=", N[propagationSeconds, 12]];
Print[
  "CHESSV2_F3_NODE_SPEEDUP_VS_HISTORICAL_MMA=",
  N[record["NodeSpeedupVersusHistoricalMathematica"], 12]
];
Print["CHESSV2_F3_MAX_ABS_ERROR=", N[maxError, 12]];
Print["CHESSV2_F3_PER_ORDER_ERRORS=", N[perOrderErrors, 12]];
Print["CHESSV2_F3_OUTPUT=", outputFile];
Print["CHESSV2_F3_BENCHMARK_PASS=", TrueQ[maxError < 10^-30]];
If[TrueQ[maxError < 10^-30], Exit[0], Exit[4]];
