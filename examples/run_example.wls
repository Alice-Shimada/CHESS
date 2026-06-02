ClearAll["Global`*"];
$HistoryLength = 0;

commandLineArgs[] := Module[{scriptArgs, cmd, pos},
  scriptArgs = If[ListQ[$ScriptCommandLine] && Length[$ScriptCommandLine] > 0,
    Rest[$ScriptCommandLine],
    {}
  ];
  If[scriptArgs =!= {}, Return[scriptArgs]];
  cmd = If[ListQ[$CommandLine], $CommandLine, {}];
  pos = FirstPosition[cmd, $InputFileName, Missing["NotFound"], {1}];
  If[MatchQ[pos, {_Integer}], Drop[cmd, First[pos]], {}]
];

scriptFile = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  First[$ScriptCommandLine]
];
scriptDir = DirectoryName[ExpandFileName[scriptFile]];
rootDir = DirectoryName[scriptDir];
dataDir = FileNameJoin[{rootDir, "data"}];
packageFile = FileNameJoin[{rootDir, "Chess.wl"}];

args = commandLineArgs[];
caseName = If[Length[args] >= 1, ToUpperCase[args[[1]]], "DP"];
nodes = If[Length[args] >= 2, ToExpression[args[[2]]], 24];
kernels = If[Length[args] >= 3, ToExpression[args[[3]]], 1];
precision = If[Length[args] >= 4, ToExpression[args[[4]]], 80];
precisionA = If[Length[args] >= 5, ToExpression[args[[5]]], 100];

ClearAll[
  requireFile,
  dpCoefficientMatrix,
  maxAbs,
  formatNumber,
  formatValue,
  unifiedAtilde,
  referenceMatrix,
  setupDP,
  setup3L5P,
  setupBanana3,
  setupCase
];

requireFile[file_] := If[
  !FileExistsQ[file],
  Print["MISSING_FILE: ", file];
  Quit[2],
  file
];

dpCoefficientMatrix[file_, prec_Integer] := Module[{raw},
  raw = Get[requireFile[file]] /. InputForm[a__] :> a;
  SetPrecision[
    Table[SeriesCoefficient[expr, {ep, 0, j}], {expr, raw[[2]]}, {j, -4, 2}],
    prec
  ]
];

maxAbs[x_] := Max[Abs[Flatten[N[x, Min[60, precision]]]]];
formatNumber[x_] := ToString[FortranForm[N[x, 16]]];
formatValue[x_] := ToString[InputForm[N[x, 20]]];

unifiedAtilde[raw_] := SparseArray[
  raw /. {
    Log[W[i_Integer]] :> logW[i],
    Log[Wtilde[i_Integer]] :> logW[100 + i],
    Log[What[i_Integer]] :> logW[200 + i]
  }
];

referenceMatrix[iqbarVector_, qValue_, prec_Integer] := Module[{qPrec},
  qPrec = SetPrecision[qValue, prec];
  Block[
    {$MaxExtraPrecision = 1000},
    N[
      Table[
        SeriesCoefficient[iqbarVector[[row]], {eps, 0, coeff}] /. qbar -> qPrec,
        {row, Length[iqbarVector]},
        {coeff, 0, 6}
      ],
      prec
    ]
  ]
];

Get[requireFile[packageFile]];

setupDP[] := Module[{dir = FileNameJoin[{dataDir, "DP"}]},
  CHESSClearCaches[];
  dLettersLine = Get[requireFile[FileNameJoin[{dir, "dLettersLine.txt"}]]];
  Atilde = SparseArray[Get[requireFile[FileNameJoin[{dir, "Atilde.txt"}]]]];
  CHESSAtildeLinearData[];
  boundaryValues = dpCoefficientMatrix[FileNameJoin[{dir, "boundary_pr1_weight6.txt"}], precision];
  referenceValues = dpCoefficientMatrix[FileNameJoin[{dir, "reference_pr2_weight6.txt"}], precision];
  regularizedEndpoints = {};
  referenceLabel = "AMFlow physical-region point";
  atFunction = nAt;
  setupLabel = "DP physical-region transport";
];

setup3L5P[family_] := Module[{dir = FileNameJoin[{dataDir, "3L5P", family}], manifest},
  CHESSClearCaches[];
  dLettersLine = Get[requireFile[FileNameJoin[{dir, "dLettersLine.wl"}]]];
  Atilde = unifiedAtilde[Get[requireFile[FileNameJoin[{dir, "Atilde.m"}]]]];
  CHESSAtildeLinearData[];
  boundaryValues = SetPrecision[
    Transpose[Get[requireFile[FileNameJoin[{dir, "BCnum.m"}]]]],
    precision
  ];
  manifest = Get[requireFile[FileNameJoin[{dir, "manifest.m"}]]];
  referenceValues = Missing["NoIndependentReferenceInRelease"];
  regularizedEndpoints = {"Left"};
  referenceLabel = "node-convergence benchmark";
  atFunction = nAt;
  setupLabel = "3L5P " <> family <> " left-regularized transport";
  targetPoint = "TargetPoint" /. manifest;
];

setupBanana3[] := Module[
  {
    dir = FileNameJoin[{dataDir, "Banana3"}],
    config,
    q0,
    q1,
    matrixLine,
    iqbarVector
  },
  CHESSClearCaches[];
  config = Get[requireFile[FileNameJoin[{dir, "config_banana3.wl"}]]];
  Atilde = Get[requireFile[FileNameJoin[{dir, "Atilde_banana3.wl"}]]];
  matrixLine = Get[requireFile[FileNameJoin[{dir, "matrixLine_banana3.wl"}]]];
  iqbarVector = Get[requireFile[FileNameJoin[{dir, "IqbarVector_banana3.wl"}]]];
  q0 = "q0" /. config;
  q1 = "q1" /. config;
  ClearAll[banana3At];
  Options[banana3At] = {"Precision" -> precisionA};
  banana3At[tvalue_, OptionsPattern[]] := Module[{p},
    p = OptionValue["Precision"];
    SparseArray[N[matrixLine /. t -> N[tvalue, p], p]]
  ];
  boundaryValues = referenceMatrix[iqbarVector, q0, precision];
  referenceValues = referenceMatrix[iqbarVector, q1, precision];
  regularizedEndpoints = {};
  referenceLabel = "Iqbar q-series reference at q1";
  atFunction = banana3At;
  setupLabel = "three-loop equal-mass banana transport";
  Quiet @ Check[DistributeDefinitions[matrixLine, banana3At, t], Null];
];

setupCase[name_] := Which[
  MemberQ[{"DP", "DPPHYSICAL", "DP_PHYSICAL"}, name],
    setupDP[],
  MemberQ[{"PBB", "BPB", "BHB", "PBP"}, name],
    setup3L5P[name],
  MemberQ[{"BANANA", "BANANA3", "ELLIPTICBANANA"}, name],
    setupBanana3[],
  True,
    Print["UNKNOWN_CASE: ", name, ". Use DP, PBB, BPB, BHB, PBP, or BANANA3."];
    Quit[2]
];

Print[
  "CHESS_RELEASE_EXAMPLE case=", caseName,
  " nodes=", nodes,
  " kernels=", kernels,
  " precision=", precision,
  " precisionA=", precisionA
];

setupCase[caseName];

Print[
  "setupLabel=", setupLabel,
  " matrixDimensions=", Dimensions[Atilde],
  " boundaryDimensions=", Dimensions[boundaryValues],
  " regularizedEndpoints=", regularizedEndpoints
];

{wallTime, result} = AbsoluteTiming[
  SpectralPropagate[
    atFunction,
    boundaryValues,
    {0, 1},
    "Nodes" -> nodes,
    "Precision" -> precision,
    "WorkingPrecisionA" -> precisionA,
    "EndpointRegularizationPrecision" -> precisionA,
    "ParallelEvaluation" -> If[kernels > 1, True, False],
    "ParallelKernels" -> kernels,
    "RegularizedEndpoints" -> regularizedEndpoints,
    "SolveStrategy" -> "SequentialEpsilon"
  ]
];

If[result === $Failed,
  Print["CHESS_RELEASE_EXAMPLE_FAILED"];
  Quit[2];
];

final = CHESSFinalState[result];
referenceError = If[ListQ[referenceValues], maxAbs[final - referenceValues], Missing["NotAvailable"]];

Print[
  "resultDimensions=", Dimensions[final],
  " wallTimeSeconds=", formatNumber[wallTime],
  " referenceLabel=", referenceLabel,
  " maxReferenceError=", If[NumberQ[referenceError], formatNumber[referenceError], "not_available"]
];

Print[
  "endpointInfo=", CHESSEndpointInfo[result],
  " finalChecksum=", formatValue[Total[Flatten[final]]]
];

Quiet[CloseKernels[]];
Quit[0];
