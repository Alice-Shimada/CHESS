ClearAll["Global`*"];
$HistoryLength = 0;

scriptFile = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  First[$ScriptCommandLine]
];
scriptDir = DirectoryName[ExpandFileName[scriptFile]];
rootDir = DirectoryName[scriptDir];
dataDir = FileNameJoin[{rootDir, "data"}];
packageFile = FileNameJoin[{rootDir, "Chess.wl"}];

args = Rest[$ScriptCommandLine];
caseName = If[Length[args] >= 1, ToUpperCase[args[[1]]], "DP"];
nodes = If[Length[args] >= 2, ToExpression[args[[2]]], 24];
kernels = If[Length[args] >= 3, ToExpression[args[[3]]], 1];
precision = If[Length[args] >= 4, ToExpression[args[[4]]], 80];
precisionA = If[Length[args] >= 5, ToExpression[args[[5]]], 100];

ClearAll[
  requireFile,
  boundaryCoefficients,
  expressionCoefficients,
  maxAbs,
  formatNumber,
  formatValue,
  unifiedAtilde,
  setupCase
];

requireFile[file_] := If[
  !FileExistsQ[file],
  Print["MISSING_FILE: ", file];
  Quit[2],
  file
];

boundaryCoefficients[file_, emin_Integer, emax_Integer] :=
  Table[SeriesCoefficient[#, {ep, 0, j}], {j, emin, emax}] & /@
    (Get[file] /. InputForm[a__] :> a)[[2]];

expressionCoefficients[file_, emin_Integer, emax_Integer] :=
  Table[SeriesCoefficient[#, {eps, 0, j}], {j, emin, emax}] & /@ Get[file];

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

Get[requireFile[packageFile]];

setupCase["DP"] := Module[{dir = FileNameJoin[{dataDir, "DP"}], emin = -4, emax = 2},
  CHESSClearCaches[];
  dLettersLine = Get[requireFile[FileNameJoin[{dir, "dLettersLine.txt"}]]];
  Atilde = Get[requireFile[FileNameJoin[{dir, "Atilde.txt"}]]];
  CHESSAtildeLinearData[];
  boundaryValues = boundaryCoefficients[
    FileNameJoin[{dir, "boundary_pr1_weight6.txt"}],
    emin,
    emax
  ];
  referenceValues = boundaryCoefficients[
    FileNameJoin[{dir, "reference_pr2_weight6.txt"}],
    emin,
    emax
  ];
  regularizedEndpoints = {};
  referenceLabel = "DP reference point";
];

setupCase["HB"] := Module[{dir = FileNameJoin[{dataDir, "HB"}], emin = -4, emax = 0},
  CHESSClearCaches[];
  dLettersLine = Get[requireFile[FileNameJoin[{dir, "dLettersLine.m"}]]];
  Atilde = Get[requireFile[FileNameJoin[{dir, "Atilde.txt"}]]];
  CHESSAtildeLinearData[];
  boundaryValues = boundaryCoefficients[
    FileNameJoin[{dir, "boundary_E1_weight6.txt"}],
    emin,
    emax
  ];
  referenceValues = expressionCoefficients[
    FileNameJoin[{dir, "reference_E1_60digits_noExpGamma.m"}],
    emin,
    emax
  ];
  regularizedEndpoints = {"Right"};
  referenceLabel = "HB independent reference";
];

setupCase["BHB"] := Module[{dir = FileNameJoin[{dataDir, "BHB"}]},
  CHESSClearCaches[];
  dLettersLine = Get[requireFile[FileNameJoin[{dir, "dLettersLine.wl"}]]];
  Atilde = unifiedAtilde[Get[requireFile[FileNameJoin[{dir, "Atilde.m"}]]]];
  CHESSAtildeLinearData[];
  boundaryValues = SetPrecision[
    Transpose[Get[requireFile[FileNameJoin[{dir, "boundary_num.m"}]]]],
    precision
  ];
  referenceValues = Missing["NoReferenceInRelease"];
  regularizedEndpoints = {"Left"};
  referenceLabel = "node-convergence example";
];

setupCase[_] := (
  Print["UNKNOWN_CASE: ", caseName, ". Use DP, HB, or BHB."];
  Quit[2];
);

Print[
  "CHESS_RELEASE_EXAMPLE case=", caseName,
  " nodes=", nodes,
  " kernels=", kernels,
  " precision=", precision,
  " precisionA=", precisionA
];

setupCase[caseName];

Print[
  "matrixDimensions=", Dimensions[Atilde],
  " boundaryDimensions=", Dimensions[boundaryValues],
  " regularizedEndpoints=", regularizedEndpoints
];

{wallTime, result} = AbsoluteTiming[
  SpectralPropagate[
    nAt,
    boundaryValues,
    {0, 1},
    "Nodes" -> nodes,
    "Precision" -> precision,
    "WorkingPrecisionA" -> precisionA,
    "EndpointRegularizationPrecision" -> precisionA,
    "ParallelEvaluation" -> If[kernels > 1, True, False],
    "ParallelKernels" -> kernels,
    "RegularizedEndpoints" -> regularizedEndpoints
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

Quit[0];
