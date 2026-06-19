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
  valueOf,
  unifiedAtilde,
  setupDP,
  setup3L5P,
  setupBhabha,
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
valueOf[rules_, key_] := key /. rules;

unifiedAtilde[raw_] := SparseArray[
  raw /. {
    Log[W[i_Integer]] :> logW[i],
    Log[Wtilde[i_Integer]] :> logW[100 + i],
    Log[What[i_Integer]] :> logW[200 + i]
  }
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
  transportInterval = {0, 1};
  systemDimensions = Dimensions[Atilde];
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
  transportInterval = {0, 1};
  systemDimensions = Dimensions[Atilde];
];

setupBhabha[] := Module[
  {
    dir = FileNameJoin[{dataDir, "BhabhaPL2"}],
    matrixDir,
    refDir,
    leftRules,
    rightRules,
    dsMatrix,
    dtMatrix,
    sampleMatrix,
    xStart,
    xEnd
  },
  CHESSClearCaches[];
  matrixDir = FileNameJoin[{dir, "matrices"}];
  refDir = FileNameJoin[{dir, "reference"}];
  leftRules = Get[requireFile[FileNameJoin[{refDir, "boundary_x0p2.wl"}]]];
  rightRules = Get[requireFile[FileNameJoin[{refDir, "reference_x1.wl"}]]];
  If[valueOf[leftRules, "Status"] =!= "ok" || valueOf[rightRules, "Status"] =!= "ok",
    Print["BAD_BHABHA_REFERENCE_STATUS"];
    Quit[2];
  ];
  dsMatrix = Get[requireFile[FileNameJoin[{matrixDir, "ds_1.m"}]]];
  dtMatrix = Get[requireFile[FileNameJoin[{matrixDir, "dt_1.m"}]]];
  xStart = SetPrecision[valueOf[leftRules, "x"], precisionA];
  xEnd = SetPrecision[valueOf[rightRules, "x"], precisionA];
  ClearAll[bhabhaAt];
  Options[bhabhaAt] = {"Precision" -> precisionA};
  bhabhaAt[xvalue_, OptionsPattern[]] := Module[{p, xx, rules},
    p = OptionValue["Precision"];
    xx = N[xvalue, p];
    rules = {s -> N[3 xx/10, p], t -> N[-7 xx/10, p], m2 -> 1};
    SparseArray @ N[
      (3/10) (dsMatrix /. rules) - (7/10) (dtMatrix /. rules),
      p
    ]
  ];
  boundaryValues = N[valueOf[leftRules, "Coefficients"], precision];
  referenceValues = N[valueOf[rightRules, "Coefficients"], precision];
  regularizedEndpoints = {};
  referenceLabel = "AMFlow endpoint at x=1";
  atFunction = bhabhaAt;
  setupLabel = "massive Bhabha planar direct-matrix transport";
  transportInterval = {xStart, xEnd};
  sampleMatrix = atFunction[(xStart + xEnd)/2, "Precision" -> precisionA];
  systemDimensions = Dimensions[sampleMatrix];
  If[systemDimensions =!= {43, 43} ||
      !FreeQ[sampleMatrix, Indeterminate | ComplexInfinity],
    Print["BAD_BHABHA_SAMPLE_MATRIX dimensions=", systemDimensions];
    Quit[2];
  ];
  Quiet @ Check[DistributeDefinitions[dsMatrix, dtMatrix, bhabhaAt, s, t, m2], Null];
];

setupCase[name_] := Which[
  MemberQ[{"DP", "DPPHYSICAL", "DP_PHYSICAL"}, name],
    setupDP[],
  MemberQ[{"PBB", "BPB", "BHB", "PBP"}, name],
    setup3L5P[name],
  MemberQ[{"BHABHA", "BHABHAPL2", "BHABHA_PL2"}, name],
    setupBhabha[],
  True,
    Print["UNKNOWN_CASE: ", name, ". Use DP, PBB, BPB, BHB, PBP, or BHABHA."];
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
  " matrixDimensions=", systemDimensions,
  " boundaryDimensions=", Dimensions[boundaryValues],
  " regularizedEndpoints=", regularizedEndpoints,
  " interval=", transportInterval
];

{wallTime, result} = AbsoluteTiming[
  SpectralPropagate[
    atFunction,
    boundaryValues,
    transportInterval,
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
