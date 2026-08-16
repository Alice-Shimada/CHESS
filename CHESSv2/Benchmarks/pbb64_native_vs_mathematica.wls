#!/usr/bin/env wolframscript

(*
  Cached native versus all-Mathematica CHESSv2 propagation on the public
  316-dimensional PBB differential equation.

  Only SpectralPropagate calls are timed.  Data loading, path construction,
  Atilde decomposition, C++ process installation, node-matrix sampling, setup
  serialization, and native LU factorization are outside those timers.  Native
  setup time is nevertheless printed separately so both the repeated-call and
  one-shot tradeoffs remain visible.

  Both routes solve exactly the same B0-only equation on the same regular path,
  with 64 decimal digits, 81 Lobatto intervals (82 points), and fixed auxiliary
  delta order 56.  Endpoint mode prevents native node-state serialization.
*)

ClearAll["Global`*"];
$HistoryLength = 0;
$MaxExtraPrecision = 5000;

packageRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
packageFile = FileNameJoin[{packageRoot, "CHESSv2.wl"}];
dataRoot = Replace[
  Environment["PBB64_DATA_ROOT"],
  value_ /; !StringQ[value] || value === "" :>
    "/tmp/All_Planar_3L5P-3e5862a2"
];
outputFile = Replace[
  Environment["CHESSV2_BENCHMARK_OUTPUT"],
  value_ /; !StringQ[value] || value === "" :>
    FileNameJoin[{$TemporaryDirectory, "chessv2-pbb64-native-vs-mma.wl"}]
];
repeats = Replace[
  Environment["CHESSV2_BENCHMARK_REPEATS"],
  value_String /; StringMatchQ[value, DigitCharacter ..] :> FromDigits[value],
  {0}
];
If[!IntegerQ[repeats] || repeats < 3, repeats = 5];
ompThreads = Replace[
  Environment["OMP_NUM_THREADS"],
  value_String /; StringMatchQ[value, DigitCharacter ..] :> FromDigits[value],
  {0}
];
If[!IntegerQ[ompThreads] || ompThreads <= 0, ompThreads = 8];
SetEnvironment["OMP_NUM_THREADS" -> ToString[ompThreads]];

requiredFiles = {
  packageFile,
  FileNameJoin[{dataRoot, "Letters.m"}],
  FileNameJoin[{dataRoot, "PBB", "Atilde.m"}],
  FileNameJoin[{dataRoot, "PBB", "BCnum.m"}]
};
If[!And @@ (FileExistsQ /@ requiredFiles),
  Print["PBB64_REQUIRED_FILE_MISSING=", Pick[
    requiredFiles, FileExistsQ /@ requiredFiles, False
  ]];
  Quit[2]
];

workingPrecision = 64;
matrixPrecision = 90;
nodes = 81;
deltaOrder = 56;
physicalEpsilon = 2275901/2275800;
pathStart = -17/10;
pathEnd = -19/10;

(* Rational pole-free path used in the earlier AMFlow/CHESS comparison. *)
s12OnPath = (48 + 36 Global`t)/(16 - Global`t^2);
eps5OnPath = 18 + Global`t s12OnPath;

(* Load the source alphabet before CHESS so W and the Atilde symbols resolve in
   the same context. *)
Get[FileNameJoin[{dataRoot, "Letters.m"}]];
sourceLetterRep = ToExpression["LetterRep"];
sourceRootDef = ToExpression["RootDef"];
sourceW = ToExpression["W"];
sourceKinematicSymbols = ToExpression /@ {"s23", "s34", "s45", "s15", "s12"};
Get[packageFile];

Atilde = Get[FileNameJoin[{dataRoot, "PBB", "Atilde.m"}]];
boundarySeries = Get[FileNameJoin[{dataRoot, "PBB", "BCnum.m"}]];
If[Dimensions[Atilde] =!= {316, 316} || Dimensions[boundarySeries] =!= {7, 316},
  Print["PBB64_INPUT_DIMENSION_FAILED=True"];
  Quit[2]
];

CHESSClearCaches[];
linearData = CHESSAtildeLinearData[];
If[linearData === $Failed, Print["PBB64_LINEAR_DATA_FAILED=True"]; Quit[2]];
usedSlots = linearData[[2]];
letterRules = Select[
  sourceLetterRep,
  Head[First[#]] === sourceW && MemberQ[usedSlots, First[#][[1]]] &
];
pathRules = Join[
  Thread[sourceKinematicSymbols -> {-2, -3, -4, -6, s12OnPath}],
  {First[sourceRootDef] -> eps5OnPath}
];
letterValues = Table[
  Module[{value = sourceW[index] /. letterRules},
    If[Head[value] === sourceW, 1, Together[value /. pathRules]]
  ],
  {index, Max[usedSlots]}
];
dLettersLine = physicalEpsilon *
  (Together[D[Log[#], Global`t]] & /@ letterValues);

boundaryVector = N[
  Sum[
    physicalEpsilon^(order - 1) boundarySeries[[order]],
    {order, Length[boundarySeries]}
  ],
  workingPrecision
];
boundary = Transpose[{boundaryVector}];
interval = N[{pathStart, pathEnd}, workingPrecision];
commonRules = {
  "Nodes" -> nodes,
  "Precision" -> workingPrecision,
  "WorkingPrecisionA" -> matrixPrecision,
  "ParallelEvaluation" -> False,
  "ParallelKernels" -> 1,
  "DeltaOrder" -> deltaOrder,
  "MaxDeltaOrder" -> deltaOrder
};

(* Preparation is explicitly outside propagation timers and is reusable for
   every boundary matrix on the same equation/path/options. *)
{prepareSeconds, nativeHandle} = AbsoluteTiming[
  CHESSNativePrepare[nAt, 316, interval, commonRules]
];
If[!AssociationQ[nativeHandle], Print["PBB64_NATIVE_PREPARE_FAILED=True"]; Quit[3]];

runMathematica[] := Module[{result},
  result = SpectralPropagate[
    {nAt}, boundary, interval,
    Sequence @@ commonRules,
    "NumericalBackend" -> "Mathematica"
  ];
  If[result === $Failed, Print["PBB64_MATHEMATICA_FAILED=True"]; Quit[4]];
  result
];

runNative[] := Module[{result},
  result = SpectralPropagate[
    {nAt}, boundary, interval,
    Sequence @@ commonRules,
    "NumericalBackend" -> "Native",
    "NativeHandle" -> nativeHandle,
    "ResultData" -> "Endpoint"
  ];
  If[result === $Failed, Print["PBB64_NATIVE_FAILED=True"]; Quit[5]];
  result
];

(* Timed cold calls warm both implementations.  Each paired repeat
   uses the fixed Mathematica-then-native order; the raw vectors are retained
   so a caller can assess timing drift instead of hiding it in one aggregate. *)
{mathematicaColdSeconds, mathematicaResult} = AbsoluteTiming[runMathematica[]];
{nativeColdSeconds, nativeResult} = AbsoluteTiming[runNative[]];
paired = Table[
  {
    AbsoluteTiming[mathematicaResult = runMathematica[]],
    AbsoluteTiming[nativeResult = runNative[]]
  },
  {repeats}
];
mathematicaTimes = paired[[All, 1, 1]];
nativeTimes = paired[[All, 2, 1]];
mathematicaMedian = Median[mathematicaTimes];
nativeMedian = Median[nativeTimes];
nativeMetadata = Association[
  Rule @@@ Partition[Drop[nativeResult[[6]], 3], 2]
];
kernelSeconds = Lookup[nativeMetadata, "KernelSeconds", Missing["Unavailable"]];
scaledError = Max[
  Abs[Flatten[nativeResult[[3]] - mathematicaResult[[3]]]]/
  (1 + Abs[Flatten[mathematicaResult[[3]]]])
];

record = <|
  "DataCommit" -> "3e5862a2c28bd1bd080d141894db85a010bce2c1",
  "Dimension" -> 316,
  "WorkingPrecision" -> workingPrecision,
  "MatrixPrecision" -> matrixPrecision,
  "Nodes" -> nodes,
  "DeltaOrder" -> deltaOrder,
  "PhysicalEpsilon" -> physicalEpsilon,
  "Path" -> {pathStart, pathEnd},
  "OMPThreads" -> ompThreads,
  "Repeats" -> repeats,
  "NativePrepareSeconds" -> prepareSeconds,
  "NativeLoadSeconds" -> nativeHandle["LoadSeconds"],
  "MathematicaColdSeconds" -> mathematicaColdSeconds,
  "NativeColdSeconds" -> nativeColdSeconds,
  "MathematicaTimesSeconds" -> mathematicaTimes,
  "NativeTimesSeconds" -> nativeTimes,
  "MathematicaMedianSeconds" -> mathematicaMedian,
  "NativeMedianSeconds" -> nativeMedian,
  "NativeLastKernelSeconds" -> kernelSeconds,
  "CachedSpeedup" -> mathematicaMedian/nativeMedian,
  "WarmMedianSetupAmortizationRatio" ->
    mathematicaMedian/(prepareSeconds + nativeMedian),
  "ColdOneShotRatio" ->
    mathematicaColdSeconds/(prepareSeconds + nativeColdSeconds),
  "ScaledEndpointDifference" -> scaledError,
  "MathematicaEndpoint" -> mathematicaResult[[3]],
  "NativeEndpoint" -> nativeResult[[3]]
|>;
Put[record, outputFile];

Print["CHESSV2_PBB64_MATHEMATICA_TIMES=", N[mathematicaTimes, 8]];
Print["CHESSV2_PBB64_NATIVE_TIMES=", N[nativeTimes, 8]];
Print["CHESSV2_PBB64_MATHEMATICA_MEDIAN=", N[mathematicaMedian, 12]];
Print["CHESSV2_PBB64_NATIVE_MEDIAN=", N[nativeMedian, 12]];
Print["CHESSV2_PBB64_NATIVE_KERNEL_LAST=", N[kernelSeconds, 12]];
Print["CHESSV2_PBB64_PREPARE_SECONDS=", N[prepareSeconds, 12]];
Print["CHESSV2_PBB64_CACHED_SPEEDUP=", N[record["CachedSpeedup"], 12]];
Print[
  "CHESSV2_PBB64_WARM_MEDIAN_SETUP_RATIO=",
  N[record["WarmMedianSetupAmortizationRatio"], 12]
];
Print[
  "CHESSV2_PBB64_COLD_ONE_SHOT_RATIO=",
  N[record["ColdOneShotRatio"], 12]
];
Print["CHESSV2_PBB64_ENDPOINT_DIFFERENCE=", N[scaledError, 12]];
Print["CHESSV2_PBB64_NATIVE_VS_MMA_PASS=", TrueQ[scaledError < 10^-45]];
Print["CHESSV2_PBB64_OUTPUT=", outputFile];

CHESSNativeUninstall[];
Quit[If[TrueQ[scaledError < 10^-45], 0, 6]];
