#!/usr/bin/env wolframscript

(* The generic FORM/FLINT adapter can be tested without shipping a 790 MB
   process-specific SLP.  This fixture writes three exact one-limb NFLOAT01
   numbers and checks the package decoder.  If the experimental evaluator paths
   are provided through the environment, it additionally performs one real
   FLINT batch-evaluation call. *)

ClearAll["Global`*"];
$HistoryLength = 0;

testRoot = DirectoryName[DirectoryName[ExpandFileName[$InputFileName]]];
Get[FileNameJoin[{testRoot, "CHESSv2.wl"}]];

failures = {};
CheckAssertion[name_String, condition_] := If[
  TrueQ[condition],
  Print["PASS: ", name],
  AppendTo[failures, name];
  Print["FAIL: ", name]
];

fixture = FileNameJoin[{
  $TemporaryDirectory,
  "chess-flint-protocol-" <> ToString[$ProcessID] <> ".bin"
}];
stream = OpenWrite[fixture, BinaryFormat -> True];
BinaryWrite[stream, ToCharacterCode["NFLOAT01"], "UnsignedInteger8"];
BinaryWrite[
  stream,
  {1, 1, 1, 3},
  "UnsignedInteger64",
  ByteOrdering -> -1
];
(* value = (-1)^sign mantissa 2^(exponent-64). *)
BinaryWrite[
  stream,
  Flatten[{
    {2, 0, 7 2^61},       (*  7/2 *)
    {2, 1, 2^63},         (* -2   *)
    {2^63, 0, 0}          (* zero sentinel exponent -2^63 *)
  }],
  "UnsignedInteger64",
  ByteOrdering -> -1
];
Close[stream];
decoded = CHESSFlintReadNFloatBinary[fixture, 60];
Quiet @ DeleteFile[fixture];
CheckAssertion[
  "flint-nfloat01-exact-decoder",
  AssociationQ[decoded] &&
  {decoded["PointCount"], decoded["ValueCount"]} === {1, 3} &&
  TrueQ @ PossibleZeroQ[
    Max[Abs[Flatten[decoded["Values"]] - {7/2, -2, 0}]]
  ]
];

(* FLINT encodes +Infinity, -Infinity, and NaN in the exponent word directly
   above WORD_MIN.  They must never be rounded into an apparently finite matrix
   entry. *)
WriteOneRecord[record_List, suffix_String] := Module[{file, localStream},
  file = FileNameJoin[{
    $TemporaryDirectory,
    "chess-flint-special-" <> ToString[$ProcessID] <> "-" <> suffix <> ".bin"
  }];
  localStream = OpenWrite[file, BinaryFormat -> True];
  BinaryWrite[
    localStream, ToCharacterCode["NFLOAT01"], "UnsignedInteger8"
  ];
  BinaryWrite[
    localStream, {1, 1, 1, 1}, "UnsignedInteger64", ByteOrdering -> -1
  ];
  BinaryWrite[
    localStream, record, "UnsignedInteger64", ByteOrdering -> -1
  ];
  Close[localStream];
  file
];
malformedRecords = {
  {2^63 + 1, 0, 0},       (* +Infinity *)
  {2^63 + 2, 1, 0},       (* -Infinity *)
  {2^63 + 3, 0, 0},       (* NaN *)
  {2^63 + 4, 0, 2^63},    (* below NFLOAT_MIN_EXP *)
  {2^61, 0, 2^63},        (* above NFLOAT_MAX_EXP *)
  {2, 2, 2^63}            (* invalid sign *)
};
specialDecoded = Table[
  specialFixture = WriteOneRecord[
    malformedRecords[[index]], ToString[index]
  ];
  value = Quiet @ CHESSFlintReadNFloatBinary[specialFixture, 60];
  Quiet @ DeleteFile[specialFixture];
  value,
  {index, Length[malformedRecords]}
];
CheckAssertion[
  "flint-nonfinite-record-fails-closed",
  specialDecoded === ConstantArray[$Failed, Length[malformedRecords]]
];

evaluator = Environment["CHESSV2_FLINT_EVALUATOR"];
slp = Environment["CHESSV2_FLINT_SLP"];
If[StringQ[evaluator] && StringQ[slp] &&
   CHESSFlintAvailableQ[evaluator, slp],
  flintResult = CHESSFlintRunSLP[
    evaluator, slp, {{1, 2, 3, 4, 5, 6, 7}}, 40, 1
  ];
  CheckAssertion[
    "flint-existing-slp-batch-smoke",
    AssociationQ[flintResult] && flintResult["PointCount"] === 1 &&
    flintResult["ValueCount"] > 0 && Dimensions[flintResult["Values"]] ===
      {1, flintResult["ValueCount"]}
  ];
  ClearAll[AdapterCoordinates, AdapterAssembly];
  AdapterCoordinates[point_, precision_] :=
    N[{1 + point, 2, 3, 4, 5, 6, 7}, precision];
  AdapterAssembly[point_, values_, precision_] :=
    N[{{{values[[1]]}}, {{values[[2]]}}}, precision];
  flintAdapter = CHESSNativeFlintMatrixAdapter[
    evaluator, slp, AdapterCoordinates, AdapterAssembly, 1
  ];
  flintAdapterValues = CHESSNodeEvaluatorBatch[flintAdapter][
    {0, 1/10}, 40, 2
  ];
  BadAssembly[point_, values_, precision_] := {IdentityMatrix[2]};
  badFlintAdapter = CHESSNativeFlintMatrixAdapter[
    evaluator, slp, AdapterCoordinates, BadAssembly, 1
  ];
  badFlintAdapterValues = Quiet @ CHESSNodeEvaluatorBatch[badFlintAdapter][
    {0}, 40, 1
  ];
  CheckAssertion[
    "flint-native-matrix-adapter",
    CHESSNodeEvaluatorQ[flintAdapter] && ListQ[flintAdapterValues] &&
    Dimensions[flintAdapterValues] === {2, 2, 1, 1} &&
    badFlintAdapterValues === $Failed
  ],
  Print["SKIP: flint-existing-slp-batch-smoke (set CHESSV2_FLINT_EVALUATOR and CHESSV2_FLINT_SLP)"]
];

If[failures === {},
  Print["CHESSV2_FLINT_EVALUATION_TESTS_PASSED"];
  Exit[0],
  Print[InputForm[{"CHESSV2_FLINT_EVALUATION_FAILURES", failures}]];
  Exit[1]
];
