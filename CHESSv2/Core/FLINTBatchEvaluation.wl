(* ::Package:: *)

(*
  FLINTBatchEvaluation.wl

  Reusable Mathematica side of the FORM -> straight-line program -> FLINT
  nfloat evaluator tested in the non-canonical phenomenology prototype.

  This module intentionally does not know any process-specific variables or
  matrix layout.  A physics example owns three small pieces:

    1. a map from a path node to the real input coordinates of its SLP;
    2. a compiled evaluator and its prepared .chslp file;
    3. assembly of the returned scalar values into B_p matrices.

  CHESS owns only the stable batch protocol.  The resulting batch function is
  passed through the already established

      CHESSNodeEvaluator[scalar, batch]

  interface, so FLINT expression evaluation and C++ spectral propagation are
  independent modules and may be enabled separately.

  Compatible FLINT evaluator executables implement

      evaluator evaluate SLP POINTS BITS THREADS OUTPUT

  and write the little-endian NFLOAT01 format used by the experimental FLINT
  backend.  This separation avoids hard-coding a particular FORM program into
  the CHESS package.
*)

ClearAll[
  CHESSFLINTAvailableQ,
  CHESSFLINTWorkingBits,
  CHESSFLINTDecimalString,
  CHESSFLINTWritePointTable,
  CHESSFLINTReadNFloatBinary,
  CHESSFLINTRunSLP
];

CHESSFLINTRunSLP::files =
  "The FLINT evaluator or prepared SLP file is missing: `1`.";
CHESSFLINTRunSLP::points =
  "Point coordinates must be a non-empty rectangular real numerical matrix.";
CHESSFLINTRunSLP::process =
  "The FLINT evaluator failed with exit code `1`: `2`.";
CHESSFLINTRunSLP::binary =
  "The FLINT evaluator returned an unreadable or dimensionally incompatible NFLOAT01 file.";
CHESSFLINTRunSLP::option =
  "GuardDigits, WorkingBits, or Timeout has an invalid value.";

CHESSFLINTAvailableQ[executable_, slpFile_] :=
  FileExistsQ[ExpandFileName[executable]] &&
  FileExistsQ[ExpandFileName[slpFile]];

(* FLINT nfloat stores complete machine limbs.  Guard digits are a tunable
   heuristic, not an error bound: cancellation, poles, or large SLP
   intermediates may require substantially more working precision. *)
CHESSFLINTWorkingBits[precision_?NumericQ] :=
  CHESSFLINTWorkingBits[precision, 20];
CHESSFLINTWorkingBits[
  precision_?NumericQ, guardDigits_Integer?NonNegative
] := Ceiling[
  (Max[20, Ceiling[N[precision]]] + guardDigits) Log[2, 10]
];

(* Fixed-point input prevents the FLINT-side parser from depending on Mathematica's
   *^ exponent notation.  NumberPadding also makes small path coordinates
   explicit instead of silently shortening their precision. *)
CHESSFLINTDecimalString[value_?NumericQ, digits_Integer?Positive] := ToString[
  NumberForm[
    N[value, digits],
    digits,
    NumberPadding -> {"", "0"},
    NumberPoint -> ".",
    ExponentFunction -> (Null &)
  ],
  OutputForm
];

CHESSFLINTWritePointTable[
  pointRows_?MatrixQ, digits_Integer?Positive, file_String
] := Module[{stream},
  stream = Quiet @ Check[OpenWrite[file, PageWidth -> Infinity], $Failed];
  If[Head[stream] =!= OutputStream, Return[$Failed]];
  Scan[
    WriteString[stream, CHESSFLINTDecimalString[#, digits], "\n"] &,
    Flatten[pointRows]
  ];
  Close[stream];
  file
];

(* Decode FLINT's explicit {exponent,sign,limbs} representation.  Unsigned
   64-bit words are converted before reconstructing the binary mantissa, so no
   machine-real conversion occurs anywhere in the import path. *)
CHESSFLINTReadNFloatBinary[file_String, outputPrecision_] := Module[
  {
    stream, magic, header, version, limbCount, pointCount, valueCount,
    rawWords, words, signedWord, zeroExponent, minimumFiniteExponent,
    maximumFiniteExponent, decode, decodedValues
  },
  stream = Quiet @ Check[OpenRead[file, BinaryFormat -> True], $Failed];
  If[stream === $Failed, Return[$Failed]];
  magic = FromCharacterCode[BinaryReadList[stream, "UnsignedInteger8", 8]];
  header = BinaryReadList[stream, "UnsignedInteger64", 4, ByteOrdering -> -1];
  If[magic =!= "NFLOAT01" || Length[header] =!= 4,
    Close[stream];
    Return[$Failed]
  ];
  {version, limbCount, pointCount, valueCount} = header;
  rawWords = BinaryReadList[stream, "UnsignedInteger64", ByteOrdering -> -1];
  Close[stream];
  If[
    version =!= 1 || limbCount <= 0 || pointCount <= 0 || valueCount <= 0 ||
    Length[rawWords] =!= pointCount valueCount (limbCount + 2),
    Return[$Failed]
  ];
  words = Partition[rawWords, limbCount + 2];
  signedWord[x_] := If[x >= 2^63, x - 2^64, x];
  zeroExponent = -2^63;
  (* These are the exact 64-bit FLINT nfloat.h constants
       NFLOAT_MIN_EXP = WORD_MIN/4, NFLOAT_MAX_EXP = WORD_MAX/4.
     Checking the complete finite interval before constructing 2^exponent also
     prevents malformed files from producing Mathematica Underflow/Overflow
     objects which could be mistaken for numerical matrix entries. *)
  minimumFiniteExponent = -2^61;
  maximumFiniteExponent = 2^61 - 1;
  decode[row_] := Module[{exponent, sign, mantissa},
    exponent = signedWord[row[[1]]];
    sign = row[[2]];
    If[!MemberQ[{0, 1}, sign], Return[$Failed]];
    If[exponent === zeroExponent, Return[0]];
    (* Every nonzero finite value must use FLINT's documented exponent range.
       This simultaneously rejects +Infinity, -Infinity, NaN, unused special
       encodings, and impossible positive exponents. *)
    If[
      exponent < minimumFiniteExponent || exponent > maximumFiniteExponent,
      Return[$Failed]
    ];
    mantissa = FromDigits[Reverse[row[[3 ;;]]], 2^64];
    If[mantissa === 0, Return[$Failed]];
    (-1)^sign mantissa 2^(exponent - 64 limbCount)
  ];
  decodedValues = decode /@ words;
  If[MemberQ[decodedValues, $Failed], Return[$Failed]];
  <|
    "Version" -> version,
    "LimbCount" -> limbCount,
    "PointCount" -> pointCount,
    "ValueCount" -> valueCount,
    "Values" -> Partition[
      N[decodedValues, outputPrecision], valueCount
    ]
  |>
];

Options[CHESSFLINTRunSLP] = {
  "GuardDigits" -> 20,
  "WorkingBits" -> Automatic,
  "Timeout" -> Infinity
};

(* Evaluate every row in one FLINT process call.  Setup/FORM compilation is outside
   this hot path: callers should prepare the SLP once, then reuse it for every
   propagation.  The returned Association keeps evaluator stdout and timing
   available for diagnostics without changing the batch evaluator's data. *)
CHESSFLINTRunSLP[
  executable_String, slpFile_String, pointRows_?MatrixQ,
  precision_?NumericQ, threadCount_ : 1, OptionsPattern[]
] := Module[
  {
    executablePath, slpPath, dimensions, threads, digits, bits, tag,
    pointsFile, valuesFile, process, seconds, decoded, cleanup,
    guardDigits, workingBits, timeout, command
  },
  executablePath = ExpandFileName[executable];
  slpPath = ExpandFileName[slpFile];
  If[!FileExistsQ[executablePath],
    Message[CHESSFLINTRunSLP::files, executablePath];
    Return[$Failed]
  ];
  If[!FileExistsQ[slpPath],
    Message[CHESSFLINTRunSLP::files, slpPath];
    Return[$Failed]
  ];
  dimensions = Dimensions[pointRows];
  If[
    Length[dimensions] =!= 2 || Times @@ dimensions <= 0 ||
    !VectorQ[Flatten[pointRows], NumericQ] ||
    !FreeQ[pointRows, _Complex],
    Message[CHESSFLINTRunSLP::points];
    Return[$Failed]
  ];
  threads = Max[1, Min[Length[pointRows], Replace[threadCount, {
    value_Integer?Positive :> value,
    _ :> 1
  }]]];
  digits = Max[20, Ceiling[N[precision]]];
  guardDigits = OptionValue["GuardDigits"];
  workingBits = OptionValue["WorkingBits"];
  timeout = OptionValue["Timeout"];
  If[!IntegerQ[guardDigits] || guardDigits < 0 ||
      !(workingBits === Automatic || IntegerQ[workingBits] && workingBits > 0) ||
      !(timeout === Infinity || NumericQ[timeout] && TrueQ[timeout > 0]),
    Message[CHESSFLINTRunSLP::option];
    Return[$Failed]
  ];
  bits = Replace[
    workingBits,
    Automatic :> CHESSFLINTWorkingBits[precision, guardDigits]
  ];
  tag = StringJoin[
    ToString[$ProcessID], "-",
    IntegerString[
      Hash[{AbsoluteTime[], $ProcessID, RandomInteger[2^31 - 1]}, "CRC32"]
    ]
  ];
  pointsFile = FileNameJoin[{
    $TemporaryDirectory, "chess-flint-points-" <> tag <> ".txt"
  }];
  valuesFile = FileNameJoin[{
    $TemporaryDirectory, "chess-flint-values-" <> tag <> ".bin"
  }];
  cleanup[] := Quiet @ Check[
    DeleteFile /@ Select[{pointsFile, valuesFile}, FileExistsQ],
    Null
  ];

  CheckAbort[
    If[
      CHESSFLINTWritePointTable[
        pointRows, digits + guardDigits, pointsFile
      ] === $Failed,
      cleanup[];
      Return[$Failed]
    ];
    command = {
      executablePath, "evaluate", slpPath, pointsFile,
      ToString[bits], ToString[threads], valuesFile
    };
    seconds = AbsoluteTiming[
      process = Quiet @ Check[
        If[
          timeout === Infinity,
          RunProcess[command],
          TimeConstrained[RunProcess[command], timeout, $Aborted]
        ],
        $Failed
      ]
    ][[1]];
    If[process === $Failed || Lookup[process, "ExitCode", 1] =!= 0,
      Message[
        CHESSFLINTRunSLP::process,
        If[AssociationQ[process], Lookup[process, "ExitCode", $Failed], $Failed],
        If[AssociationQ[process], Lookup[process, "StandardError", ""], ""]
      ];
      cleanup[];
      Return[$Failed]
    ];
    decoded = CHESSFLINTReadNFloatBinary[valuesFile, digits],
    cleanup[];
    Abort[]
  ];
  cleanup[];
  If[
    !AssociationQ[decoded] ||
    Lookup[decoded, "PointCount", -1] =!= Length[pointRows],
    Message[CHESSFLINTRunSLP::binary];
    Return[$Failed]
  ];
  Join[
    decoded,
    <|
      "Bits" -> bits,
      "Threads" -> threads,
      "FLINTProcessSeconds" -> seconds,
      "StandardOutput" -> Lookup[process, "StandardOutput", ""]
    |>
  ]
];
