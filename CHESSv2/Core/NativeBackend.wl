(* ::Package:: *)

(*
  NativeBackend.wl

  Thin Mathematica interface to the optional C++ spectral-propagation backend.
  The backend owns only floating-point work:

    - multiprecision LU factorization of the scalar Lobatto block;
    - sparse node-matrix/vector products;
    - repeated sequential epsilon or fake-delta layer solves;
    - direct active-support LU for regular mixed polynomial systems.

  Mathematica retains the public API, route classification, matrix-evaluator
  semantics, automatic fake-delta stopping rule, and result assembly.  This is
  deliberately a capability module rather than a replacement package.  The
  direct mixed-polynomial path mirrors the active-support Mathematica algorithm
  and remains restricted to regular endpoints; singular endpoint
  regularization stays with the verified Mathematica core.

  The executable is not committed.  Build it with

      make -C CHESSv2/Native/Propagation

  CHESSNativePrepare performs matrix evaluation and setup once.  Its returned
  handle can be supplied to repeated propagations so the expensive node-matrix
  construction and LU factorization are not charged on every boundary vector.
*)

ClearAll[
  CHESSNativeAvailableQ,
  CHESSNativeInstall,
  CHESSNativeUninstall,
  CHESSNativeClear,
  CHESSNativePrepare,
  CHESSNativePolynomialPrepare,
  CHESSNativeRun,
  CHESSNativePolynomialRun,
  CHESSNativeFetch,
  CHESSNativeCanonicalResult,
  CHESSNativePolynomialResult,
  CHESSNativePolynomialPropagate,
  CHESSNativeFakeDeltaRunAtOrder,
  CHESSNativeFakeDeltaAssembleResult,
  CHESSNativeFakeDeltaPropagateSingle,
  CHESSNativeFakeDeltaPropagate,
  CHESSNativeNumberString,
  CHESSNativeWriteComplex,
  CHESSNativeParseNumber,
  CHESSNativeParseComplexValues,
  CHESSNativeValuePrecisionQ,
  CHESSNativeDataPrecisionQ,
  CHESSNativeMatrixPrecisionQ,
  CHESSNativeHandleValidQ,
  CHESSNativeHandleCompatibleQ,
  CHESSNativePolynomialHandleCompatibleQ,
  CHESSNativeWriteSparseMatrix,
  CHESSNativeDefinitionSnapshot,
  CHESSNativeReferencedSymbols,
  CHESSNativeDefinitionClosure,
  CHESSNativeEvaluatorToken,
  CHESSNativeOptionValue
];

CHESSNativeInstall::missing =
  "The optional native executable was not found at `1`. Run make in Native/Propagation first.";
CHESSNativePrepare::matrix =
  "Native preparation could not evaluate numerical `1` by `1` operators on all ordinary Lobatto nodes.";
CHESSNativePrepare::option =
  "Invalid native preparation option values: Nodes=`1`, Precision=`2`, WorkingPrecisionA=`3`, NativeGuardDigits=`4`.";
CHESSNativePrepare::precision =
  "Native `1` data contain nonzero approximate values below the required precision `2`.";
CHESSNativePolynomialPrepare::matrix =
  "Native polynomial preparation could not obtain degree `1` numerical `2` by `2` coefficient matrices at every ordinary Lobatto node.";
CHESSNativePolynomialPrepare::active =
  "Direct native polynomial transport requires a non-empty B0 active support.";
CHESSNativeRun::handle =
  "The native handle is no longer active. Preparing another system invalidates earlier handles because one backend process caches one collocation system.";
CHESSNativeRun::boundary =
  "Native boundary tensor has dimensions `1`; expected {columns,layers,`2`}.";
CHESSNativeRun::protocol =
  "The native backend returned an invalid or incompatible response.";
CHESSNativeRun::precision =
  "Native boundary data contain nonzero approximate values below the requested output precision `1`.";

Module[{packageDirectory},
  packageDirectory = DirectoryName[
    DirectoryName[ExpandFileName[$InputFileName]]
  ];
  $CHESSNativeExecutable = FileNameJoin[{
    packageDirectory, "Native", "Propagation", "chess_native_link"
  }];
];

(* Get may be called more than once in a notebook.  Preserve a live WSTP link
   and its monotonically increasing generation instead of losing ownership of
   the old process and allowing generation numbers to repeat. *)
If[!ValueQ[$CHESSNativeLink], $CHESSNativeLink = None];
If[!IntegerQ[$CHESSNativeGeneration], $CHESSNativeGeneration = 0];

CHESSNativeAvailableQ[] := FileExistsQ[$CHESSNativeExecutable];

CHESSNativeInstall[] := (
  If[MatchQ[$CHESSNativeLink, _LinkObject], Return[$CHESSNativeLink]];
  If[!CHESSNativeAvailableQ[],
    Message[CHESSNativeInstall::missing, $CHESSNativeExecutable];
    Return[$Failed]
  ];
  $CHESSNativeLink = Quiet @ Check[Install[$CHESSNativeExecutable], $Failed];
  If[!MatchQ[$CHESSNativeLink, _LinkObject],
    $CHESSNativeLink = None;
    Return[$Failed]
  ];
  $CHESSNativeLink
);

CHESSNativeUninstall[] := (
  If[MatchQ[$CHESSNativeLink, _LinkObject],
    Quiet @ Check[CHESSNativeBackend`Link`Private`CHESSNativeClear[], Null];
    Quiet @ Check[Uninstall[$CHESSNativeLink], Null]
  ];
  $CHESSNativeLink = None;
  $CHESSNativeGeneration++;
  Null
);

(* Clear only the cached collocation system, leaving the WSTP process alive so
   a later Prepare does not pay process startup again. *)
CHESSNativeClear[] := (
  If[MatchQ[$CHESSNativeLink, _LinkObject],
    Quiet @ Check[CHESSNativeBackend`Link`Private`CHESSNativeClear[], Null]
  ];
  $CHESSNativeGeneration++;
  Null
);

CHESSNativeNumberString[value_, precision_] :=
  ToString[FortranForm[N[value, precision]]];

CHESSNativeWriteComplex[stream_, value_, precision_] := WriteString[
  stream,
  CHESSNativeNumberString[Re[value], precision], " ",
  CHESSNativeNumberString[Im[value], precision], "\n"
];

CHESSNativeWriteSparseMatrix[stream_, matrix_, precision_] := Module[
  {rules, entries},
  rules = ArrayRules[SparseArray[matrix]];
  entries = Cases[
    rules,
    HoldPattern[{row_Integer, column_Integer} -> value_] :>
      {row - 1, column - 1, value}
  ];
  WriteString[stream, Length[entries], "\n"];
  Scan[
    Function[entry,
      WriteString[stream, entry[[1]], " ", entry[[2]], " "];
      CHESSNativeWriteComplex[stream, entry[[3]], precision]
    ],
    entries
  ];
  Length[entries]
];

CHESSNativeParseNumber[token_, precision_] := SetPrecision[
  ToExpression[StringReplace[token, {"e" -> "*^", "E" -> "*^"}], InputForm],
  precision
];

CHESSNativeParseComplexValues[tokens_List, precision_] := MapThread[
  Complex,
  Transpose @ Partition[CHESSNativeParseNumber[#, precision] & /@ tokens, 2]
];

(* Never manufacture precision by serializing a machine number with a longer
   decimal format.  Exact numbers and numerical zeros are safe at any requested
   precision; every other value must carry enough genuine arbitrary-precision
   information before it may enter the Native backend. *)
CHESSNativeValuePrecisionQ[value_?NumberQ, requiredPrecision_Integer] :=
  TrueQ[value == 0] || Precision[value] === Infinity ||
  TrueQ[Precision[value] >= requiredPrecision];
CHESSNativeValuePrecisionQ[_, _] := False;

CHESSNativeDataPrecisionQ[data_, requiredPrecision_Integer] := And @@ (
  CHESSNativeValuePrecisionQ[#, requiredPrecision] & /@ Flatten[data]
);

(* SparseArray's final ArrayRules item is the default zero.  Check only stored
   nonzero entries so approximate 0. values do not reject an otherwise exact
   sparse operator. *)
CHESSNativeMatrixPrecisionQ[matrix_, requiredPrecision_Integer] := Module[
  {values},
  values = Cases[
    Most[ArrayRules[SparseArray[matrix]]],
    HoldPattern[_ -> value_] :> value
  ];
  And @@ (
    CHESSNativeValuePrecisionQ[#, requiredPrecision] & /@ values
  )
];

CHESSNativeOptionValue[rules_List, name_, default_] :=
  Replace[name /. rules, name -> default];

(* Keep symbols unevaluated while recording their definitions.  Hashing only a
   symbol name would miss the common notebook operation Clear[f]; f[x_]:=... . *)
SetAttributes[CHESSNativeDefinitionSnapshot, HoldAllComplete];
CHESSNativeDefinitionSnapshot[symbol_] := {
  HoldComplete[symbol],
  OwnValues[symbol],
  DownValues[symbol],
  SubValues[symbol],
  UpValues[symbol],
  NValues[symbol],
  DefaultValues[symbol],
  Options[Unevaluated[symbol]],
  Attributes[symbol]
};

(* Return held names so inspecting a dependency never triggers its OwnValues.
   Heads -> True is important because function names occur as expression heads,
   not ordinary arguments. *)
CHESSNativeReferencedSymbols[expression_] := DeleteDuplicates @ Cases[
  Unevaluated[expression],
  symbol_Symbol /; Context[Unevaluated[symbol]] =!= "System`" :>
    HoldComplete[symbol],
  Infinity,
  Heads -> True
];

(* Follow helper functions to a fixed point.  A one-level snapshot is unsafe:
   f[x_]:=g[x] does not change when g is redefined, although the sampled
   matrices do.  The queue is cycle-safe, so mutually recursive definitions do
   not loop.  This is deliberately conservative: an irrelevant user symbol in
   a held definition may invalidate a handle, but a stale handle is never used. *)
CHESSNativeDefinitionClosure[rootSymbols_List] := Module[
  {queue, seen, current, snapshot, dependencies, closure},
  queue = DeleteDuplicates[rootSymbols];
  seen = {};
  closure = {};
  While[Length[queue] > 0,
    current = First[queue];
    queue = Rest[queue];
    If[MemberQ[seen, current], Continue[]];
    AppendTo[seen, current];
    snapshot = current /. HoldComplete[symbol_Symbol] :>
      CHESSNativeDefinitionSnapshot[symbol];
    AppendTo[closure, snapshot];
    dependencies = CHESSNativeReferencedSymbols[snapshot];
    queue = Join[
      queue,
      Select[dependencies, ! MemberQ[seen, #] && ! MemberQ[queue, #] &]
    ];
  ];
  closure
];

(* Bind an opaque prepared handle to the evaluator expression and the current
   transitive definition closure of every non-System symbol referenced by it.
   Following helper functions is necessary because changing g in f[x_]:=g[x]
   changes the sampled equation without changing f's DownValues.  This catches
   both a different evaluator and ordinary redefinitions without re-evaluating
   all nodes.  $CHESSCacheGeneration additionally tracks the package's explicit
   cache invalidation operation.  Arbitrary mutable external state which is
   invisible in Wolfram definitions remains the caller's responsibility. *)
CHESSNativeEvaluatorToken[matrixSpec_] := Module[
  {held, symbols, definitions, cacheGeneration},
  held = HoldComplete[matrixSpec];
  symbols = CHESSNativeReferencedSymbols[held];
  definitions = CHESSNativeDefinitionClosure[symbols];
  cacheGeneration = If[
    ValueQ[$CHESSCacheGeneration], $CHESSCacheGeneration, Missing["Absent"]
  ];
  Hash[{held, definitions, cacheGeneration}, "SHA256"]
];

CHESSNativeHandleValidQ[handle_] := AssociationQ[handle] &&
  Lookup[handle, "Generation", Missing["Generation"]] === $CHESSNativeGeneration &&
  MatchQ[$CHESSNativeLink, _LinkObject];

CHESSNativeHandleCompatibleQ[
  handle_, matrixSpec_, dimension_Integer, interval : {_, _},
  canonicalRules_List
] := Module[{nodes, precision, precisionA, guardDigits, handleInterval},
  If[!CHESSNativeHandleValidQ[handle], Return[False]];
  nodes = CHESSNativeOptionValue[canonicalRules, "Nodes", 48];
  precision = CHESSNativeOptionValue[canonicalRules, "Precision", 160];
  precisionA = CHESSNativeOptionValue[
    canonicalRules, "WorkingPrecisionA", 220
  ];
  guardDigits = CHESSNativeOptionValue[
    canonicalRules, "NativeGuardDigits", 20
  ];
  If[!IntegerQ[guardDigits] || guardDigits < 0, Return[False]];
  handleInterval = Lookup[handle, "Interval", Missing["Interval"]];
  Lookup[handle, "Dimension", Missing["Dimension"]] === dimension &&
  Lookup[handle, "NodesOption", Missing["NodesOption"]] === nodes &&
  Lookup[handle, "Precision", Missing["Precision"]] === precision &&
  Lookup[handle, "NativeGuardDigits", Missing["NativeGuardDigits"]] ===
    guardDigits &&
  Lookup[handle, "NativePrecision", Missing["NativePrecision"]] ===
    precision + guardDigits &&
  Lookup[handle, "WorkingPrecisionA", Missing["WorkingPrecisionA"]] ===
    precisionA &&
  Lookup[handle, "EvaluatorToken", Missing["EvaluatorToken"]] ===
    CHESSNativeEvaluatorToken[matrixSpec] &&
  ListQ[handleInterval] && Length[handleInterval] === 2 &&
  And @@ MapThread[
    TrueQ[PossibleZeroQ[N[#1 - #2, precision]]] &,
    {handleInterval, interval}
  ]
];

CHESSNativePolynomialHandleCompatibleQ[
  handle_, matrixSpec_, dimension_Integer, interval : {_, _}, rules_List,
  degree_Integer
] := CHESSNativeHandleCompatibleQ[
    handle, matrixSpec, dimension, interval, rules
  ] &&
  Lookup[handle, "Mode", Missing["Mode"]] === "Polynomial" &&
  Lookup[handle, "Degree", Missing["Degree"]] === degree;

(* Prepare one regular-point canonical collocation system.  The common matrix
   adapter is used here as well, so native and Mathematica paths accept exactly
   the same evaluator signatures and enforce the same n x n dimensions. *)
CHESSNativePrepare[
  matrixSpec_, dimension_Integer?Positive, interval : {_, _},
  canonicalRules_List
] := Module[
  {
    nodesOption, precision, nativePrecision, precisionA, nativePrecisionA,
    guardDigits, parallelSpec, kernelSpec, threads,
    collocation, nodes, publicNodes, scalarMatrix, evaluator, nodeMatrices,
    setupFile,
    stream, rules, entries, loadSeconds, loadResult, tag, Cleanup
  },
  nodesOption = CHESSNativeOptionValue[canonicalRules, "Nodes", 48];
  precision = CHESSNativeOptionValue[canonicalRules, "Precision", 160];
  precisionA = CHESSNativeOptionValue[
    canonicalRules, "WorkingPrecisionA", 220
  ];
  guardDigits = CHESSNativeOptionValue[
    canonicalRules, "NativeGuardDigits", 20
  ];
  nativePrecision = precision + guardDigits;
  (* Ask matrix evaluators for an additional ten digits when the guarded
     transport precision would otherwise determine the sampling precision.
     Ordinary arithmetic may lose a few digits before the completed matrix is
     returned; the nonzero-entry check below still requires the full guarded
     transport precision actually to survive. *)
  nativePrecisionA = Max[precisionA, nativePrecision + 10];
  parallelSpec = CHESSNativeOptionValue[
    canonicalRules, "ParallelEvaluation", Automatic
  ];
  kernelSpec = CHESSNativeOptionValue[canonicalRules, "ParallelKernels", 1];
  If[
    !IntegerQ[nodesOption] || nodesOption <= 0 ||
    !IntegerQ[precision] || precision <= 0 ||
    !IntegerQ[precisionA] || precisionA < precision ||
    !IntegerQ[guardDigits] || guardDigits < 0,
    Message[
      CHESSNativePrepare::option,
      nodesOption, precision, precisionA, guardDigits
    ];
    Return[$Failed]
  ];
  If[!CHESSNativeDataPrecisionQ[interval, precision],
    Message[CHESSNativePrepare::precision, "interval", precision];
    Return[$Failed]
  ];
  If[!MatchQ[CHESSNativeInstall[], _LinkObject], Return[$Failed]];

  collocation = ChebyshevLobattoData[
    nodesOption, interval, nativePrecisionA
  ];
  nodes = N[collocation[[1]], nativePrecisionA];
  publicNodes = N[nodes, precision];
  scalarMatrix = N[
    Normal @ BaseScalarCollocationMatrix[collocation[[2]]], nativePrecision
  ];
  evaluator = CHESSCanonicalNodeEvaluator[
    matrixSpec, {dimension, dimension}
  ];
  threads = CHESSBatchThreadCount[parallelSpec, nodesOption, kernelSpec];
  nodeMatrices = Quiet @ Check[
    CHESSNodeEvaluatorBatch[evaluator][
      Rest[nodes], nativePrecisionA, threads
    ],
    $Failed
  ];
  If[
    nodeMatrices === $Failed || !ListQ[nodeMatrices] ||
    Length[nodeMatrices] =!= nodesOption ||
    !AllTrue[
      nodeMatrices,
      CHESSMatrixMatchesDimensionsQ[#, {dimension, dimension}] &&
        CHESSNativeMatrixPrecisionQ[#, nativePrecision] &
    ],
    If[
      ListQ[nodeMatrices] && Length[nodeMatrices] === nodesOption &&
      AllTrue[
        nodeMatrices,
        CHESSMatrixMatchesDimensionsQ[#, {dimension, dimension}] &
      ],
      Message[
        CHESSNativePrepare::precision, "matrix", nativePrecision
      ],
      Message[CHESSNativePrepare::matrix, dimension]
    ];
    Return[$Failed]
  ];

  tag = IntegerString[
    Hash[{AbsoluteTime[], $ProcessID, RandomInteger[2^31 - 1]}, "CRC32"]
  ];
  setupFile = FileNameJoin[{
    $TemporaryDirectory,
    "chess-native-" <> ToString[$ProcessID] <> "-" <> tag <> ".dat"
  }];
  stream = None;
  Cleanup[] := (
    If[Head[stream] === OutputStream,
      Quiet @ Check[Close[stream], Null];
      stream = None
    ];
    If[StringQ[setupFile] && FileExistsQ[setupFile],
      Quiet @ Check[DeleteFile[setupFile], Null]
    ]
  );
  CheckAbort[
    stream = Quiet @ Check[
      OpenWrite[setupFile, PageWidth -> Infinity], $Failed
    ];
    If[Head[stream] =!= OutputStream,
      Cleanup[];
      Return[$Failed]
    ];
    WriteString[
      stream, "CHESSCPP1 ", nativePrecision, " ", dimension, " ",
      Length[nodes], "\n"
    ];
    Scan[
      CHESSNativeWriteComplex[stream, #, nativePrecision] &,
      Flatten[scalarMatrix]
    ];
    Do[
      rules = ArrayRules[SparseArray[nodeMatrices[[node]]]];
      entries = Cases[
        rules,
        HoldPattern[{row_Integer, column_Integer} -> value_] :>
          {row - 1, column - 1, value}
      ];
      WriteString[stream, Length[entries], "\n"];
      Scan[
        Function[entry,
          WriteString[stream, entry[[1]], " ", entry[[2]], " "];
          CHESSNativeWriteComplex[stream, entry[[3]], nativePrecision]
        ],
        entries
      ],
      {node, Length[nodeMatrices]}
    ];
    Close[stream];
    stream = None;
    {loadSeconds, loadResult} = AbsoluteTiming[
      CHESSNativeBackend`Link`Private`CHESSNativeLoad[setupFile]
    ],
    Cleanup[];
    Abort[]
  ];
  Cleanup[];
  If[loadResult =!= Null, Return[$Failed]];
  $CHESSNativeGeneration++;
  <|
    "Generation" -> $CHESSNativeGeneration,
    "Dimension" -> dimension,
    "Nodes" -> publicNodes,
    "NodesOption" -> nodesOption,
    "Precision" -> precision,
    "NativePrecision" -> nativePrecision,
    "NativeGuardDigits" -> guardDigits,
    "WorkingPrecisionA" -> precisionA,
    "NativeWorkingPrecisionA" -> nativePrecisionA,
    "EvaluatorToken" -> CHESSNativeEvaluatorToken[matrixSpec],
    "Interval" -> interval,
    "LoadSeconds" -> loadSeconds
  |>
];

(* Prepare the direct active-support operator for a regular mixed polynomial
   system.  Matrix evaluation remains on the Mathematica/FORM/FLINT side.  The
   8 x 97 F3 active block, for example, is assembled once and factorized by the
   C++ backend with FLINT nfloat; positive epsilon powers are stored only as
   sparse node operators for sequential RHS construction. *)
CHESSNativePolynomialPrepare[
  matrixSpec_CHESSNodeEvaluator, boundary_?MatrixQ, interval : {_, _},
  rules_List
] := Module[
  {
    nodesOption, precision, nativePrecision, precisionA, nativePrecisionA,
    guardDigits, parallelSpec, kernelSpec, degree,
    dimension, threads, collocation, nodes, publicNodes, scalarMatrix,
    ordinaryValues,
    matricesByNode, zeroMatrices, b0Nodes, active, flintBits, setupFile, stream,
    loadSeconds, loadResult, tag, Cleanup
  },
  nodesOption = CHESSNativeOptionValue[rules, "Nodes", 48];
  precision = CHESSNativeOptionValue[rules, "Precision", 160];
  precisionA = CHESSNativeOptionValue[rules, "WorkingPrecisionA", 220];
  guardDigits = CHESSNativeOptionValue[rules, "NativeGuardDigits", 20];
  nativePrecision = precision + guardDigits;
  nativePrecisionA = Max[precisionA, nativePrecision + 10];
  parallelSpec = CHESSNativeOptionValue[
    rules, "ParallelEvaluation", Automatic
  ];
  kernelSpec = CHESSNativeOptionValue[rules, "ParallelKernels", 1];
  degree = CHESSNativeOptionValue[rules, "EpsilonDegree", Automatic];
  dimension = First[Dimensions[boundary]];
  If[
    !IntegerQ[nodesOption] || nodesOption <= 0 ||
    !IntegerQ[precision] || precision <= 0 ||
    !IntegerQ[precisionA] || precisionA < precision ||
    !IntegerQ[guardDigits] || guardDigits < 0 ||
    !IntegerQ[degree] || degree <= 0,
    Message[
      CHESSNativePrepare::option,
      nodesOption, precision, precisionA, guardDigits
    ];
    Return[$Failed]
  ];
  If[!CHESSNativeDataPrecisionQ[interval, precision],
    Message[CHESSNativePrepare::precision, "interval", precision];
    Return[$Failed]
  ];
  If[!MatchQ[CHESSNativeInstall[], _LinkObject], Return[$Failed]];

  collocation = ChebyshevLobattoData[
    nodesOption, interval, nativePrecisionA
  ];
  nodes = N[collocation[[1]], nativePrecisionA];
  publicNodes = N[nodes, precision];
  scalarMatrix = N[
    Normal @ BaseScalarCollocationMatrix[collocation[[2]]], nativePrecision
  ];
  threads = CHESSBatchThreadCount[parallelSpec, nodesOption, kernelSpec];
  ordinaryValues = Quiet @ Check[
    CHESSNodeEvaluatorBatch[matrixSpec][
      Rest[nodes], nativePrecisionA, threads
    ],
    $Failed
  ];
  If[
    ordinaryValues === $Failed || !ListQ[ordinaryValues] ||
    Length[ordinaryValues] =!= nodesOption ||
    !And @@ Map[
      Function[nodeMatrices,
        ListQ[nodeMatrices] && Length[nodeMatrices] >= degree + 1 &&
        And @@ (
          MatrixQ[#, NumericQ] &&
            Dimensions[#] === {dimension, dimension} &&
            CHESSNativeMatrixPrecisionQ[#, nativePrecision] & /@
              Take[nodeMatrices, degree + 1]
        )
      ],
      ordinaryValues
    ],
    Message[CHESSNativePolynomialPrepare::matrix, degree, dimension];
    Return[$Failed]
  ];
  ordinaryValues = Map[
    Function[nodeMatrices,
      SparseArray[N[#, nativePrecisionA]] & /@
        Take[nodeMatrices, degree + 1]
    ],
    ordinaryValues
  ];
  zeroMatrices = ConstantArray[
    SparseArray[{}, {dimension, dimension}], degree + 1
  ];
  matricesByNode = Join[{zeroMatrices}, ordinaryValues];

  b0Nodes = matricesByNode[[All, 1]];
  active = CHESSNonCanonicalActiveSupport[b0Nodes];
  If[active === {},
    Message[CHESSNativePolynomialPrepare::active];
    Return[$Failed]
  ];
  flintBits = Ceiling[(nativePrecisionA + 20) Log[2, 10]];

  tag = IntegerString[
    Hash[{AbsoluteTime[], $ProcessID, RandomInteger[2^31 - 1]}, "CRC32"]
  ];
  setupFile = FileNameJoin[{
    $TemporaryDirectory,
    "chess-native-polynomial-" <> ToString[$ProcessID] <> "-" <> tag <>
      ".dat"
  }];
  stream = None;
  Cleanup[] := (
    If[Head[stream] === OutputStream,
      Quiet @ Check[Close[stream], Null];
      stream = None
    ];
    If[StringQ[setupFile] && FileExistsQ[setupFile],
      Quiet @ Check[DeleteFile[setupFile], Null]
    ]
  );
  CheckAbort[
    stream = Quiet @ Check[
      OpenWrite[setupFile, PageWidth -> Infinity], $Failed
    ];
    If[Head[stream] =!= OutputStream,
      Cleanup[];
      Return[$Failed]
    ];
    WriteString[
      stream, "CHESSCPP2 ", nativePrecision, " ", dimension, " ",
      Length[nodes], "\n"
    ];
    Scan[
      CHESSNativeWriteComplex[stream, #, nativePrecision] &,
      Flatten[scalarMatrix]
    ];
    WriteString[
      stream, degree, " ", Length[active], " ", flintBits, "\n",
      StringRiffle[ToString /@ (active - 1), " "], "\n"
    ];
    Do[
      CHESSNativeWriteSparseMatrix[
        stream, matricesByNode[[node, 1]], nativePrecisionA + 20
      ],
      {node, 2, Length[nodes]}
    ];
    Do[
      CHESSNativeWriteSparseMatrix[
        stream, matricesByNode[[node, power + 1]], nativePrecisionA + 20
      ],
      {power, 1, degree}, {node, 2, Length[nodes]}
    ];
    Close[stream];
    stream = None;
    {loadSeconds, loadResult} = AbsoluteTiming[
      CHESSNativeBackend`Link`Private`CHESSNativeLoad[setupFile]
    ],
    Cleanup[];
    Abort[]
  ];
  Cleanup[];
  If[loadResult =!= Null, Return[$Failed]];
  $CHESSNativeGeneration++;
  <|
    "Generation" -> $CHESSNativeGeneration,
    "Mode" -> "Polynomial",
    "Degree" -> degree,
    "ActiveCount" -> Length[active],
    "Dimension" -> dimension,
    "Nodes" -> publicNodes,
    "NodesOption" -> nodesOption,
    "Precision" -> precision,
    "NativePrecision" -> nativePrecision,
    "NativeGuardDigits" -> guardDigits,
    "WorkingPrecisionA" -> precisionA,
    "NativeWorkingPrecisionA" -> nativePrecisionA,
    "EvaluatorToken" -> CHESSNativeEvaluatorToken[matrixSpec],
    "Interval" -> interval,
    "LoadSeconds" -> loadSeconds
  |>
];

(* boundaryTensor has shape {physical columns,layers,dimension}.  The native
   protocol stores endpoint coefficients and all node states from this call;
   Fetch can subsequently return either coefficient-resolved or delta-summed
   node data without recomputing the propagation. *)
CHESSNativeRun[
  handle_Association, boundaryTensor_List, cacheStates_: True
] := Module[
  {
    dimension, precision, outputPrecision, dimensions, columns, layers, boundaryText,
    raw, tokens, header, valueTokens, values, endpointLayers, callSeconds
  },
  If[!CHESSNativeHandleValidQ[handle],
    Message[CHESSNativeRun::handle];
    Return[$Failed]
  ];
  If[!BooleanQ[cacheStates],
    Message[CHESSNativeRun::protocol];
    Return[$Failed]
  ];
  dimension = Lookup[handle, "Dimension"];
  outputPrecision = Lookup[handle, "Precision"];
  precision = Lookup[
    handle, "NativePrecision", outputPrecision
  ];
  dimensions = Dimensions[boundaryTensor];
  If[Length[dimensions] =!= 3 || Last[dimensions] =!= dimension,
    Message[CHESSNativeRun::boundary, dimensions, dimension];
    Return[$Failed]
  ];
  If[!CHESSNativeDataPrecisionQ[boundaryTensor, outputPrecision],
    Message[CHESSNativeRun::precision, outputPrecision];
    Return[$Failed]
  ];
  {columns, layers} = Take[dimensions, 2];
  callSeconds = AbsoluteTiming[
    boundaryText = StringRiffle[
      Join[
        {ToString[dimension], ToString[layers], ToString[columns]},
        Flatten[
          ({
              CHESSNativeNumberString[Re[#], precision],
              CHESSNativeNumberString[Im[#], precision]
            } &) /@ Flatten[boundaryTensor]
        ]
      ],
      " "
    ];
    raw = CHESSNativeBackend`Link`Private`CHESSNativeRun[
      boundaryText, layers, columns, If[TrueQ[cacheStates], 1, 0]
    ];
  ][[1]];
  If[!StringQ[raw],
    Message[CHESSNativeRun::protocol];
    Return[$Failed]
  ];
  tokens = StringSplit[raw];
  If[Length[tokens] < 5 || First[tokens] =!= "CHESSNATIVE1",
    Message[CHESSNativeRun::protocol];
    Return[$Failed]
  ];
  header = ToExpression[#, InputForm] & /@ tokens[[2 ;; 5]];
  valueTokens = tokens[[6 ;;]];
  If[
    header[[2 ;;]] =!= {layers, columns, dimension} ||
    Length[valueTokens] =!= 2 columns layers dimension,
    Message[CHESSNativeRun::protocol];
    Return[$Failed]
  ];
  values = CHESSNativeParseComplexValues[valueTokens, precision];
  endpointLayers = Partition[Partition[values, dimension], layers];
  <|
    "EndpointLayers" -> endpointLayers,
    "KernelSeconds" -> header[[1]],
    "BackendCallSeconds" -> callSeconds,
    "Layers" -> layers,
    "Columns" -> columns
  |>
];

CHESSNativePolynomialRun[
  handle_Association, boundaryTensor_List, cacheStates_: True
] := Module[
  {
    dimension, precision, outputPrecision, dimensions, columns, layers,
    boundaryText, raw,
    tokens, header, valueTokens, values, endpointLayers, callSeconds,
    activeCount
  },
  If[
    !CHESSNativeHandleValidQ[handle] ||
    Lookup[handle, "Mode", Missing["Mode"]] =!= "Polynomial",
    Message[CHESSNativeRun::handle];
    Return[$Failed]
  ];
  If[!BooleanQ[cacheStates],
    Message[CHESSNativeRun::protocol];
    Return[$Failed]
  ];
  dimension = Lookup[handle, "Dimension"];
  outputPrecision = Lookup[handle, "Precision"];
  precision = Lookup[
    handle, "NativePrecision", outputPrecision
  ];
  dimensions = Dimensions[boundaryTensor];
  If[Length[dimensions] =!= 3 || Last[dimensions] =!= dimension,
    Message[CHESSNativeRun::boundary, dimensions, dimension];
    Return[$Failed]
  ];
  If[!CHESSNativeDataPrecisionQ[boundaryTensor, outputPrecision],
    Message[CHESSNativeRun::precision, outputPrecision];
    Return[$Failed]
  ];
  {columns, layers} = Take[dimensions, 2];
  callSeconds = AbsoluteTiming[
    boundaryText = StringRiffle[
      Join[
        {ToString[dimension], ToString[layers], ToString[columns]},
        Flatten[
          ({
              CHESSNativeNumberString[Re[#], precision],
              CHESSNativeNumberString[Im[#], precision]
            } &) /@ Flatten[boundaryTensor]
        ]
      ],
      " "
    ];
    raw = CHESSNativeBackend`Link`Private`CHESSNativePolynomialRun[
      boundaryText, layers, columns, If[TrueQ[cacheStates], 1, 0]
    ];
  ][[1]];
  If[!StringQ[raw],
    Message[CHESSNativeRun::protocol];
    Return[$Failed]
  ];
  tokens = StringSplit[raw];
  If[Length[tokens] < 6 || First[tokens] =!= "CHESSPOLY1",
    Message[CHESSNativeRun::protocol];
    Return[$Failed]
  ];
  header = ToExpression[#, InputForm] & /@ tokens[[2 ;; 6]];
  activeCount = Lookup[handle, "ActiveCount"];
  valueTokens = tokens[[7 ;;]];
  If[
    header[[2 ;;]] =!= {layers, columns, dimension, activeCount} ||
    Length[valueTokens] =!= 2 columns layers dimension,
    Message[CHESSNativeRun::protocol];
    Return[$Failed]
  ];
  values = CHESSNativeParseComplexValues[valueTokens, precision];
  endpointLayers = Partition[Partition[values, dimension], layers];
  <|
    "EndpointLayers" -> endpointLayers,
    "KernelSeconds" -> header[[1]],
    "BackendCallSeconds" -> callSeconds,
    "Layers" -> layers,
    "Columns" -> columns,
    "ActiveCount" -> activeCount
  |>
];

CHESSNativePolynomialResult[
  handle_Association, run_Association, resultData_
] := Module[
  {layers, outputPrecision, endpoint, fetched, states, nodeValues},
  layers = Lookup[run, "Layers"];
  outputPrecision = Lookup[handle, "Precision"];
  endpoint = N[
    Transpose[First[Lookup[run, "EndpointLayers"]]], outputPrecision
  ];
  If[resultData === "Endpoint",
    nodeValues = Missing["NotRequested"],
    fetched = CHESSNativeFetch[handle, 0, layers - 1];
    If[fetched === $Failed, Return[$Failed]];
    states = First[Lookup[fetched, "Data"]];
    nodeValues = N[
      Table[
        Flatten[states[[All, node]], 1],
        {node, Length[Lookup[handle, "Nodes"]]}
      ],
      outputPrecision
    ]
  ];
  {
    Lookup[handle, "Nodes"],
    nodeValues,
    endpoint,
    {},
    {},
    {
      "PolynomialEpsilonActiveSupport", Lookup[handle, "ActiveCount"],
      "Native", "DirectCoupledLU",
      "KernelSeconds", Lookup[run, "KernelSeconds"],
      "BackendCallSeconds", Lookup[run, "BackendCallSeconds"],
      "NativeWorkingPrecision", Lookup[handle, "NativePrecision"],
      "ResultData", resultData
    }
  }
];

CHESSNativePolynomialPropagate[
  matrixSpec_CHESSNodeEvaluator, boundary_?MatrixQ, interval_, rules_List,
  nativeHandle_,
  resultData_
] := Module[{dimension, degree, handle, run, tensor},
  If[!MemberQ[{"Full", "Endpoint"}, resultData],
    Message[SpectralPropagate::nativeresult, resultData];
    Return[$Failed]
  ];
  dimension = First[Dimensions[boundary]];
  degree = CHESSNativeOptionValue[rules, "EpsilonDegree", Automatic];
  If[degree <= 0, Return[$Failed]];
  handle = If[
    nativeHandle === Automatic,
    CHESSNativePolynomialPrepare[matrixSpec, boundary, interval, rules],
    nativeHandle
  ];
  If[
    handle === $Failed ||
    !CHESSNativePolynomialHandleCompatibleQ[
      handle, matrixSpec, dimension, interval, rules, degree
    ],
    Message[SpectralPropagate::nativehandle];
    Return[$Failed]
  ];
  tensor = {Transpose[boundary]};
  run = CHESSNativePolynomialRun[
    handle, tensor, resultData === "Full"
  ];
  If[run === $Failed, Return[$Failed]];
  CHESSNativePolynomialResult[handle, run, resultData]
];

CHESSNativeFetch[
  handle_Association, mode_Integer, selectedOrder_Integer?NonNegative
] := Module[
  {
    precision, dimension, raw, tokens, header, valueTokens, values,
    layers, columns, returnedDimension, nodeCount, expectedCount, data
  },
  If[!CHESSNativeHandleValidQ[handle],
    Message[CHESSNativeRun::handle];
    Return[$Failed]
  ];
  precision = Lookup[
    handle, "NativePrecision", Lookup[handle, "Precision"]
  ];
  dimension = Lookup[handle, "Dimension"];
  raw = CHESSNativeBackend`Link`Private`CHESSNativeFetch[mode, selectedOrder];
  If[!StringQ[raw], Return[$Failed]];
  tokens = StringSplit[raw];
  If[Length[tokens] < 7 || First[tokens] =!= "CHESSFETCH1", Return[$Failed]];
  header = ToExpression[#, InputForm] & /@ tokens[[2 ;; 7]];
  {layers, columns, returnedDimension, nodeCount} = header[[3 ;; 6]];
  If[
    header[[1]] =!= mode || header[[2]] =!= selectedOrder ||
    returnedDimension =!= dimension ||
    nodeCount =!= Length[Lookup[handle, "Nodes"]],
    Return[$Failed]
  ];
  valueTokens = tokens[[8 ;;]];
  expectedCount = If[
    mode === 0,
    2 columns (selectedOrder + 1) nodeCount dimension,
    2 columns nodeCount dimension
  ];
  If[Length[valueTokens] =!= expectedCount, Return[$Failed]];
  values = CHESSNativeParseComplexValues[valueTokens, precision];
  data = If[
    mode === 0,
    Partition[
      Partition[Partition[values, dimension], nodeCount],
      selectedOrder + 1
    ],
    Partition[Partition[values, dimension], nodeCount]
  ];
  <|
    "Mode" -> mode,
    "SelectedOrder" -> selectedOrder,
    "Data" -> data
  |>
];

(* Assemble the ordinary canonical six-part CHESS result.  Endpoint-only mode
   deliberately leaves node values as Missing; callers selecting it are asking
   for the high-throughput repeated-boundary path and must not assume the full
   historical result payload was transferred through WSTP. *)
CHESSNativeCanonicalResult[
  handle_Association, run_Association, resultData_
] := Module[
  {
    layers, outputPrecision, endpoint, fetched, states, nodeValues
  },
  layers = Lookup[run, "Layers"];
  outputPrecision = Lookup[handle, "Precision"];
  endpoint = N[
    Transpose[First[Lookup[run, "EndpointLayers"]]], outputPrecision
  ];
  If[resultData === "Endpoint",
    nodeValues = Missing["NotRequested"],
    fetched = CHESSNativeFetch[handle, 0, layers - 1];
    If[fetched === $Failed, Return[$Failed]];
    states = First[Lookup[fetched, "Data"]];
    nodeValues = N[
      Table[
        Flatten[states[[All, node]], 1],
        {node, Length[Lookup[handle, "Nodes"]]}
      ],
      outputPrecision
    ]
  ];
  {
    Lookup[handle, "Nodes"],
    nodeValues,
    endpoint,
    {},
    {},
    {
      "SequentialEpsilon", "Native",
      "KernelSeconds", Lookup[run, "KernelSeconds"],
      "BackendCallSeconds", Lookup[run, "BackendCallSeconds"],
      "NativeWorkingPrecision", Lookup[handle, "NativePrecision"],
      "ResultData", resultData
    }
  }
];

(* One B0 fake-delta call transports all physical boundary columns together.
   The returned endpoint tensor has shape {physical columns,layers,n}; transpose
   each layer to the n x q coefficient matrices used by the existing automatic
   order estimator. *)
CHESSNativeFakeDeltaRunAtOrder[
  handle_Association, boundary_?MatrixQ, order_Integer?NonNegative,
  cacheStates_: True
] := Module[{n, physicalColumns, tensor, run, endpointLayers, coefficients},
  n = Dimensions[boundary][[1]];
  physicalColumns = Dimensions[boundary][[2]];
  tensor = ConstantArray[0, {physicalColumns, order + 1, n}];
  tensor[[All, 1, All]] = Transpose[boundary];
  run = CHESSNativeRun[handle, tensor, cacheStates];
  If[run === $Failed, Return[$Failed]];
  endpointLayers = Lookup[run, "EndpointLayers"];
  coefficients = Table[
    Transpose[endpointLayers[[All, layer, All]]],
    {layer, order + 1}
  ];
  <|"Run" -> run, "CoefficientMatrices" -> coefficients|>
];

CHESSNativeFakeDeltaAssembleResult[
  handle_Association, nativeRun_Association,
  coefficientMatrices_List, selectedOrder_Integer?NonNegative,
  tailEstimate_, tolerance_, selectionMode_, computedOrder_Integer?NonNegative,
  resultData_, preserveGuardState_ : False
] := Module[
  {
    outputPrecision, nativeFinalState, finalState, fetched, columnNodeData,
    nodeMatrices, nodeValues
  },
  outputPrecision = Lookup[handle, "Precision"];
  nativeFinalState = Total[
    Take[coefficientMatrices, selectedOrder + 1]
  ];
  finalState = If[
    TrueQ[preserveGuardState],
    nativeFinalState,
    N[nativeFinalState, outputPrecision]
  ];
  If[resultData === "Endpoint",
    nodeValues = Missing["NotRequested"],
    fetched = CHESSNativeFetch[handle, 1, selectedOrder];
    If[fetched === $Failed, Return[$Failed]];
    columnNodeData = Lookup[fetched, "Data"];
    nodeMatrices = Table[
      Transpose[columnNodeData[[All, node, All]]],
      {node, Length[Lookup[handle, "Nodes"]]}
    ];
    nodeValues = N[
      Flatten[Transpose[#], 1] & /@ nodeMatrices,
      outputPrecision
    ]
  ];
  {
    Lookup[handle, "Nodes"],
    nodeValues,
    finalState,
    {},
    {},
    {
      "FakeDelta", selectionMode, "Native",
      "SelectedOrder", selectedOrder,
      "ComputedOrder", computedOrder,
      "TailEstimate", tailEstimate,
      "Tolerance", tolerance,
      "GuardLayers", If[selectionMode === "Automatic", 1, 0],
      "KernelSeconds", Lookup[Lookup[nativeRun, "Run"], "KernelSeconds"],
      "BackendCallSeconds",
        Lookup[Lookup[nativeRun, "Run"], "BackendCallSeconds"],
      "NativeWorkingPrecision", Lookup[handle, "NativePrecision"],
      "ResultData", resultData
    }
  }
];

CHESSNativeFakeDeltaPropagateSingle[
  b0_, boundary_?MatrixQ, interval_, canonicalRules_List,
  deltaOrder_, tolerance_, maxOrder_Integer, nativeHandle_, resultData_,
  preserveGuardState_ : False
] := Module[
  {
    dimension, handle, currentOrder, nativeRun, coefficientMatrices,
    accepted, selectedOrder, tailEstimate, metrics, estimates
  },
  If[!MemberQ[{"Full", "Endpoint"}, resultData],
    Message[SpectralPropagate::nativeresult, resultData];
    Return[$Failed]
  ];
  dimension = Dimensions[boundary][[1]];
  handle = If[
    nativeHandle === Automatic,
    CHESSNativePrepare[b0, dimension, interval, canonicalRules],
    nativeHandle
  ];
  If[
    handle === $Failed ||
    !CHESSNativeHandleCompatibleQ[
      handle, b0, dimension, interval, canonicalRules
    ],
    Message[SpectralPropagate::nativehandle];
    Return[$Failed]
  ];

  If[IntegerQ[deltaOrder],
    nativeRun = CHESSNativeFakeDeltaRunAtOrder[
      handle, boundary, deltaOrder, resultData === "Full"
    ];
    If[nativeRun === $Failed, Return[$Failed]];
    coefficientMatrices = Lookup[nativeRun, "CoefficientMatrices"];
    metrics = CHESSDeltaCoefficientMetrics[coefficientMatrices];
    estimates = CHESSDeltaTailEstimates[metrics];
    tailEstimate = If[
      deltaOrder + 1 <= Length[estimates],
      estimates[[deltaOrder + 1]],
      Infinity
    ];
    Return @ CHESSNativeFakeDeltaAssembleResult[
      handle, nativeRun, coefficientMatrices, deltaOrder, tailEstimate,
      tolerance, "Fixed", deltaOrder, resultData, preserveGuardState
    ]
  ];

  currentOrder = Min[8, maxOrder];
  While[True,
    nativeRun = CHESSNativeFakeDeltaRunAtOrder[
      handle, boundary, currentOrder, resultData === "Full"
    ];
    If[nativeRun === $Failed, Return[$Failed]];
    coefficientMatrices = Lookup[nativeRun, "CoefficientMatrices"];
    accepted = CHESSFindAcceptedDeltaOrder[
      coefficientMatrices, tolerance
    ];
    selectedOrder = accepted[[1]];
    tailEstimate = accepted[[2]];
    If[!MissingQ[selectedOrder],
      Return @ CHESSNativeFakeDeltaAssembleResult[
        handle, nativeRun, coefficientMatrices, selectedOrder, tailEstimate,
        tolerance, "Automatic", currentOrder, resultData,
        preserveGuardState
      ]
    ];
    If[currentOrder >= maxOrder,
      Message[
        CHESSFakeDeltaPropagate::noconvergence,
        maxOrder, tailEstimate, tolerance
      ];
      Return[$Failed]
    ];
    currentOrder = Min[2 currentOrder, maxOrder]
  ]
];

(* Segment composition mirrors the Mathematica FakeDelta module.  A prepared
   handle is interval-specific, so it is accepted only for a single segment;
   segmented native propagation prepares one cached system per subinterval. *)
CHESSNativeFakeDeltaPropagate[
  b0_, boundary_?MatrixQ, {x0_, x1_}, canonicalRules_List,
  deltaOrder_, tolerance_, maxOrder_Integer, segments_Integer?Positive,
  nativeHandle_, resultData_
] := Module[
  {
    precision, precisionA, guardDigits, nativePrecision, nativePrecisionA,
    segmentPoints, state, result, allNodes, allValues, segmentInfo,
    localHandle, segmentFailure
  },
  If[segments > 1 && nativeHandle =!= Automatic,
    Message[SpectralPropagate::nativehandle];
    Return[$Failed]
  ];
  precision = CHESSNativeOptionValue[canonicalRules, "Precision", 160];
  precisionA = CHESSNativeOptionValue[
    canonicalRules, "WorkingPrecisionA", 220
  ];
  guardDigits = CHESSNativeOptionValue[
    canonicalRules, "NativeGuardDigits", 20
  ];
  If[!IntegerQ[guardDigits] || guardDigits < 0,
    Message[
      CHESSNativePrepare::option,
      CHESSNativeOptionValue[canonicalRules, "Nodes", 48],
      precision,
      CHESSNativeOptionValue[canonicalRules, "WorkingPrecisionA", 220],
      guardDigits
    ];
    Return[$Failed]
  ];
  nativePrecision = precision + guardDigits;
  nativePrecisionA = Max[precisionA, nativePrecision + 10];
  segmentPoints = N[
    Subdivide[x0, x1, segments], nativePrecisionA
  ];
  state = boundary;
  allNodes = {};
  allValues = If[resultData === "Endpoint", Missing["NotRequested"], {}];
  segmentInfo = {};
  segmentFailure = False;
  Do[
    localHandle = If[segments === 1, nativeHandle, Automatic];
    result = CHESSNativeFakeDeltaPropagateSingle[
      b0,
      state,
      {segmentPoints[[segment]], segmentPoints[[segment + 1]]},
      canonicalRules,
      deltaOrder,
      tolerance,
      maxOrder,
      localHandle,
      resultData,
      segments > 1
    ];
    If[result === $Failed,
      segmentFailure = True;
      Break[]
    ];
    If[segment === 1,
      allNodes = result[[1]];
      If[resultData === "Full", allValues = result[[2]]],
      allNodes = Join[allNodes, Rest[result[[1]]]];
      If[resultData === "Full",
        allValues = Join[allValues, Rest[result[[2]]]]
      ]
    ];
    state = result[[3]];
    AppendTo[segmentInfo, result[[6]]],
    {segment, 1, segments}
  ];
  If[TrueQ[segmentFailure], Return[$Failed]];
  {
    allNodes,
    allValues,
    N[state, precision],
    {},
    {},
    If[
      segments === 1,
      First[segmentInfo],
      {"FakeDeltaSegments", "Native", "Segments", segments,
       "SegmentInfo", segmentInfo, "ResultData", resultData}
    ]
  }
];
