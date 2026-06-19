ClearAll["Global`*"];
$HistoryLength = 0;

scriptFile = If[StringQ[$InputFileName] && $InputFileName =!= "",
  $InputFileName,
  First[$ScriptCommandLine]
];
scriptDir = DirectoryName[ExpandFileName[scriptFile]];
runExample = FileNameJoin[{scriptDir, "run_example.wls"}];

tests = {
  {"DP", 8, 1, 60, 80, 10^-8},
  {"PBB", 8, 1, 60, 80, Missing["NoIndependentReference"]},
  {"BPB", 8, 1, 60, 80, Missing["NoIndependentReference"]},
  {"BHB", 8, 1, 60, 80, Missing["NoIndependentReference"]},
  {"PBP", 8, 1, 60, 80, Missing["NoIndependentReference"]},
  {"BHABHA", 24, 1, 80, 100, 10^-7}
};

ClearAll[extractError, formatValue, runOne];

extractError[text_] := Module[{hit},
  hit = StringCases[
    text,
    RegularExpression["maxReferenceError=([^\\s]+)"] -> "$1"
  ];
  If[hit === {} || First[hit] === "not_available",
    Missing["NotAvailable"],
    Quiet @ Check[
      ToExpression[StringReplace[First[hit], RegularExpression["[eE]"] -> "*^"]],
      Missing["ParseFailed"]
    ]
  ]
];

formatValue[x_] := If[NumberQ[x], ToString[FortranForm[N[x, 12]]], ToString[x]];

runOne[{case_, nodes_, kernels_, precision_, precisionA_, threshold_}] := Module[
  {cmd, proc, out, err, exit, referenceError, pass},
  cmd = {
    "wolframscript",
    "-file", runExample,
    case,
    ToString[nodes],
    ToString[kernels],
    ToString[precision],
    ToString[precisionA]
  };
  proc = RunProcess[cmd];
  out = proc["StandardOutput"];
  err = proc["StandardError"];
  exit = proc["ExitCode"];
  referenceError = extractError[out];
  pass = exit === 0 &&
    !StringContainsQ[out <> err, "CHESS_RELEASE_EXAMPLE_FAILED"] &&
    If[NumberQ[threshold],
      NumberQ[referenceError] && referenceError < threshold,
      True
    ];
  Print[
    "CHESS_RELEASE_TEST case=", case,
    " exitCode=", exit,
    " maxReferenceError=", formatValue[referenceError],
    " threshold=", formatValue[threshold],
    " status=", If[pass, "PASS", "FAIL"]
  ];
  If[!pass,
    Print["--- stdout ---"];
    Print[out];
    Print["--- stderr ---"];
    Print[err];
  ];
  pass
];

passes = runOne /@ tests;
If[And @@ passes,
  Print["CHESS_RELEASE_TESTS_PASSED"];
  Quit[0],
  Print["CHESS_RELEASE_TESTS_FAILED"];
  Quit[1]
];
