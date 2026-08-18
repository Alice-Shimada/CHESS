(* ::Package:: *)

(*
  NativeMatrixAdapter.wl

  Small conversion helpers for users who already have differential-equation
  matrices in Wolfram Language form.  The output is always the existing
  CHESSNodeEvaluator[scalar,batch] contract; no new backend or dispatcher is
  introduced here.

  Two paths are intentionally separate:

    CHESSNativeMatrixAdapter
      evaluates ordinary coefficient functions or path-variable expressions
      in Mathematica, then hands the sampled matrices to NativeBackend.

    CHESSNativeFlintMatrixAdapter
      evaluates a prepared FORM straight-line program with FLINT and lets a
      problem-specific assembler turn each returned value row into
      {B0,B1,...}.  This path also moves node expression evaluation out of
      Mathematica.
*)

ClearAll[
  CHESSNativeMatrixAdapter,
  CHESSNativeFlintMatrixAdapter,
  CHESSNativeEvaluateCoefficientList,
  CHESSNativeEvaluateExpressionList,
  CHESSNativeValidatePolynomialMatrices,
  CHESSNativeFlintBatchEvaluate
];

CHESSNativeMatrixAdapter::spec =
  "Expected a non-empty coefficient-matrix list and a positive dimension.";
CHESSNativeMatrixAdapter::expr =
  "Every expression coefficient must have dimensions `1`.";
CHESSNativeFlintMatrixAdapter::coord =
  "The coordinate function did not produce one rectangular real row per node.";
CHESSNativeFlintMatrixAdapter::assemble =
  "The FLINT assembly function returned invalid coefficient matrices at node `1`.";
CHESSNativeFlintMatrixAdapter::spec =
  "Expected executable and SLP paths, coordinate and assembly functions, and a positive dimension.";

CHESSNativeValidatePolynomialMatrices[matrices_, dimensions_List] :=
  ListQ[matrices] && matrices =!= {} &&
  And @@ (CHESSMatrixMatchesDimensionsQ[#, dimensions] & /@ matrices);

CHESSNativeEvaluateCoefficientList[
  specs_List, point_, precision_, dimensions_List
] := Module[{values},
  values = CHESSUnifiedEvaluateMatrix[
    #, point, precision, dimensions
  ] & /@ specs;
  If[MemberQ[values, $Failed], $Failed, values]
];

CHESSNativeEvaluateExpressionList[
  expressions_List, variable_Symbol, point_, precision_, dimensions_List
] := Module[{values},
  values = N[expressions /. variable -> point, precision];
  If[
    CHESSNativeValidatePolynomialMatrices[values, dimensions],
    SparseArray /@ values,
    $Failed
  ]
];

(* Convert {B0,B1,...}, where each item is a constant matrix or an evaluator
   accepted by CHESSUnifiedEvaluateMatrix, into the polynomial batch contract
   required by direct mixed Native transport. *)
CHESSNativeMatrixAdapter[
  specs_List, dimension_Integer?Positive
] /; specs =!= {} := With[
  {localSpecs = specs, localDimensions = {dimension, dimension}},
  CHESSNodeEvaluator[
    Function[{point, precision},
      CHESSNativeEvaluateCoefficientList[
        localSpecs, point, precision, localDimensions
      ]
    ],
    (* Keep threads in the shared batch signature; this built-in wrapper is
       serial because arbitrary Wolfram evaluators may not be thread-safe. *)
    Function[{points, precision, threads}, Module[{values},
      values = CHESSNativeEvaluateCoefficientList[
        localSpecs, #, precision, localDimensions
      ] & /@ points;
      If[MemberQ[values, $Failed], $Failed, values]
    ]]
  ]
];

(* Convert explicit matrices containing one path variable.  Dimensions are
   checked before constructing the evaluator, so malformed symbolic input
   fails immediately. *)
CHESSNativeMatrixAdapter[
  expressions_List, variable_Symbol, dimension_Integer?Positive
] /; expressions =!= {} := Module[{dimensions},
  dimensions = {dimension, dimension};
  If[!And @@ (Dimensions[#] === dimensions & /@ expressions),
    Message[CHESSNativeMatrixAdapter::expr, dimensions];
    Return[$Failed]
  ];
  With[
    {
      localExpressions = expressions,
      localVariable = variable,
      localDimensions = dimensions
    },
    CHESSNodeEvaluator[
      Function[{point, precision},
        CHESSNativeEvaluateExpressionList[
          localExpressions, localVariable, point, precision, localDimensions
        ]
      ],
      (* Expression substitution is intentionally serial for the same contract
         and thread-safety reasons as the coefficient-list adapter above. *)
      Function[{points, precision, threads}, Module[{values},
        values = CHESSNativeEvaluateExpressionList[
          localExpressions, localVariable, #, precision, localDimensions
        ] & /@ points;
        If[MemberQ[values, $Failed], $Failed, values]
      ]]
    ]
  ]
];

CHESSNativeMatrixAdapter[___] := (
  Message[CHESSNativeMatrixAdapter::spec];
  $Failed
);

Options[CHESSNativeFlintMatrixAdapter] = Options[CHESSFlintRunSLP];

CHESSNativeFlintBatchEvaluate[
  executable_String,
  slpFile_String,
  coordinateFunction_,
  assembleFunction_,
  dimensions_List,
  flintOptions_List,
  points_List,
  precision_,
  threads_
] := Module[{coordinates, result, assembled, badPoint},
  coordinates = Quiet @ Check[
    coordinateFunction[#, precision] & /@ points,
    $Failed
  ];
  If[
    coordinates === $Failed || !MatrixQ[coordinates, NumericQ] ||
    !FreeQ[coordinates, _Complex] || Length[coordinates] =!= Length[points],
    Message[CHESSNativeFlintMatrixAdapter::coord];
    Return[$Failed]
  ];
  result = CHESSFlintRunSLP[
    executable, slpFile, coordinates, precision, threads,
    Sequence @@ flintOptions
  ];
  If[!AssociationQ[result], Return[$Failed]];
  assembled = MapThread[
    Quiet @ Check[
      assembleFunction[#1, #2, precision],
      $Failed
    ] &,
    {points, result["Values"]}
  ];
  If[Length[assembled] =!= Length[points],
    Message[CHESSNativeFlintMatrixAdapter::assemble, Missing["Unknown"]];
    Return[$Failed]
  ];
  badPoint = FirstCase[
    MapThread[
      If[
        CHESSNativeValidatePolynomialMatrices[#1, dimensions],
        Nothing,
        #2
      ] &,
      {assembled, points}
    ],
    _, Missing["Absent"]
  ];
  If[!MissingQ[badPoint],
    Message[CHESSNativeFlintMatrixAdapter::assemble, badPoint];
    Return[$Failed]
  ];
  Map[SparseArray /@ # &, assembled]
];

(* coordinateFunction[point,precision] returns one real SLP input row.
   assembleFunction[point,values,precision] returns {B0,B1,...}. *)
CHESSNativeFlintMatrixAdapter[
  executable_String,
  slpFile_String,
  coordinateFunction_,
  assembleFunction_,
  dimension_Integer?Positive,
  opts : OptionsPattern[]
] := Module[{dimensions, flintOptions, batch},
  dimensions = {dimension, dimension};
  flintOptions = FilterRules[{opts}, Options[CHESSFlintRunSLP]];
  batch = Function[{points, precision, threads},
    CHESSNativeFlintBatchEvaluate[
      executable, slpFile, coordinateFunction, assembleFunction, dimensions,
      flintOptions, points, precision, threads
    ]
  ];
  CHESSNodeEvaluator[
    Function[{point, precision}, Module[{value},
      value = batch[{point}, precision, 1];
      If[ListQ[value] && Length[value] === 1, First[value], $Failed]
    ]],
    batch
  ]
];

CHESSNativeFlintMatrixAdapter[___] := (
  Message[CHESSNativeFlintMatrixAdapter::spec];
  $Failed
);
