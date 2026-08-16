(* ::Package:: *)

(*
  MatrixAdapters.wl

  The two historical CHESS solvers accept slightly different matrix-function
  conventions:

    canonical:     B[t, "Precision" -> p]
    non-canonical: B[t, p] or B[t]

  The unified front end must not force existing callers to rewrite their matrix
  evaluators.  The helpers below provide the smallest common adapter.  They do
  not cache, simplify, or symbolically inspect user functions; they only try the
  established call conventions and require a numerical square matrix.
*)

ClearAll[
  CHESSNumericMatrixQ,
  CHESSMatrixMatchesDimensionsQ,
  CHESSUnifiedEvaluateMatrix,
  CHESSDimensionCheckedNodeEvaluator,
  CHESSCanonicalNodeEvaluator,
  CHESSCheckedCanonicalPropagate,
  CHESSStructuralZeroCoefficientQ,
  CHESSNormalizePolynomialSpec,
  CHESSFilterOptionRules
];

CHESSUnifiedEvaluateMatrix::matrix =
  "The coefficient evaluator did not return a numerical square matrix at point `1`.";
CHESSUnifiedEvaluateMatrix::dims =
  "The coefficient evaluator returned matrix dimensions `1` at point `2`; expected `3` from the boundary data.";
CHESSNormalizePolynomialSpec::dims =
  "A constant coefficient matrix has dimensions `1`; expected `2`.";

(* A coefficient matrix may be dense or SparseArray, real or complex.  The
   square-dimension check is important because both solvers later multiply the
   matrix by a vector of master integrals without another shape check. *)
CHESSNumericMatrixQ[value_] := Module[{dims},
  If[!MatrixQ[value, NumericQ], Return[False]];
  dims = Dimensions[value];
  Length[dims] === 2 && dims[[1]] === dims[[2]]
];

(* Automatic is useful for standalone adapter calls.  Every public propagation
   route supplies the boundary-derived n x n dimensions, because merely being
   square is not sufficient for a valid differential operator. *)
CHESSMatrixMatchesDimensionsQ[value_, Automatic] := CHESSNumericMatrixQ[value];
CHESSMatrixMatchesDimensionsQ[value_, expectedDimensions_List] :=
  CHESSNumericMatrixQ[value] && Dimensions[value] === expectedDimensions;

(* Evaluate one coefficient using the three conventions already accepted by
   the two source packages.  A constant numerical matrix is also accepted; it
   is useful for tests and for genuinely constant differential equations. *)
CHESSUnifiedEvaluateMatrix[
  spec_, point_, precision_, expectedDimensions_ : Automatic
] := Module[{value},
  If[CHESSNumericMatrixQ[spec],
    value = spec,

    (* Try the conventions lazily.  Building a list of all three attempts
       would evaluate a valid user function three times at every node, which
       is both wasteful and surprising for evaluators with memoization or other
       side effects. *)
    value = Quiet @ Check[
      spec[point, "Precision" -> precision], $Failed
    ];
    If[!CHESSMatrixMatchesDimensionsQ[value, expectedDimensions],
      value = Quiet @ Check[spec[point, precision], $Failed]
    ];
    If[!CHESSMatrixMatchesDimensionsQ[value, expectedDimensions],
      value = Quiet @ Check[spec[point], $Failed]
    ]
  ];
  If[value === $Failed,
    Message[CHESSUnifiedEvaluateMatrix::matrix, point];
    Return[$Failed]
  ];
  If[!CHESSNumericMatrixQ[value],
    Message[CHESSUnifiedEvaluateMatrix::matrix, point];
    Return[$Failed]
  ];
  If[!CHESSMatrixMatchesDimensionsQ[value, expectedDimensions],
    Message[
      CHESSUnifiedEvaluateMatrix::dims,
      Dimensions[value], point, expectedDimensions
    ];
    Return[$Failed]
  ];
  SparseArray @ N[value, precision]
];

(* The canonical core already has a batch-node contract.  Wrapping an explicit
   B_p coefficient in CHESSNodeEvaluator lets the existing collocation code use
   the common evaluator above without modifying a single canonical hot loop.
   The scalar member is retained for completeness; ordinary nodes are evaluated
   by the batch member in one Mathematica call. *)
CHESSDimensionCheckedNodeEvaluator[
  spec_CHESSNodeEvaluator, expectedDimensions_List
] := CHESSNodeEvaluator[
  Function[
    point,
    CHESSUnifiedEvaluateMatrix[
      CHESSNodeEvaluatorScalar[spec], point, 200, expectedDimensions
    ]
  ],
  Function[{points, precision, threads}, Module[{values, checked},
    values = Quiet @ Check[
      CHESSNodeEvaluatorBatch[spec][points, precision, threads],
      $Failed
    ];
    If[values === $Failed || !ListQ[values] ||
        Length[values] =!= Length[points],
      Return[$Failed]
    ];
    checked = MapThread[
      Function[{value, point},
        If[CHESSMatrixMatchesDimensionsQ[value, expectedDimensions],
          SparseArray @ N[value, precision],
          If[CHESSNumericMatrixQ[value],
            Message[
              CHESSUnifiedEvaluateMatrix::dims,
              Dimensions[value], point, expectedDimensions
            ],
            Message[CHESSUnifiedEvaluateMatrix::matrix, point]
          ];
          $Failed
        ]
      ],
      {values, points}
    ];
    If[MemberQ[checked, $Failed], $Failed, checked]
  ]]
];

CHESSCanonicalNodeEvaluator[
  spec_, expectedDimensions_ : Automatic
] := If[
  CHESSNodeEvaluatorQ[spec] && ListQ[expectedDimensions],
  CHESSDimensionCheckedNodeEvaluator[spec, expectedDimensions],
  CHESSNodeEvaluator[
  Function[
    point,
    CHESSUnifiedEvaluateMatrix[spec, point, 200, expectedDimensions]
  ],
  Function[{points, precision, threads}, Module[{values},
    values = CHESSUnifiedEvaluateMatrix[
      spec, #, precision, expectedDimensions
    ] & /@ points;
    If[MemberQ[values, $Failed], $Failed, values]
  ]]
  ]
];

(* The copied canonical core emits SpectralPropagate::solvefail when an epsilon
   layer fails, but its historical Return sits inside Do and can be swallowed by
   that loop.  Keep the source file byte-identical and repair the public route at
   one boundary: any such message makes the whole solve fail closed. *)
CHESSCheckedCanonicalPropagate[
  evaluator_, boundary_, interval_, canonicalRules_List
] := Check[
  SpectralPropagateSequential[
    evaluator, boundary, interval, Sequence @@ canonicalRules
  ],
  $Failed,
  {SpectralPropagate::solvefail}
];

(* Automatic routing must be structural, not tolerance based.  A user function
   which happens to be tiny or happens to vanish at one sampled point is still
   a genuine coefficient.  Only literal zero or an exactly zero constant matrix
   is classified as absent. *)
CHESSStructuralZeroCoefficientQ[spec_] := Module[{rules},
  If[TrueQ[PossibleZeroQ[spec]], Return[True]];
  If[!MatrixQ[spec], Return[False]];
  rules = Most[ArrayRules[SparseArray[spec]]];
  rules === {} || And @@ (TrueQ[PossibleZeroQ[#]] & /@ rules[[All, 2]])
];

(* The old non-canonical solver expects every list entry to be callable using
   f[t,p] or f[t].  Turn zero and constant matrices into two-argument functions,
   and wrap genuine functions with the common evaluator so that the canonical
   f[t,"Precision"->p] convention works on mixed routes as promised by the
   unified interface.  The boundary determines the required n x n dimensions;
   no numerical probing is used for routing. *)
CHESSNormalizePolynomialSpec[specs_List, boundary_?MatrixQ] := Module[
  {n, dims, normalized, entry, normalizationFailure},
  n = Dimensions[boundary][[1]];
  dims = {n, n};
  normalizationFailure = False;
  normalized = Table[
    entry = specs[[p]];
    Which[
      CHESSStructuralZeroCoefficientQ[entry],
        With[{localDims = dims},
          Function[{point, precision}, SparseArray[{}, localDims]]
        ],
      MatrixQ[entry],
        If[Dimensions[entry] =!= dims,
          Message[CHESSNormalizePolynomialSpec::dims, Dimensions[entry], dims];
          normalizationFailure = True;
          $Failed,
          With[{matrix = entry},
            Function[{point, precision}, SparseArray @ N[matrix, precision]]
          ]
        ],
      True,
        With[{coefficient = entry, expectedDimensions = dims},
          Function[{point, precision},
            CHESSUnifiedEvaluateMatrix[
              coefficient, point, precision, expectedDimensions
            ]
          ]
        ]
    ],
    {p, Length[specs]}
  ];
  If[TrueQ[normalizationFailure], Return[$Failed]];
  normalized
];

(* Route only options understood by the selected historical solver.  This
   avoids Mathematica option warnings and, more importantly, makes it explicit
   which layer owns each option. *)
CHESSFilterOptionRules[rules_List, optionDefaults_List] := Module[{names},
  names = First /@ optionDefaults;
  Select[
    rules,
    (MatchQ[Unevaluated[#], _Rule | _RuleDelayed] &&
       MemberQ[names, First[Unevaluated[#]]]) &
  ]
];
