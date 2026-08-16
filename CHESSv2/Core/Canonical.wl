(* ::Package:: *)

(* Chess: release-ready Chebyshev pseudospectral transport package.
   The caller supplies dLettersLine, Atilde, and boundary values for each case.
   Atilde must use only the unified letter markers Log[W[i]] or logW[i]. *)

(*
  Release implementation policy:
  - use ordinary positional lists instead of keyed containers;
  - keep matrix data sparse whenever possible;
  - keep case-specific input conversions in the examples, not in the core.
  The result is less decorative, but faster to distribute to subkernels and
  cheaper to index in the hot loops.
*)

ClearAll[
  nAt,
  CHESSNormalizeAtildeInput,
  CHESSBuildAtildeLinearData,
  CHESSAtildeLinearData,
  CHESSAtildeLinearCombination,
  $CHESSAtildeLinearHash,
  $CHESSAtildeLinearData,
  $CHESSCacheGeneration,
  $CHESSVersion,
  CHESSNodeEvaluator,
  CHESSNodeEvaluatorQ,
  CHESSNodeEvaluatorScalar,
  CHESSNodeEvaluatorBatch,
  CHESSBatchThreadCount
];

$CHESSVersion = "0.2.0";
$CHESSCacheGeneration = 0;

(* Public messages for malformed prepared input. *)
nAt::nodata = "Set global dLettersLine and Atilde, or a cached $CHESSAtildeLinearData, before calling nAt.";
CHESSBuildAtildeLinearData::badlog =
  "Atilde contains logarithms outside the unified Log[W[i]]/logW[i] interface. Convert case-specific letter formats before building the package data.";
CHESSAtildeLinearCombination::slots =
  "dLettersLine has length `1`, but Atilde uses logW slot `2`.";

nAt::usage = "nAt[t] evaluates the cached one-dimensional differential-equation matrix at path parameter t, using global Atilde and dLettersLine.";
SpectralPropagate::usage = "SpectralPropagate[nAfun, y0, {x0, x1}, opts] transports epsilon coefficients from x0 to x1 by Chebyshev--Lobatto collocation.";
CHESSFinalState::usage = "CHESSFinalState[result] returns the transported epsilon-coefficient matrix at the right endpoint.";
CHESSEndpointInfo::usage = "CHESSEndpointInfo[result] returns endpoint regularization diagnostics stored in a CHESS result.";
CHESSClearCaches::usage = "CHESSClearCaches[] clears cached Atilde decompositions, endpoint data, and node operators.";
CHESSAtildeLinearData::usage = "CHESSAtildeLinearData[] returns the cached sparse letter-slot decomposition {zero, slots, matrices} of the current Atilde.";
CHESSAtildeLinearCombination::usage = "CHESSAtildeLinearCombination[values, prec] combines cached Atilde letter matrices with dLog values.";
CHESSNodeEvaluator::usage = "CHESSNodeEvaluator[scalar,batch] supplies both a scalar matrix evaluator and batch[nodes,precision,threads] to CHESS.";

(* A positional wrapper advertises a batch node backend without changing the
   long-standing scalar evaluator interface or introducing keyed containers. *)
CHESSNodeEvaluatorQ[expr_] := MatchQ[Unevaluated[expr], CHESSNodeEvaluator[_, _]];
CHESSNodeEvaluatorScalar[CHESSNodeEvaluator[scalar_, _]] := scalar;
CHESSNodeEvaluatorBatch[CHESSNodeEvaluator[_, batch_]] := batch;

(* Native evaluators manage their own threads and therefore need the requested
   count without launching Wolfram subkernels at this stage. *)
CHESSBatchThreadCount[parallelSpec_, m_Integer, kernelSpec_] := Which[
  parallelSpec === False || parallelSpec === None, 1,
  parallelSpec === Automatic && m < 8, 1,
  True, CHESSRequestedKernelCount[kernelSpec, m]
];

(* Normalize the public Atilde input convention to internal logW[i] markers. *)
CHESSNormalizeAtildeInput[mat_] :=
  mat /. Log[W[i_Integer]] :> logW[i];

(* Split Atilde into sparse matrices multiplying logW slots.  Terms independent
   of logW are constants in Atilde and disappear after differentiating. *)
CHESSBuildAtildeLinearData[mat_] := Module[
  {normalized, sparse, dims, rules, harvested, slotEntries, slots, matrices},

  (* The package has one public letter convention.  Any family-specific heads
     must be converted by the caller before this decomposition is built. *)
  normalized = CHESSNormalizeAtildeInput[mat];
  If[!FreeQ[normalized, _Log],
    Message[CHESSBuildAtildeLinearData::badlog];
    Return[$Failed];
  ];

  (* Work only with explicit nonzero entries.  This avoids scanning the many
     structural zeros in the differential-equation matrix. *)
  sparse = SparseArray[normalized];
  dims = Dimensions[sparse];
  rules = Most[ArrayRules[sparse]];

  (* For each nonzero entry, extract only coefficients of logW slots.  Reap/Sow
     groups entries by slot without keyed containers. *)
  harvested = Reap[
    Do[
      Module[{pos, expr, exprVars, coeff},
        pos = rules[[r, 1]];
        expr = rules[[r, 2]];

        (* Keep the original logW symbol together with its integer slot so that
           Coefficient can extract the exact multiplier. *)
        exprVars = DeleteDuplicates @ Cases[
          expr,
          lw : logW[i_] :> {Round[i], lw},
          {0, Infinity}
        ];

        (* Every letter slot has its own sparse coefficient matrix.  Pure
           constants have no logW slot and contribute nothing to B(t). *)
        Do[
          coeff = Coefficient[expr, exprVars[[slot, 2]]];
          If[coeff =!= 0, Sow[pos -> coeff, exprVars[[slot, 1]]]],
          {slot, Length[exprVars]}
        ];
      ],
      {r, Length[rules]}
    ],
    _,
    List
  ][[2]];
  slotEntries = harvested;
  slots = If[slotEntries === {}, {}, slotEntries[[All, 1]]];

  (* The result is positional: {zero matrix, letter slots, slot matrices}.
     Keeping this as plain lists avoids keyed-container overhead in the hot path. *)
  matrices = If[slotEntries === {}, {}, SparseArray[Flatten[#[[2]], 1], dims] & /@ slotEntries];
  {SparseArray[{}, dims], slots, matrices}
];

(* Cache the sparse linear decomposition of the current global Atilde. *)
CHESSAtildeLinearData[] := Module[{hash},

  (* Some callers may distribute only the prebuilt linear data to subkernels. *)
  If[!ValueQ[Atilde] && ValueQ[$CHESSAtildeLinearData],
    Return[$CHESSAtildeLinearData];
  ];
  If[!ValueQ[Atilde], Return[$Failed]];

  (* Rebuild only when the global Atilde expression actually changes. *)
  hash = Hash[Atilde];
  If[
    ValueQ[$CHESSAtildeLinearHash] && ValueQ[$CHESSAtildeLinearData] &&
      $CHESSAtildeLinearHash === hash,
    Return[$CHESSAtildeLinearData];
  ];
  $CHESSAtildeLinearHash = hash;
  $CHESSAtildeLinearData = CHESSBuildAtildeLinearData[Atilde]
];

Options[CHESSAtildeLinearCombination] = {};

(* Evaluate the cached Atilde decomposition after dLog letter values are known. *)
CHESSAtildeLinearCombination[values_, prec_, OptionsPattern[]] := Module[
  {data, slots, matrices, acc, i, maxSlot, dims},
  data = CHESSAtildeLinearData[];
  If[data === $Failed, Return[$Failed]];

  dims = Dimensions[data[[1]]];
  slots = data[[2]];
  matrices = data[[3]];
  If[
    slots === {},
    SparseArray[{}, dims],
    maxSlot = Max[slots];

    (* A short dLettersLine would silently pick wrong data via Part errors; fail
       early with the largest required slot. *)
    If[Length[values] < maxSlot,
      Message[CHESSAtildeLinearCombination::slots, Length[values], maxSlot];
      Return[$Failed];
    ];

    (* Accumulate in place rather than materializing a list of sparse matrices
       at every node evaluation. *)
    acc = SparseArray[{}, dims];
    Do[
      acc = acc + N[values[[slots[[i]]]], prec] matrices[[i]],
      {i, Length[slots]}
    ];
    N[acc, prec]
  ]
];

Options[nAt] = {"Precision" -> 200};

(* Numerical differential matrix B(t) = sum_i dLog(W_i(t)) A_i. *)
nAt[tvalue_, OptionsPattern[]] := Module[{ndLettersLine, p},

  (* nAt is intentionally thin: case-specific letter preparation belongs to the
     notebook/script, while this package only substitutes t and combines slots. *)
  If[!ValueQ[dLettersLine] || (!ValueQ[Atilde] && !ValueQ[$CHESSAtildeLinearData]),
    Message[nAt::nodata];
    Return[$Failed];
  ];
  p = OptionValue["Precision"];

  (* Evaluate every prepared dLog letter at the requested precision before
     forming the sparse matrix combination. *)
  ndLettersLine = N[dLettersLine /. t -> N[tvalue, p], p];
  CHESSAtildeLinearCombination[ndLettersLine, p]
];


(* ::Section:: *)
(* Chebyshev pseudospectral propagation *)


ClearAll[
  ChebyshevLobattoData,
  BaseScalarCollocationMatrix,
  LiftedScalarCollocationMatrix,
  SpectralPropagate,
  SpectralPropagateSequential,
  CHESSResultPart,
  CHESSNodes,
  CHESSStateVectors,
  CHESSFinalState,
  CHESSEndpointSubBlocks,
  CHESSEndpointInfo,
  EndpointDirection,
  EndpointShiftedExpression,
  FiniteNumericValueQ,
  OneSidedLaurentData,
  FastEndpointLetterData,
  CHESSLettersHash,
  CHESSAtildeDataHash,
  EndpointLetterPoleData,
  EndpointLetterPoleDataCached,
  EndpointMatrixPoleData,
  EndpointMatrixPoleDataCached,
  RegularizedEndpointOperatorData,
  RegularizedEndpointOperatorDataCached,
  EndpointRegularityResidual,
  NormalizeRegularizedEndpoints,
  ActiveRegularizedEndpoints,
  LeftEndpointRegularizedQ,
  EndpointPointIndex,
  CHESSParallelAvailableQ,
  CHESSResolveParallelQ,
  CHESSPrepareParallel,
  CHESSPrepareEndpointParallel,
  CHESSPrepareNodeParallel,
  BuildNodeOperatorData,
  EndpointSequentialSource,
  NodeSequentialSource,
  CHESSMaxLocalKernelCount,
  CHESSRequestedKernelCount,
  CHESSValidSolutionMatrixQ,
  CHESSScalarCollocationSolve,
  CHESSBuildLinearSolver,
  CHESSClearCaches
];

SpectralPropagate::solvefail =
  "Linear solve failed at epsilon layer `1`; returned head `2` with dimensions `3`.";
SpectralPropagate::badbatch =
  "The batch node evaluator returned `1`; expected one matrix for each of `2` ordinary nodes.";

SpectralPropagate::usage = "SpectralPropagate[nAfun, y0, {x0, x1}, opts] transports epsilon coefficients from x0 to x1 by Chebyshev--Lobatto collocation.";
CHESSFinalState::usage = "CHESSFinalState[result] returns the transported epsilon-coefficient matrix at the right endpoint.";
CHESSEndpointInfo::usage = "CHESSEndpointInfo[result] returns endpoint regularization diagnostics stored in a CHESS result.";
CHESSClearCaches::usage = "CHESSClearCaches[] clears cached Atilde decompositions, endpoint data, and node operators.";

(* Result accessors keep the result as a compact positional list without keyed data. *)
CHESSResultPart[result_, "Nodes"] := result[[1]];
CHESSResultPart[result_, "StateVectors"] := result[[2]];
CHESSResultPart[result_, "FinalState"] := result[[3]];
CHESSResultPart[result_, "EndpointSubBlocks"] := result[[4]];
CHESSResultPart[result_, "EndpointInfo"] := result[[5]];
CHESSResultPart[result_, "SolveStrategy"] := result[[6]];
CHESSResultPart[result_, _] := $Failed;

CHESSNodes[result_] := result[[1]];
CHESSStateVectors[result_] := result[[2]];
CHESSFinalState[result_] := result[[3]];
CHESSEndpointSubBlocks[result_] := result[[4]];
CHESSEndpointInfo[result_] := result[[5]];

(*
  Spurious endpoint pole handling.

  Near a regularized endpoint c the pulled-back operator is expanded as
      B(t) = R/(t-c) + Mfinite + O(t-c).
  The active sequential epsilon solver stores only the sparse pair
      {R, Mfinite}.
  Finite endpoint sources are built weight-by-weight by applying Mfinite to
  known solution vectors first and then applying residue actions.
*)

EndpointDirection::notendpoint = "Point `1` is not an endpoint of interval `2`.";

(* Return the one-sided direction used to expand around an interval endpoint. *)
EndpointDirection[a_, {x0_, x1_}] := Which[
  (* The local coordinate s is always positive; left and right endpoints differ
     only by whether t=a+s or t=a-s. *)
  PossibleZeroQ[N[a - x0, 200]], "FromAbove",
  PossibleZeroQ[N[a - x1, 200]], "FromBelow",
  True,
    Message[EndpointDirection::notendpoint, a, {x0, x1}];
    $Failed
];

(* Recenter an expression at the endpoint using a positive local coordinate s. *)
EndpointShiftedExpression[expr_, a_, dir_, s_Symbol] := Switch[
  (* FromBelow uses a minus sign, which later flips the residue back to the
     derivative with respect to the global coordinate t. *)
  dir,
  "FromAbove", expr /. t -> a + s,
  "FromBelow", expr /. t -> a - s
];

(* True only for finite numerical endpoint values. *)
FiniteNumericValueQ[x_] :=
  NumericQ[x] && FreeQ[x, _DirectedInfinity | ComplexInfinity | Indeterminate];

(* Extract residue and finite part of dLog(W) at a one-sided endpoint.
   SeriesCoefficient is used directly because it is much faster than building
   full series for the letters encountered in this workflow. *)
OneSidedLaurentData[expr_, a_, interval : {x0_, x1_}, prec_ : 200] := Module[
  {dir, s, shifted, residueS, residue, regular, signFactor, cleanLocalCoefficient},
  dir = EndpointDirection[a, interval];
  s = Unique["s"];

  (* Work in a one-sided local coordinate.  This keeps branch assumptions local
     and avoids asking Simplify to reason about the full transport interval. *)
  shifted = EndpointShiftedExpression[expr, a, dir, s];
  signFactor = If[dir === "FromAbove", 1, -1];

  (* The local coordinate is one-sided and positive.  Most coefficients are
     already free of s; only simplify the rare ones that retain branch forms
     such as s/Sqrt[s^2]. *)
  cleanLocalCoefficient[x_] := If[
    FreeQ[x, s],
    x,
    Quiet @ Assuming[s > 0, Simplify[x]]
  ];

  residueS = cleanLocalCoefficient @ Quiet @ Check[SeriesCoefficient[shifted, {s, 0, -1}], 0];
  regular = cleanLocalCoefficient @ Quiet @ Check[SeriesCoefficient[shifted, {s, 0, 0}], 0];

  (* SeriesCoefficient was taken in s.  Convert the residue to the global t
     convention so B(t)=R/(t-a)+Mfinite+... holds for both endpoints. *)
  residue = signFactor residueS;

  {N[residue, prec], N[regular, prec], "SeriesCoefficient"}
];

(* Fast endpoint letter extraction: use direct substitution for regular letters
   and fall back to Laurent coefficients only when a pole is present. *)
FastEndpointLetterData[expr_, a_, interval : {x0_, x1_}, prec_ : 200] := Module[
  {direct},

  (* Most letters are regular at a chosen endpoint.  Direct substitution is far
     cheaper than a Laurent expansion, so use SeriesCoefficient only on poles. *)
  direct = Quiet @ Check[N[expr /. t -> N[a, prec], prec], $Failed];
  If[
    direct =!= $Failed && FiniteNumericValueQ[direct],
    {0, direct, "Direct"},
    OneSidedLaurentData[expr, a, interval, prec]
  ]
];

(* Hash helpers for memoized endpoint data.  They make it safe to run several
   families in the same kernel without reusing residue data from a previous
   dLettersLine/Atilde pair. *)
CHESSLettersHash[] :=
  If[ValueQ[dLettersLine], Hash[dLettersLine], 0];

(* Build/update the cached Atilde decomposition before deriving the data key. *)
CHESSAtildeDataHash[] := Module[{data},
  data = CHESSAtildeLinearData[];
  If[data === $Failed, 0, If[ValueQ[$CHESSAtildeLinearHash], $CHESSAtildeLinearHash, Hash[data]]]
];

EndpointLetterPoleData[a_, interval : {x0_, x1_}, prec_ : 200, parallelQ_ : False] :=
  EndpointLetterPoleDataCached[$CHESSCacheGeneration, CHESSLettersHash[], a, interval, prec, parallelQ];

(* Compute endpoint residue/finite vectors for every letter in dLettersLine. *)
EndpointLetterPoleDataCached[generation_, lettersHash_, a_, interval : {x0_, x1_}, prec_, parallelQ_] :=
  EndpointLetterPoleDataCached[generation, lettersHash, a, interval, prec, parallelQ] = Module[
    {data, methods},
    (* Each letter is independent, so this is a safe parallel stage. *)
    data = If[
      TrueQ[parallelQ],
      Quiet @ ParallelMap[
        FastEndpointLetterData[#, a, interval, prec] &,
        dLettersLine,
        Method -> "FinestGrained"
      ],
      FastEndpointLetterData[#, a, interval, prec] & /@ dLettersLine
    ];
    methods = data[[All, 3]];

    (* Return residue vector, finite vector, and diagnostic method counts. *)
    {data[[All, 1]], data[[All, 2]], Count[methods, "Direct"], Count[methods, "SeriesCoefficient"]}
  ];

(* Convert letter endpoint data into residue and finite matrices for A(t). *)
EndpointMatrixPoleData[a_, interval : {x0_, x1_}, prec_ : 200, parallelQ_ : False] :=
  EndpointMatrixPoleDataCached[$CHESSCacheGeneration, CHESSLettersHash[], CHESSAtildeDataHash[], a, interval, prec, parallelQ];

EndpointMatrixPoleDataCached[generation_, lettersHash_, atildeHash_, a_, interval : {x0_, x1_}, prec_, parallelQ_] :=
  EndpointMatrixPoleDataCached[generation, lettersHash, atildeHash, a, interval, prec, parallelQ] = Module[
    {letterData},
    letterData = EndpointLetterPoleData[a, interval, prec, parallelQ];
    {
      (* Feed the residue/finite letter vectors through the same Atilde linear
         combination used at ordinary nodes. *)
      N[CHESSAtildeLinearCombination[letterData[[1]], prec], prec],
      N[CHESSAtildeLinearCombination[letterData[[2]], prec], prec],
      letterData[[3]],
      letterData[[4]]
    }
  ];

(* Current sequential-weight propagator only needs the unstacked residue/finite
   pair {R, Mfinite}; higher-weight endpoint source terms are built on demand. *)
RegularizedEndpointOperatorData[a_, interval : {x0_, x1_}, prec_ : 200, parallelQ_ : False] :=
  RegularizedEndpointOperatorDataCached[$CHESSCacheGeneration, CHESSLettersHash[], CHESSAtildeDataHash[], a, interval, prec, parallelQ];

RegularizedEndpointOperatorDataCached[generation_, lettersHash_, atildeHash_, a_, interval : {x0_, x1_}, prec_, parallelQ_] :=
  RegularizedEndpointOperatorDataCached[generation, lettersHash, atildeHash, a, interval, prec, parallelQ] = Module[
    {opData},
    opData = EndpointMatrixPoleData[a, interval, prec, parallelQ];
    {
      (* Store endpoint operators sparsely; the source recursion acts on vectors
         and should not form dense matrix powers. *)
      SparseArray @ N[opData[[1]], prec],
      SparseArray @ N[opData[[2]], prec]
    }
  ];

(* Start a fresh problem without reloading the package.  Existing memoized
   endpoint data is left in old generations and will no longer be reused. *)
CHESSClearCaches[] := (
  Clear[$CHESSAtildeLinearHash, $CHESSAtildeLinearData];
  $CHESSCacheGeneration = If[ValueQ[$CHESSCacheGeneration], $CHESSCacheGeneration + 1, 1];
  Null
);

(* Check the regularity condition R.y(a)=0 for supplied endpoint values. *)
EndpointRegularityResidual[a_, yAtEndpoint_, interval : {x0_, x1_}, prec_ : 200] := Module[
  {residue},
  (* This is a diagnostic helper: for a genuine spurious pole, boundary data
     should lie close to the kernel of the residue matrix. *)
  residue = EndpointMatrixPoleData[a, interval, prec][[1]];
  Max @ Abs @ Flatten @ N[residue . yAtEndpoint, prec]
];

(* Normalize endpoint options such as "Left", "Right", All, and None. *)
NormalizeRegularizedEndpoints[spec_, {x0_, x1_}] := Module[
  {raw},
  (* The public option accepts symbolic names and direct endpoint coordinates. *)
  raw = Which[
    spec === None || spec === False || spec === {},
      {},
    spec === All,
      {x0, x1},
    True,
      Flatten @ {spec}
  ];

  (* Replace names after flattening, then remove duplicates while preserving the
     caller's order. *)
  DeleteDuplicates @ Replace[raw, {"Left" -> x0, "Right" -> x1}, {1}]
];

(* Endpoints active in ordinary collocation exclude the left boundary row. *)
ActiveRegularizedEndpoints[endpoints_, nodes_] :=
  Select[endpoints, EndpointPointIndex[#, Rest[nodes]] > 0 &];

(* The left endpoint needs special lifted unknowns z=(y-y0)/(t-x0). *)
LeftEndpointRegularizedQ[endpoints_, x0_] := EndpointPointIndex[x0, endpoints] > 0;

(* Return the 1-based position of x in an endpoint list, or 0 if absent. *)
EndpointPointIndex[x_, endpoints_] := Module[{pos},
  (* PossibleZeroQ keeps exact endpoints and high-precision numeric endpoints
     comparable without forcing machine precision. *)
  pos = FirstPosition[PossibleZeroQ[x - #] & /@ endpoints, True, {0}];
  pos[[1]]
];

(* Chebyshev-Lobatto nodes and derivative matrix mapped to {a,b}.
   Nodes are returned in increasing t order, so the first node is the left endpoint. *)
ChebyshevLobattoData[m_Integer?Positive, {a_, b_}, prec_ : 100] := Module[
  {xCheb, c, dX, xGrid, dXMat, dTheta, xMapped, order, dMapped},
  (* Start with the standard descending Lobatto grid on [-1,1]. *)
  xCheb = N[Table[Cos[Pi j/m], {j, 0, m}], prec];
  c = N[Join[{2}, ConstantArray[1, m - 1], {2}] * (-1)^Range[0, m], prec];

  (* Barycentric differentiation matrix on the reference interval. *)
  xGrid = ConstantArray[xCheb, m + 1];
  dX = xGrid - Transpose[xGrid];
  dTheta = (Outer[Times, c, 1/c])/(dX + IdentityMatrix[m + 1]);
  dXMat = dTheta - DiagonalMatrix[Total[dTheta, {2}]];

  (* Map to {a,b}.  The minus sign appears because the construction above uses
     the descending reference order before the final reversal. *)
  xMapped = N[(a + b)/2 + (b - a)/2 xCheb, prec];
  dMapped = N[-(2/(b - a)) dXMat, prec];

  (* Return nodes in increasing transport direction, from left to right. *)
  order = Reverse[Range[m + 1]];
  {xMapped[[order]], dMapped[[order, order]]}
];

(* Scalar collocation matrix for y'(t)=source(t), with y(x0) fixed in row 1. *)
BaseScalarCollocationMatrix[derivMat_] := Module[{m, rules},
  m = Length[derivMat] - 1;

  (* Row 1 is the boundary condition y(x0)=given.  All other rows are the
     derivative collocation equations. *)
  rules = Join[
    {{1, 1} -> 1},
    Flatten[
      Table[
        {j, l} -> derivMat[[j, l]],
        {j, 2, m + 1},
        {l, 1, m + 1}
      ],
      1
    ]
  ];
  SparseArray[rules, {m + 1, m + 1}]
];

(* Scalar collocation matrix for left endpoint lifting:
      y(t) = y0 + (t-x0) z(t), so y' = z + (t-x0) z'. *)
LiftedScalarCollocationMatrix[derivMat_, nodes_, x0_] :=
  (* This matrix acts on z-values at all nodes, including the left endpoint. *)
  SparseArray[IdentityMatrix[Length[nodes]] + DiagonalMatrix[nodes - x0] . derivMat];

(* Maximum local Wolfram subkernels configured for this machine. *)
CHESSMaxLocalKernelCount[] := CHESSMaxLocalKernelCount[] = Module[
  {counts, name},
  (* Prefer Mathematica's configured local kernel cap when available. *)
  counts = Quiet @ Check[
    Table[
      name = ToLowerCase[ToString[Quiet @ Check[($ConfiguredKernels[[i]])["Name"], ""]]];
      If[
        name == "local" || name == "localhost",
        Quiet @ Check[($ConfiguredKernels[[i]])["KernelCount"], Nothing],
        Nothing
      ],
      {i, Length[$ConfiguredKernels]}
    ],
    {}
  ];
  counts = Cases[counts, _Integer?Positive];

  (* Fall back to processor count on systems without an explicit configuration. *)
  If[counts === {}, Max[1, $ProcessorCount], Max[counts]]
];

(* Interpret "ParallelKernels": Automatic defaults to one; All/Infinity uses the local cap. *)
CHESSRequestedKernelCount[spec_, m_Integer] := Module[
  {maxKernels, requested},
  maxKernels = CHESSMaxLocalKernelCount[];

  (* Automatic stays conservative for reproducibility; All/Infinity means the
     largest local subkernel count Mathematica exposes. *)
  requested = Which[
    spec === Automatic, 1,
    spec === All || spec === Infinity, maxKernels,
    IntegerQ[spec], spec,
    True, 1
  ];
  Min[maxKernels, Max[1, requested]]
];

(* Launch enough local subkernels for parallel node/endpoint preparation. *)
CHESSParallelAvailableQ[kernelSpec_, m_Integer] := Module[{target},
  target = CHESSRequestedKernelCount[kernelSpec, m];
  If[target <= 1, Return[False]];

  (* Launch only the missing kernels; reuse already-open kernels when possible. *)
  Quiet @ Check[
    If[Length[Kernels[]] < target, LaunchKernels[target - Length[Kernels[]]]];
    Length[Kernels[]] > 0,
    False
  ]
];

(* Decide whether the current problem is large enough to benefit from parallel setup. *)
CHESSResolveParallelQ[spec_, m_Integer, kernelSpec_] := Which[
  (* Small node counts often lose to launch/distribution overhead. *)
  spec === False || spec === None, False,
  spec === Automatic, m >= 8 && CHESSParallelAvailableQ[kernelSpec, m],
  True, CHESSParallelAvailableQ[kernelSpec, m]
];

(* Send endpoint helper definitions to subkernels before ParallelMap. *)
CHESSPrepareEndpointParallel[] := Module[{dir},
  dir = Directory[];

  (* Keep relative file behavior identical on subkernels for user-defined
     letter expressions that may contain Get/Import side effects. *)
  Quiet @ Check[ParallelEvaluate[SetDirectory[dir]], Null];

  (* Only endpoint extraction helpers are needed for ParallelMap over letters. *)
  Quiet @ Check[
    DistributeDefinitions[
      dLettersLine,
      t,
      EndpointDirection,
      EndpointShiftedExpression,
      FiniteNumericValueQ,
      OneSidedLaurentData,
      FastEndpointLetterData,
      EndpointPointIndex
    ],
    Null
  ];
];

(* Send nAt and its cached Atilde decomposition to subkernels before ParallelTable. *)
CHESSPrepareNodeParallel[nAfun_] := Module[{dir},
  dir = Directory[];

  (* Build the linear Atilde data once on the master before distributing it. *)
  If[ValueQ[Atilde], CHESSAtildeLinearData[]];
  Quiet @ Check[ParallelEvaluate[SetDirectory[dir]], Null];

  (* Node evaluation needs both dLettersLine and the cached Atilde decomposition. *)
  Quiet @ Check[
    DistributeDefinitions[
      dLettersLine,
      t,
      logW,
      nAt,
      CHESSNormalizeAtildeInput,
      CHESSBuildAtildeLinearData,
      CHESSAtildeLinearData,
      CHESSAtildeLinearCombination,
      $CHESSAtildeLinearHash,
      $CHESSAtildeLinearData
    ],
    Null
  ];
  Quiet @ Check[DistributeDefinitions[nAfun], Null];
];

(* Prepare both endpoint and node evaluation helpers. *)
CHESSPrepareParallel[nAfun_] := (
  CHESSPrepareEndpointParallel[];
  CHESSPrepareNodeParallel[nAfun]
);

(* Build per-node operators.
   A regularized endpoint stores {1,R,Mfinite}; an ordinary node stores {0,B(t_j)}. *)
BuildNodeOperatorData[
  nAfun_, nodes_, endpoints_, endpointOperators_, prec_, precA_, parallelQ_,
  firstIndex_ : 2, batchThreads_ : 1
] := Module[
  {
    m, buildOne, taskCount, batchQ, indices, ordinaryIndices,
    ordinaryValues, ordinaryCursor, result, endpointIndex, position
  },
  m = Length[nodes] - 1;
  taskCount = m + 2 - firstIndex;
  batchQ = CHESSNodeEvaluatorQ[nAfun];

  (* A batch evaluator receives all ordinary nodes in one call.  Endpoint rows
     remain owned by the existing regularization code and are spliced back into
     their original positions afterwards. *)
  If[batchQ,
    indices = Range[firstIndex, m + 1];
    ordinaryIndices = Select[
      indices,
      EndpointPointIndex[nodes[[#]], endpoints] === 0 &
    ];
    ordinaryValues = If[
      ordinaryIndices === {},
      {},
      Quiet @ Check[
        CHESSNodeEvaluatorBatch[nAfun][
          nodes[[ordinaryIndices]],
          precA,
          batchThreads
        ],
        $Failed
      ]
    ];
    If[ordinaryValues === $Failed || !ListQ[ordinaryValues] ||
        Length[ordinaryValues] =!= Length[ordinaryIndices],
      Message[
        SpectralPropagate::badbatch,
        If[ordinaryValues === $Failed, $Failed, Head[ordinaryValues]],
        Length[ordinaryIndices]
      ];
      Return[$Failed]
    ];

    result = ConstantArray[Null, Length[indices]];
    ordinaryCursor = 1;
    Do[
      position = indices[[j]];
      endpointIndex = EndpointPointIndex[nodes[[position]], endpoints];
      If[endpointIndex > 0,
        result[[j]] = {
          1,
          endpointOperators[[endpointIndex, 1]],
          endpointOperators[[endpointIndex, 2]]
        },
        result[[j]] = {
          0,
          SparseArray[N[ordinaryValues[[ordinaryCursor]], prec]]
        };
        ordinaryCursor++
      ],
      {j, Length[indices]}
    ];
    Return[result]
  ];

  (* The left boundary row is skipped in ordinary propagation but included for
     the lifted left-regularized problem. *)
  buildOne[j_] := Module[{endpointIndex},
    endpointIndex = EndpointPointIndex[nodes[[j]], endpoints];
    If[
      endpointIndex > 0,
      (* Endpoint nodes use the finite regularized pair {R,Mfinite}. *)
      {1, endpointOperators[[endpointIndex, 1]], endpointOperators[[endpointIndex, 2]]},
      (* Ordinary nodes use the numerical differential matrix B(t_j). *)
      {0, SparseArray[N[nAfun[nodes[[j]], "Precision" -> precA], prec]]}
    ]
  ];

  (* CoarsestGrained keeps each subkernel working on a block of nodes, reducing
     scheduling overhead for expensive nAt evaluations. *)
  If[
    TrueQ[parallelQ] && taskCount >= 8,
    Quiet @ ParallelTable[buildOne[j], {j, firstIndex, m + 1}, Method -> "CoarsestGrained"],
    Table[buildOne[j], {j, firstIndex, m + 1}]
  ]
];

(* Source term at a regularized endpoint for the current epsilon coefficient.
   Always form Mfinite.y first, then apply residue powers to the resulting vector.
   This ordering is intentional: in these systems Mfinite.y can acquire extra zeros
   that would be hidden or delayed by forming matrix products such as R.Mfinite. *)
EndpointSequentialSource[residue_, regular_, coeffSolutions_, coeff_Integer, nodeIndex_Integer, n_Integer] := Module[
  {source, offset, y, regularAction, v, repeat},
  source = ConstantArray[0, n];
  Do[
    (* For offset r, use y_{coeff-r} at the same endpoint node. *)
    y = coeffSolutions[[coeff - offset, nodeIndex]];

    (* Parentheses are intentional: compute Mfinite.y before any residue action so
       structural zeros in the vector are exposed as early as possible. *)
    regularAction = regular . (y);
    v = regularAction;

    (* Apply R^(offset-1) by repeated sparse matrix-vector products.  Do not
       precompute R powers; they can become less sparse. *)
    Do[
      v = residue . (v),
      {repeat, 1, offset - 1}
    ];
    source = source + v;,
    {offset, 1, coeff - 1}
  ];
  source
];

(* Source term at either an ordinary node or a regularized endpoint. *)
NodeSequentialSource[nodeData_, coeffSolutions_, coeff_Integer, nodeIndex_Integer, n_Integer] := If[
  nodeData[[1]] === 1,
  (* Regularized endpoint: finite source built from residue/regular pieces. *)
  EndpointSequentialSource[nodeData[[2]], nodeData[[3]], coeffSolutions, coeff, nodeIndex, n],
  (* Ordinary node: use A(t_j).y_{previous weight}. *)
  nodeData[[2]] . coeffSolutions[[coeff - 1, nodeIndex]]
];

(* Factor the small scalar collocation matrix once and reuse it for all weights. *)
CHESSBuildLinearSolver[matrix_, Automatic] :=
  (* Let LinearSolve pick the best method for the small scalar collocation block. *)
  Quiet[LinearSolve[matrix]];
CHESSBuildLinearSolver[matrix_, method_] :=
  (* Expose Method for experiments without changing the propagation code. *)
  Quiet[LinearSolve[matrix, Method -> method]];

(* Reject failed LinearSolve output before it contaminates later weight layers. *)
CHESSValidSolutionMatrixQ[solution_, dims : {_Integer?Positive, _Integer?Positive}] :=
  MatrixQ[solution] && Dimensions[solution] == dims && FreeQ[solution, _LinearSolve];

(* Solve the scalar collocation system against n right-hand sides at once.
   The full vector is interpreted as a node-by-component matrix of size
   (nodes) x (master integrals). *)
CHESSScalarCollocationSolve[scalarSolverOrMatrix_, rhs_, n_Integer?Positive, method_, prec_] := Module[
  {rhsMatrix, scalarSolver, solutionMatrix},

  (* New code passes rhs already as a node-by-master matrix.  The Partition path
     keeps compatibility with older callers that passed a flattened vector. *)
  rhsMatrix = SetPrecision[
    If[MatrixQ[rhs], rhs, Partition[rhs, n]],
    prec
  ];

  (* Accept either a pre-factored LinearSolveFunction or a raw matrix. *)
  scalarSolver = If[
    Head[scalarSolverOrMatrix] === LinearSolveFunction,
    scalarSolverOrMatrix,
    CHESSBuildLinearSolver[scalarSolverOrMatrix, method]
  ];
  If[Head[scalarSolver] =!= LinearSolveFunction, Return[$Failed]];

  (* Solving all master components at once reuses the same scalar factorization. *)
  solutionMatrix = Quiet[scalarSolver[rhsMatrix]];
  If[
    !CHESSValidSolutionMatrixQ[solutionMatrix, Dimensions[rhsMatrix]],
    Return[$Failed]
  ];
  SetPrecision[solutionMatrix, prec]
];

Options[SpectralPropagate] = {
  (* Chebyshev-Lobatto order; the actual collocation grid has Nodes+1 points. *)
  "Nodes" -> 48,

  (* Precision for the collocation solve and for numerical A(t) evaluation. *)
  "Precision" -> 160,
  "WorkingPrecisionA" -> 220,

  (* LinearSolve method for the scalar collocation block. *)
  "LinearSolverMethod" -> Automatic,

  (* Accepts {}, None, "Left", "Right", All, or explicit endpoint coordinates. *)
  "RegularizedEndpoints" -> {},

  (* Residue and finite-part extraction can use a separate working precision. *)
  "EndpointRegularizationPrecision" -> Automatic,

  (* Parallel stages are endpoint-letter extraction and node-matrix construction. *)
  "ParallelEvaluation" -> Automatic,
  "ParallelKernels" -> 1,

  (* Kept for older callers; the only implementation is the sequential-weight solver. *)
  "SolveStrategy" -> "SequentialEpsilon"
};

(* Main propagation routine.
   Epsilon coefficients are solved one layer at a time. The Chebyshev matrix is
   scalar in the master-integral index, so this avoids the huge Kronecker system
   and solves a small node matrix with many right-hand sides. *)
SpectralPropagateSequential[nAfun_, y0_, {x0_, x1_}, opts : OptionsPattern[SpectralPropagate]] := Module[
  {
    m, prec, precA, endpointPrec, collocation, nodes, derivMat, n, k,
    scalarBaseBlock, scalarSolver, solverMethod, coeffSolutions, rhs, solutionMat,
    values, finalValue, coeff, j, offset, regularizedEndpoints, activeRegularizedEndpoints, endpointSubBlocks,
    endpointInfo, parallelQ, nodeParallelQ, nodeSubBlocks, source,
    leftEndpointRegularized, zValues, batchEvaluatorQ, batchThreadCount,
    parallelSpec, kernelSpec
  },
  (* Read options once so the inner loops do not repeatedly call OptionValue. *)
  m = OptionValue["Nodes"];
  prec = OptionValue["Precision"];
  precA = OptionValue["WorkingPrecisionA"];
  endpointPrec = Replace[OptionValue["EndpointRegularizationPrecision"], Automatic -> precA];
  solverMethod = OptionValue["LinearSolverMethod"];

  (* Build collocation nodes and the scalar derivative matrix once. *)
  collocation = ChebyshevLobattoData[m, {x0, x1}, prec];
  nodes = N[collocation[[1]], prec];
  derivMat = N[collocation[[2]], prec];

  (* y0 is n x k: n master integrals and k epsilon coefficients. *)
  n = Dimensions[y0][[1]];
  k = Dimensions[y0][[2]];

  (* Left endpoint regularization changes the unknown from y to z=(y-y0)/(t-x0). *)
  regularizedEndpoints = NormalizeRegularizedEndpoints[OptionValue["RegularizedEndpoints"], {x0, x1}];
  leftEndpointRegularized = LeftEndpointRegularizedQ[regularizedEndpoints, x0];
  activeRegularizedEndpoints = If[
    TrueQ[leftEndpointRegularized],
    (* In the lifted formulation the left endpoint is an unknown z(0), so it
       must remain active rather than being replaced by a boundary row. *)
    Select[regularizedEndpoints, EndpointPointIndex[#, nodes] > 0 &],
    (* In the ordinary formulation the left endpoint row is the boundary
       condition, so only endpoints among Rest[nodes] are active. *)
    ActiveRegularizedEndpoints[regularizedEndpoints, nodes]
  ];

  parallelSpec = OptionValue["ParallelEvaluation"];
  kernelSpec = OptionValue["ParallelKernels"];
  batchEvaluatorQ = CHESSNodeEvaluatorQ[nAfun];
  batchThreadCount = If[
    batchEvaluatorQ,
    CHESSBatchThreadCount[parallelSpec, m, kernelSpec],
    1
  ];

  (* A native batch evaluator receives a thread count directly.  Wolfram
     subkernels are launched only if endpoint extraction itself needs them. *)
  parallelQ = If[
    batchEvaluatorQ && activeRegularizedEndpoints === {},
    False,
    CHESSResolveParallelQ[parallelSpec, m, kernelSpec]
  ];
  If[parallelQ, CHESSPrepareEndpointParallel[]];

  (* Convert each active spurious pole into sparse residue/finite matrices. *)
  endpointSubBlocks = Table[
    RegularizedEndpointOperatorData[point, {x0, x1}, endpointPrec, parallelQ],
    {point, activeRegularizedEndpoints}
  ];

  (* Public diagnostic metadata: {endpoint, directLetterCount, poleLetterCount}. *)
  endpointInfo = Table[
    With[{data = EndpointMatrixPoleData[point, {x0, x1}, endpointPrec, parallelQ]},
      {point, data[[3]], data[[4]]}
    ],
    {point, activeRegularizedEndpoints}
  ];

  nodeParallelQ = !batchEvaluatorQ && parallelQ &&
    (m + 2 - If[TrueQ[leftEndpointRegularized], 1, 2] >= 8);
  If[nodeParallelQ, CHESSPrepareNodeParallel[nAfun]];

  (* Build B(t_j) for ordinary nodes and endpoint operator data for endpoint nodes. *)
  nodeSubBlocks = BuildNodeOperatorData[
    nAfun,
    nodes,
    activeRegularizedEndpoints,
    endpointSubBlocks,
    prec,
    precA,
    nodeParallelQ,
    (* Ordinary propagation starts at node 2 because node 1 is the boundary row.
       Lifted left regularization starts at node 1 because z is unknown there. *)
    If[TrueQ[leftEndpointRegularized], 1, 2],
    batchThreadCount
  ];
  If[nodeSubBlocks === $Failed, Return[$Failed]];

  (* The scalar collocation matrix is the same for every master integral. *)
  scalarBaseBlock = N[
    If[
      TrueQ[leftEndpointRegularized],
      LiftedScalarCollocationMatrix[derivMat, nodes, x0],
      BaseScalarCollocationMatrix[derivMat]
    ],
    prec
  ];
  scalarSolver = CHESSBuildLinearSolver[scalarBaseBlock, solverMethod];
  If[Head[scalarSolver] =!= LinearSolveFunction, Return[$Failed]];

  (* coeffSolutions[[c,j]] is the vector of all master integrals at epsilon
     coefficient c and node j.  Previous coefficients drive later ones. *)
  coeffSolutions = ConstantArray[0, {k, m + 1, n}];
  If[
    TrueQ[leftEndpointRegularized],
    (* Left endpoint: solve for z; reconstruct y = y0 + (t-x0) z. *)
    Do[
      (* For c=1 the source is zero because there is no lower epsilon layer. *)
      rhs = ConstantArray[0, {m + 1, n}];
      If[
        coeff > 1,
        Do[
          (* Every row is a collocation equation for z, including the left
             endpoint finite equation. *)
          source = NodeSequentialSource[nodeSubBlocks[[j]], coeffSolutions, coeff, j, n];
          rhs[[j]] = source;,
          {j, 1, m + 1}
        ];
      ];
      rhs = SetPrecision[rhs, prec];
      solutionMat = CHESSScalarCollocationSolve[scalarSolver, rhs, n, solverMethod, prec];
      If[
        solutionMat === $Failed,
        Message[SpectralPropagate::solvefail, coeff, Head[solutionMat], Dimensions[solutionMat]];
        Return[$Failed];
      ];
      zValues = solutionMat;

      (* Convert solved z-values back to physical y-values immediately so the
         next epsilon layer can use the same coeffSolutions layout as usual. *)
      coeffSolutions[[coeff]] = Table[
        SetPrecision[y0[[All, coeff]] + (nodes[[j]] - x0) zValues[[j]], prec],
        {j, 1, m + 1}
      ];,
      {coeff, 1, k}
    ];,
    (* Ordinary case: row 1 imposes the boundary value at x0, other rows impose y'. *)
    Do[
      rhs = ConstantArray[0, {m + 1, n}];

      (* The first row of BaseScalarCollocationMatrix is the boundary row. *)
      rhs[[1]] = y0[[All, coeff]];
      If[
        coeff > 1,
        Do[
          (* Other rows are driven by B(t_j).y_{c-1}(t_j), or by the finite
             endpoint source if the node is a regularized right endpoint. *)
          source = NodeSequentialSource[nodeSubBlocks[[j - 1]], coeffSolutions, coeff, j, n];
          rhs[[j]] = source;,
          {j, 2, m + 1}
        ];
      ];
      rhs = SetPrecision[rhs, prec];
      solutionMat = CHESSScalarCollocationSolve[scalarSolver, rhs, n, solverMethod, prec];
      If[
        solutionMat === $Failed,
        Message[SpectralPropagate::solvefail, coeff, Head[solutionMat], Dimensions[solutionMat]];
        Return[$Failed];
      ];
      coeffSolutions[[coeff]] = solutionMat;,
      {coeff, 1, k}
    ];
  ];

  (* Repackage node values and final epsilon series in the legacy result shape. *)
  (* Each node value is flattened as {all eps coefficients for master 1, ...},
     matching older downstream scripts. *)
  values = Table[
    Flatten[coeffSolutions[[All, j]], 1],
    {j, 1, m + 1}
  ];

  (* Final state is n x k at the right endpoint. *)
  finalValue = Transpose[coeffSolutions[[All, -1]]];

  (* Result is a compact positional list:
     {nodes, nodeValues, finalState, endpointSubBlocks, endpointInfo, strategy}. *)
  {nodes, values, finalValue, Transpose[{activeRegularizedEndpoints, endpointSubBlocks}], endpointInfo, "SequentialEpsilon"}
];

(* Public entry point. *)
SpectralPropagate[nAfun_, y0_, interval_, opts : OptionsPattern[]] :=
  SpectralPropagateSequential[nAfun, y0, interval, opts];


(*
  Diagnostics:
  If the left endpoint is claimed to be a spurious pole, this should be tiny
  after boundary values have been loaded by the calling notebook:
     EndpointRegularityResidual[0, boundaryValues, {0, 1}, 80]
*)
