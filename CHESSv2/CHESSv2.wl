(* ::Package:: *)

(*
  CHESSv2 -- unified pure-Mathematica spectral propagation package.

  This file is intentionally only a loader.  The numerical implementations are
  separated by responsibility so that a reader can inspect the canonical,
  non-canonical, fake-delta, and dispatch logic independently:

    Core/Canonical.wl     original canonical CHESS implementation;
    Core/NonCanonical.wl  original polynomial-epsilon extension;
    Core/FakeDelta.wl     B0-only auxiliary-delta propagation and order choice;
    Core/Dispatch.wl      the single public SpectralPropagate entry point.

  The canonical and non-canonical cores are loaded first.  Dispatch.wl is
  loaded last because it deliberately replaces their two public entry points
  with one small routing layer while continuing to call the original solver
  functions underneath.
*)

CHESSv2::load = "Failed to load the CHESSv2 module `1` from `2`.";

Module[{packageDirectory, moduleFiles, moduleFile, loaded},
  packageDirectory = If[
    StringQ[$InputFileName] && $InputFileName =!= "",
    DirectoryName[ExpandFileName[$InputFileName]],
    Directory[]
  ];

  moduleFiles = {
    {"Core", "Canonical.wl"},
    {"Core", "NonCanonical.wl"},
    {"Core", "MatrixAdapters.wl"},
    {"Core", "FakeDelta.wl"},
    {"Core", "Dispatch.wl"}
  };

  Do[
    moduleFile = FileNameJoin[Prepend[moduleFiles[[i]], packageDirectory]];
    loaded = Quiet @ Check[Get[moduleFile], $Failed];
    If[loaded === $Failed,
      Message[CHESSv2::load, StringRiffle[moduleFiles[[i]], "/"], moduleFile];
      Abort[]
    ],
    {i, Length[moduleFiles]}
  ];
];

$CHESSv2Version = "1.0.0";
