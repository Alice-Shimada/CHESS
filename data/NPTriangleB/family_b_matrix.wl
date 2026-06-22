ClearAll[
  NPTriangleBPeriodData,
  NPTriangleBMatrix,
  NPTriangleBSubsectorMatrix
];

NPTriangleBPeriodData[yIn_, wp_] := Module[
  {y, r8, u1, u2, u3, u4, k2, psi0, wronskian, jacobian, yInvariant},
  y = SetPrecision[yIn, wp];
  r8 = Sqrt[1 - 8 y];
  u1 = -1;
  u2 = -((r8 + 1)^2)/4;
  u3 = -((r8 - 1)^2)/4;
  u4 = 0;
  k2 = (1 - 4 y - 8 y^2 - r8)/(1 - 4 y - 8 y^2 + r8);
  psi0 = N[(2/Pi) y^2 EllipticK[k2]/Sqrt[(u3 - u1) (u4 - u2)], wp];
  wronskian = y^3/((1 - 8 y) (1 + y));
  jacobian = psi0^2/wronskian;
  yInvariant = psi0^3/((1 - 8 y) wronskian^2);
  {psi0, wronskian, jacobian, yInvariant}
];

Options[NPTriangleBMatrix] = {"Precision" -> 100};
NPTriangleBMatrix[yIn_, OptionsPattern[]] := Module[
  {
    wp, y, sq4, psi0, wronskian, jacobian, yInvariant,
    f11, f21, a31, a33, g2, g3, mat, denom
  },
  wp = OptionValue["Precision"];
  y = SetPrecision[yIn, wp];
  sq4 = Sqrt[1 - 4 y];
  {psi0, wronskian, jacobian, yInvariant} = NPTriangleBPeriodData[y, wp];

  f11 = ((28 y^2 + 2 y + 1) psi0^2)/(3 y^4);
  f21 = ((2 y - 1) (88 y^3 + 84 y^2 + 66 y - 11) psi0^4)/(3 y^8);
  a31 = (64 (y + 1)^2 (8 y - 1) psi0^3)/(27 y^6);
  a33 = (2 (y + 1) (8 y - 1) psi0^2)/(3 y^4);

  g2 = (2/(1 + y)) {
      0, 10, -10, -8, -8, (7 - 8 y)/sq4, -9, 10, 12, 8,
      -4, -7, -58, -30, -8
    };
  g3 = (16/(3 y (y + 1))) {
      0, -4, 4, 2, 2, -4 (y + 1)/sq4, 0, (4 y + 2)/(1 - y),
      0, 4, -2, 1, -8 (y + 2)/(y - 1), -6 (y + 1)/(y - 1),
      (5 y + 4)/(y - 1)
    };

  mat = ConstantArray[0, {18, 18}];
  mat[[1, 1]] = f11/jacobian;
  mat[[1, 2]] = 1/jacobian;
  mat[[2, 1]] = f21/jacobian;
  mat[[2, 2]] = f11/jacobian;
  mat[[2, 3]] = yInvariant/jacobian;
  Do[mat[[2, 3 + j]] = (yInvariant/jacobian) g2[[j]], {j, 1, 15}];
  mat[[3, 1]] = a31/jacobian;
  mat[[3, 3]] = a33/jacobian;
  Do[mat[[3, 3 + j]] = g3[[j]], {j, 1, 15}];

  mat[[4, 10]] = -2/y;
  mat[[4, 11]] = 3/(2 y);
  mat[[4, 12]] = 1/y;
  mat[[4, 16]] = 1/y;
  mat[[4, 18]] = -3/(2 y);

  denom = y (1 + y);
  mat[[5, 5]] = 2/denom;
  mat[[5, 6]] = -2/denom;
  mat[[5, 13]] = -2/denom;
  mat[[5, 14]] = 1/denom;
  mat[[5, 15]] = -1/denom;
  mat[[5, 16]] = -2/denom;
  mat[[5, 17]] = -2/denom;

  mat[[6, 5]] = 1/y;
  mat[[6, 6]] = -1/y;
  mat[[6, 15]] = -1/(2 y);
  mat[[6, 16]] = -3/y;
  mat[[6, 17]] = -1/y;

  mat[[7, 7]] = -(1 + 3 y)/(y (1 + y));
  mat[[7, 8]] = -2/(1 + y);
  mat[[7, 11]] = (2 + 3 y)/(y (1 + y));
  mat[[7, 12]] = 2/y;
  mat[[7, 15]] = -1/(y (1 + y));
  mat[[7, 16]] = 2/y;
  mat[[7, 18]] = -(5 + 6 y)/(2 y (1 + y));

  mat[[8, 7]] = 1/y;
  mat[[8, 11]] = -3/(2 y);
  mat[[8, 12]] = -1/y;
  mat[[8, 16]] = -1/y;
  mat[[8, 18]] = 3/(2 y);

  mat[[9, 9]] = (3 - 4 y)/(y (1 - 4 y));
  mat[[9, 10]] = -3/(y sq4);
  mat[[9, 14]] = 1/(y sq4);
  mat[[9, 18]] = 1/(y sq4);
  mat[[10, 9]] = 1/(y sq4);
  mat[[10, 10]] = -1/y;

  mat[[11, 11]] = -1/y;
  mat[[11, 12]] = -2/y;
  mat[[11, 18]] = 1/y;

  mat[[12, 11]] = 1/(y (1 - y));
  mat[[12, 12]] = 2/y;
  mat[[12, 16]] = 4/(y (1 - y));
  mat[[12, 17]] = 4/(y (1 - y) (1 + y));
  mat[[12, 18]] = -3/(2 y (1 - y));

  mat[[13, 13]] = 2/y;
  mat[[14, 14]] = 1/y;
  mat[[16, 16]] = -1/y;
  mat[[16, 17]] = -2/(y (1 + y));
  mat[[17, 16]] = 3/y;
  mat[[17, 17]] = 4/(y (1 + y));

  SparseArray[N[mat, wp]]
];

Options[NPTriangleBSubsectorMatrix] = {"Precision" -> 100};
NPTriangleBSubsectorMatrix[yIn_, OptionsPattern[]] := Module[{mat},
  mat = NPTriangleBMatrix[yIn, "Precision" -> OptionValue["Precision"]];
  mat[[4 ;; 18, 4 ;; 18]]
];
