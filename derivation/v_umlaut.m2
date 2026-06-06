-- Symbolic derivation of the closed-form fundamental-matrix entries (F_coords)
-- used by the V-Umlaut 5-point solver; see compute_f_coords in Solver.cpp.
-- Run with Macaulay2: https://macaulay2.com/
--
-- Paper: Linear Fundamental Matrix Estimation from 7 or 5 Points, CVPR 2026.

FF = frac(QQ[x_1..x_3,y_1..y_3,
        s_(1,1)..s_(2,2)]);
R = FF[f_(2,1),f_(3,1),f_(1,2),
       f_(1,3),f_(2,3)];
F = matrix{
       {0,      f_(1,2), f_(1,3)},
       {f_(2,1),      0, f_(2,3)},
       {f_(3,1),      1, 0      }
       };
I = id_(R^3);
xs = apply(3, i -> I_{i});
xs = xs | {sum xs}
ys = xs;
comb = (pts, i, j, c) ->
       c*pts#i + (1-c) * pts#j;
xs = xs | {
       comb(xs,0,1,s_(1,1)),
       comb(xs,0,2,s_(1,2)),
       matrix{{x_1},{x_2},{x_3}}};
ys = ys | {
       comb(ys,0,1,s_(2,1)),
       comb(ys,0,2,s_(2,2)),
       matrix{{y_1},{y_2},{y_3}}};
Ilin = ideal apply(xs, ys, (x,y) ->
       transpose y * F * x);
netList Ilin_*
eq4 = (Ilin_*)#4;
eq5 = (Ilin_*)#5;
factor((det F)%(ideal(eq4, eq5)))
c1s = s_(1,1)*s_(2,2)*
       (1-s_(1,2))*(1-s_(2,1))
c2s = s_(1,2)*s_(2,1)*
       (1-s_(1,1))*(1-s_(2,2))
I = Ilin + ideal(c1s*f_(2,3)+c2s)
netList apply(gens R, g -> g => g%I)
