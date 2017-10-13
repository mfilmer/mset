%% Set up SET parameters
e=1.60217e-19;

SET.Cs = 30e-18;
SET.Cd = 30e-18;
SET.Cg = 0.1e-18;
SET.Gs = 1e-6;
SET.Gd = 1e-6;
SET.T = 0;
SET.DeltaL = 3.4e-4*e;
SET.DeltaI = 3.4e-4*e;

Bias.Vs = 0;
Bias.Vd = 0;
Bias.Vg = 0;

%% Simulate
[G, vds, vgs] = basicset(SET, Bias);

%% Plot
pcolor(vgs, vds, abs(G));
shading flat;
colormap gray;
colorbar; xlabel('V_{gs} [V]');
ylabel('V_{ds} [V]');
title('SISIS at 0 K');