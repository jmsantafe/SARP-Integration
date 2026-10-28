%% Analysis for Airbrakes Structural Spars: Buckling Failure
% Group: SARP, University of Washington
clear all;
clc;
%% Geometry Parameters and Loads
D = 6; % diameter of rocket
L = 11.5; % Length of spar
t = 0.05:0.05:1.0; % Thickness of spar
w = 0.1:0.05:1.5; % Width of spar (radially)
[T, W] = meshgrid(t, w);
y0 = 3 - W./2; % offset of spar from center
A = T .* W; % Finding the Area of the cross section
Ix = (1/12).*T.*W.^3 + A.*y0.^2; % Calcuating moment of inertia of the beam relative to rocket axis
MaxMoment = 4276; % Expected Aerodynamic Bending Load (lbf-in)
MaxStress = (MaxMoment .* y0) ./ Ix;
MaxLoad   = MaxStress .* A; % On Spar
%% Buckling critical load
E = 10e6; % Elasticity Modulus (6061-T6) (psi)
yield = 40000; % Tensile Yield Strength (6061-T6) (psi)
K = 0.65; % effective length factor
Le = K * L; % effective length
I = (1/12).*T.*W.^3; % Calcuating moment of inertia of the beam
r = sqrt(I ./ A); % radius of gyration
S = Le ./ r; % slenderness ratio
Scomp = sqrt(2*pi^2*E/yield); % Determination of whether or not to use Euler's Buckling Theorem or Johnson's Formula (Calculating the Slenderness)
% Element-wise selection: use Euler where S > Scomp, Johnson otherwise
useEuler = S > Scomp;
F = nan(size(S));
F(useEuler)  = pi^2 * E .* I(useEuler) ./ Le.^2;                                       % Euler's Buckling theorem
F(~useEuler) = A(~useEuler) .* (yield - ((yield .* Le)./(2*pi.*r(~useEuler))).^2 .* (1/E)); % Johnson's Formula
FOS = F ./ MaxLoad; % Calculating the Factor of Safety
%% Plot
figure
contourf(T, W, FOS, 20)
colorbar
xlabel('Thickness t (in)')
ylabel('Width w (in)')
title('FOS Contours vs Spar Thickness and Width')
% Highlights FOS = 1 boundary
hold on
contour(T, W, FOS, [1 1], 'k', 'LineWidth', 2)
text(0.3,0.4,'FOS = 1','Color','k','FontWeight','bold')
% Mark FOS at the design point: t = 0.5", w = 0.6"
t_mark = 0.5;
w_mark = 0.6;
FOS_mark = interp2(T, W, FOS, t_mark, w_mark);   % FOS at the chosen geometry
plot(t_mark, w_mark, 'rp', 'MarkerSize', 14, ...
    'MarkerFaceColor', 'r', 'LineWidth', 1.5)
text(t_mark, w_mark, ...
    sprintf('  FOS = %.2f\n  (t = %.2f", w = %.2f")', FOS_mark, t_mark, w_mark), ...
    'Color', 'w', 'FontWeight', 'bold', ...
    'VerticalAlignment', 'bottom')
hold off
%% Console summary at the design point
fprintf('--- Design point (t = %.2f", w = %.2f") ---\n', t_mark, w_mark);
fprintf('Max applied load   : %.2f lbf\n', interp2(T, W, MaxLoad, t_mark, w_mark));
fprintf('Buckling capacity  : %.2f lbf\n', interp2(T, W, F,        t_mark, w_mark));
fprintf('Factor of Safety   : %.2f\n',     FOS_mark);
%% Export data for sheets
Tcol = T(:);
Wcol = W(:);
FOScol = FOS(:);
MaxLoadCol = MaxLoad(:);
BucklingLoadCol = F(:);
DataTable = table(Tcol, Wcol, FOScol, MaxLoadCol, BucklingLoadCol, ...
    'VariableNames', {'Thickness_in','Width_in','FOS','MaxLoad_lbf','BucklingLoad_lbf'});
writetable(DataTable, 'spar_FOS_results.csv');
