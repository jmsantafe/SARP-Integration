%% Analysis for Airbrakes Structural Spars: Buckling Failure
% Group: SARP, University of Washington
clear all;
clc;
%% Assumptions and Limitations
% Load case is max bending moment @ max q plus liftoff thrust (max axial).
% Most conservative case, as these two don't occur at the same instant
% All calculations are for a single spar assuming a constant stress over
% the face. Cases where the bending moment is not applied directly over a
% pair of spars. i.e. the rocket is slightly rotated, are not accounted
% for. 
%% Geometry Parameters and Loads
D = 6; % diameter of rocket
L = 15; % Length of spar
t = 0.05:0.05:1.0; % Thickness range  of spar
w = 0.1:0.05:1.5; % Width range of spar (radially)
[T, W] = meshgrid(t, w); % T: each row is a copy of t. W: each column is a copy of. Allows us to test all combos of t, w
y0 = 3 - W./2; % offset of spar from center
A = T .* W; % Finding the Area of the cross section
Ix = (1/12).*T.*W.^3 + A.*y0.^2; % Calcuating area moment of inertia of the beam relative to rocket longitudinal axis
MaxMoment = 2663; % Expected Aerodynamic Bending Load (lbf-in)
MaxStress = (MaxMoment .* y0) ./ Ix; % Stress caused by aerodynamic bending @ max Q
MaxAxial = 1325;  % Max axial load at liftoff (lbf)
MaxLoad   = MaxAxial + (MaxStress .* A); % Force on Spar (different for each possible size)
%% Buckling critical load
E = 10e6; % Elasticity Modulus (6061-T6) (psi)
yield = 40000; % Tensile Yield Strength (6061-T6) (psi)
K = 0.65; % effective length factor; k = .5 for fixed-fixed, but .65 is recommended
Le = K * L; % effective length
I_t = (1/12).*T.*W.^3; % Calculating moment of inertia of the beam about thickness axis
I_w = (1/12).*T.^3 .*W; % Moment of inertia for beam about radial width axis

for i=1:size(I_t,1) % rows
    for j=1:size(I_t, 2) % cols
        I(i,j) = min(I_t(i,j), I_w(i,j)); % smallest I for each possible shape
    end
end
r = sqrt(I ./ A); % radius of gyration
S = Le ./ r; % slenderness ratio
Scomp = sqrt(2*pi^2*E/yield); % Determination of whether or not to use Euler's Buckling Theorem or Johnson's Formula (Calculating the Slenderness)
% Element-wise selection: use Euler where S > Scomp, Johnson otherwise
% F is critical load
useEuler = S > Scomp; % binary matrix that shows where statement is true (1)
F = nan(size(S)); % empty matrix of all NaN; size matches meshgrid
F(useEuler)  = pi^2 * E .* I(useEuler) ./ Le.^2; % Euler's Buckling theorem; (useEuler) auto applies to elements where Euler is true
F(~useEuler) = A(~useEuler) .* (yield - ((yield .* Le)./(2*pi.*r(~useEuler))).^2 .* (1/E)); % Johnson's Formula (applies where Euler is false)
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
% Mark FOS at the design point: t = ", w = "
t_mark = 0.48;
w_mark = 0.63;
FOS_mark = interp2(T, W, FOS, t_mark, w_mark);   % FOS at the chosen geometry
plot(t_mark, w_mark, 'rp', 'MarkerSize', 14, ...
    'MarkerFaceColor', 'r', 'LineWidth', 1.5)
text(t_mark, w_mark, ...
    sprintf('  FOS = %.2f\n  (t = %.2f", w = %.2f")', FOS_mark, t_mark, w_mark), ...
    'Color', 'w', 'FontWeight', 'bold', ...
    'VerticalAlignment', 'bottom')
hold off
MOS_mark = FOS_mark - 2; % Margin of safety, using 2 as SARP baseline for analytical analysis only
%% Console summary at the design point
fprintf('--- Design point (t = %.2f", w = %.2f") ---\n', t_mark, w_mark);
fprintf('Max applied load   : %.2f lbf\n', interp2(T, W, MaxLoad, t_mark, w_mark));
fprintf('Buckling capacity  : %.2f lbf\n', interp2(T, W, F,        t_mark, w_mark));
fprintf('Factor of Safety   : %.2f\n',     FOS_mark);
fprintf('Margin of Safety   : %.2f\n', MOS_mark);
%% Export data for sheets
Tcol = T(:);
Wcol = W(:);
FOScol = FOS(:);
MaxLoadCol = MaxLoad(:);
BucklingLoadCol = F(:);
DataTable = table(Tcol, Wcol, FOScol, MaxLoadCol, BucklingLoadCol, ...
    'VariableNames', {'Thickness_in','Width_in','FOS','MaxLoad_lbf','BucklingLoad_lbf'});
writetable(DataTable, 'spar_FOS_results.csv');
%% References
% https://mechanicalc.com/calculators/column-buckling/validation gives the
% equation to calculate the critical slenderness ratio

