close all;
clear vars;
clc;
% This function outputs the magnitude of maximum bending moment and
% its station for 80 angles of attack from 0 -> 20 degrees.

% Vars to be inputted in metric units

R = 0.075; % rocket's inner radius
V = 520.51; % max velocity
rho = 1.2250; % air density
D_n = 0.15; % diameter of nose cone in m
alpha_d = linspace(0,20,80); % angle of attack in deg
alpha_r = alpha_d * ((2*pi)/180); % AOA in rad
A = 0.25 * pi * (2*R)^2; % reference area
q = 0.5 * rho * V^2; % dynamic pressure
C_N_n = 2; % normal force coefficient of nosecone (depends on nosecone shape)
x_1 = 2.02755; % distance between the nosecone's center of pressure and the rocket's CG
x_2 = 0.79485; % distance between the rocket's CG and the fins' center of pressure

% fin vars
l_t = 0.09375; % tip chord
l_r = 0.275; % root chord
l_s = 0.19375; % fin semi-span
theta = atan((1/l_r) * ((l_r - l_t) + .5*(l_t - l_r))); %sweep angle (rad) of fin
l_m = l_s/cos(theta); % mid chord
n = 4; % number of fins
A_f = 2 * (l_s * l_t); % area of trapezoidal fin
f = 1; % interference coefficient (1 for 3 or 4 fins, 0.5 for 6 fins)
X_B = 2.9972; % distance from nosecone tip to beginning of root chord
X_R = 0.1905; % parellel distance from leading edge to beginning of tip chord

% Fin Calcs based on their geometry and quantity
C_N_f = (1 + ((f*R)/(l_s+R))) * ((4*n*(l_s/(2*R))^2)/(1+sqrt(1+((2*l_m)/(l_r+l_t))^2))); % normal force coefficient of fins
X_CPF = X_B + ((X_R/3) * ((l_r + 2*l_t)/(l_r + l_t))) + (1/6) * ((l_r + l_t) - ((l_r*l_t)/(l_r+l_t))); % distance to fins center of pressure from nosecone tip

M_max = zeros(size(alpha_d));
X_max = zeros(size(alpha_d));
N_nose_data = zeros(size(alpha_d));
N_fins_data = zeros(size(alpha_d));
for k = 1:80
    % calculate normal force on nose cone and fins (depends on alpha)
    N_nose = q * A * alpha_r(k) * C_N_n; % normal force on nosecone in newtons
    N_nose_data(k) = N_nose;
    N_fins = q * A * alpha_r(k) * C_N_f; % normal force on fins in newtons
    N_nose_data(k) = N_fins;
    % distributed loads acting above and below the rocket's CG
    
    w2 = (N_fins*(2*x_2 + x_1) - (N_nose*x_1)) / (x_2^2 + (x_1*x_2)); % distributed load below CG
    w1 = (N_nose + N_fins - (w2*x_2)) / x_1; % distributed load above CG

    % shear force, max bending, max bending location; Nakka eqs, assuming
    % that max bending is aft of CG
    
    V1 = N_nose - (w1*x_1); % shear force
    L = x_1 + x_2;
    M_max(k) = (V1^2/(2*w2)) - V1*(L-x_1) + w2*(0.5*L^2 + 0.5*x_1^2 - L*x_1); % max bending moment in newtons/m
    X_max(k) = (V1/w2) + x_1;  % max bending location in m
end
%coeffs = polyfit(alpha_d, M_max,1)

% fprintf("Max bending moment is %f Nm, located %f meters down from the nosecone \n", M_max, X_max);