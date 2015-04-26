% MechSolids Airfoil Project
% Lindsey Andrade, Ambika Goel, Caleb Kissel

c = 1; %cord length meters

% Constants 
rho =  1.225; % density of air kg/m3
V = 200; % airspeed m/s
F = 100; % Force on wing N

alpha = ; % angle of attack

l = F*cos(alpha); %lift

C_lift = l/(.5 * rho * V^2 * c); %coefficient of lift

% Generalization equations 

