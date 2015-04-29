% MechSolids Airfoil Project
% Lindsey Andrade, Ambika Goel, Caleb Kissel

function airfoil()
close all
c = 1; %cord length meters

% Constants 
rho =  1.225; % density of air kg/m3
V = 200; % airspeed m/s
A = 8; % cross section of the wing (for lift)
G = 24e9; %N/m^2 
J = 1300;
L = 100; % length of wing (m)

alpha_initial = .1; %initial angle of attack in radians

LiftForce = 2*pi*(1/2 * rho * V^2 * A); %lift Lift = 2pi(1/2 * rho * V^2 * A); 

C = LiftForce/(.5 * rho * V^2 * c); %coefficient of lift

X = linspace(0, L, 50);
T = zeros(1, length(X));
Phi = zeros(1, length(X));

for i = 1 : length(X)
T(i) = C * alpha_initial * exp((C*X(i))/(J*G)); 
Phi(i) = alpha_initial * (exp(C*X(i)/(J*G)) - 1); 
end

plot(X, T)

figure
plot(X, Phi)
end





