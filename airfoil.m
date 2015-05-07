% MechSolids Airfoil Project
% Lindsey Andrade, Ambika Goel, Caleb Kissel

function airfoil()
close all
c = 1; %cord length meters

% Constants 
rho =  1.225; % density of air kg/m3
V = 300; % airspeed m/s
A = 20.4124; % cross section of the wing (for lift)
G = 24e9; %N/m^2 
J = 0.0046;
L = 10; % length of wing (m)
Izz =  0.0046;
E = 68.9e9; 

alpha_initial = 0.2; %initial angle of attack in radians

LiftConst = 2*pi*(1/2 * rho * V^2 * A); %lift bullshit Lift = 2pi(1/2 * rho * V^2 * A); 

C = LiftConst/(.5 * rho * V^2 * c); %coefficient of lift

X = linspace(0, L, 50);
T = zeros(1, length(X));
Phi = zeros(1, length(X));

for i = 1 : length(X)
T(i) = C * alpha_initial * exp((C*X(i))/(J*G)); 
Phi(i) = alpha_initial * (exp(C*X(i)/(J*G)) - 1); 
end

PhiTotal = cumsum(Phi); 
AngleOfAttack = alpha_initial*ones(1, length(PhiTotal)) + PhiTotal; 


LiftForce = C * AngleOfAttack;

% Find bending displacement 
disp = zeros(1, length(X)); 

for i = 1 : length(X)
    x = X(i); 
    P = LiftForce(i); 
    disp(i) = P*x^2 / (6*E*Izz) * (3*L - x); 
end

plot(X, T)
title('Torque acting on wing')

figure
plot(X, Phi)
title('Ammount of deflection long wing')

figure
plot(X, PhiTotal)
title('Total Angluar Displacement')

figure
plot(X, LiftForce)
title('Lift Force Over Wing')

figure
plot(X, disp)
title('Bending Displacement')
end





