% MechSolids Airfoil Project
% Lindsey Andrade, Ambika Goel, Caleb Kissel

function airfoil()

c = 1; %cord length meters

% Constants 
rho =  1.225; % density of air kg/m3
V = 200; % airspeed m/s
A = ; % cross section of the wing (for lift)

alpha_initial = .1; %initial angle of attack in radians

LiftForce = 2*pi*(1/2 * rho * V^2 * A); %lift Lift = 2pi(1/2 * rho * V^2 * A); 

C_lift = LiftForce/(.5 * rho * V^2 * c); %coefficient of lift

%     function res = integralFunction(phi)
%         T = C_lift * (alpha_initial + phi);
%         res = T/(G*J);
%     end
% 
%     function PHI = fsolveFunction(x)
% 


end




