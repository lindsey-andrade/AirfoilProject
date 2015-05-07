% MechSolids Airfoil Project
% Lindsey Andrade, Ambika Goel, Caleb Kissel

function airfoil()
close all
% c = 1; %cord length meters

% Constants 
rho =  1.225; % density of air kg/m3
V = 245; % airspeed m/s
%A = 20.4124; % cross section of the wing (for lift)
G = 24e9; %N/m^2 
%J = 0.0046;
%L = 10; % length of wing (m)
%Izz =  0.0046;
E = 68.9e9; 

alpha_initial = 0.1; %initial angle of attack in radians

% solid wing
[X1, T1, Phi1, PhiTotal1, LiftForce1, disp1] = calculation(20.4124, 0.0046, 0.0046, 10, 1);

% 2cm shelled hollow wing
[X2, T2, Phi2, PhiTotal2, LiftForce2, disp2] = calculation(20.4124, 0.0028, 0.0028, 10, 1); 

% 1cm shelled hollow wing
[X3, T3, Phi3, PhiTotal3, LiftForce3, disp3] = calculation(20.4123, 0.0016, 0.0016, 10, 1);

% 6cm shelled 3m cord
[X4, T4, Phi4, PhiTotal4, LiftForce4, disp4] = calculation(61.2347, 0.2275, 0.2275, 10, 3);

% 3m cord, solid wing
[X5, T5, Phi5, PhiTotal5, LiftForce5, disp5] = calculation(61.2347, 0.3731, 0.3731, 10, 3)

% 5cm shelled 3m cord
[X6, T6, Phi6, PhiTotal6, LiftForce6, disp6] = calculation(61.2347, 0.1987,  0.1987, 10, 3);

% 4cm shelled 3m cord
[X7, T7, Phi7, PhiTotal7, LiftForce7, disp7] = calculation(61.2347, 0.1665,  0.1665, 10, 3);

% 3cm shelled 3m cord
[X8, T8, Phi8, PhiTotal8, LiftForce8, disp8] = calculation(61.2347, 0.1304,  0.1304, 10, 3);

% 2cm shelled 3m cord
[X9, T9, Phi9, PhiTotal9, LiftForce9, disp9] = calculation(61.2347, 0.0912,  0.0912, 10, 3);

figure
hold all
plot(X1, disp1)
plot(X2, disp2)
plot(X3, disp3)
title('Displacement from Bending of Wings')
legend('Solid Wing', '2cm Shelled Wing', '1cm Shelled Wing', 'Location', 'Best')
xlabel('Position over length of wing (m)')
ylabel('Linear Displacement of the Wing (m)')

figure
hold all
plot(X1, PhiTotal1*180/pi)
plot(X2, PhiTotal2*180/pi)
plot(X3, PhiTotal3*180/pi)
title('Twist of Wings, 1m Cord Length') 
legend('Solid Wing','2cm Shelled Wing','1cm Shelled Wing', 'Location', 'Best')
xlabel('Position over length of wing (m)')
ylabel('Angular Displacement of the Wing (degrees)')

figure 
hold all
plot(X5, PhiTotal5*180/pi)
plot(X4, PhiTotal4*180/pi)
plot(X6, PhiTotal6*180/pi)
plot(X7, PhiTotal7*180/pi)
plot(X8, PhiTotal8*180/pi)
plot(X9, PhiTotal9*180/pi)
title('Twist of Wings, 3m Cord Length') 
legend( 'Solid Wing', '6cm Shelled Wing','5cm Shelled Wing', '4cm Shelled Wing', '3cm Shelled Wing', '2cm Shelled Wing', 'Location', 'Best')
xlabel('Position over length of wing (m)')
ylabel('Angular Displacement of the Wing (degrees)')


    function [X, T, Phi, PhiTotal, LiftForce, disp] = calculation(A, J, Izz, L, cord)

            LiftConst = (1/2 * rho * V^2 * A); %lift bullshit Lift = 2pi(1/2 * rho * V^2 * A); 
            f = 2*pi; % lift fudge factor
            C = LiftConst * f; %coefficient of lift
            
            d = abs(0.4679 - cord/4);

            X = linspace(0, L, 50);
            T = zeros(1, length(X));
            Phi = zeros(1, length(X));

            for i = 1 : length(X)
            T(i) = d * C * alpha_initial * exp((d*C*X(i))/(J*G)); 
            Phi(i) = alpha_initial * (exp(d*C*X(i)/(J*G)) - 1); 
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

    end
% plot(X, T)
% title('Torque acting on wing')
% 
% figure
% plot(X, Phi)
% title('Ammount of deflection long wing')
% 
% figure
% plot(X, PhiTotal)
% title('Total Angluar Displacement')
% 
% figure
% plot(X, LiftForce)
% title('Lift Force Over Wing')
% 
% figure
% plot(X, disp)
% title('Bending Displacement')
end





