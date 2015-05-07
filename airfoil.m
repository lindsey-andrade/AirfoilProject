% MechSolids Airfoil Project
% Lindsey Andrade, Ambika Goel, Caleb Kissel

function airfoil()
close all
% c = 1; %cord length meters

% Constants 
rho =  1.225; % density of air kg/m3
V = 245; % airspeed m/s
G = 24e9; %N/m^2 
E = 68.9e9; 

alpha_initial = 0.1; %initial angle of attack in radians

% % solid wing
% [X1, T1, Phi1, PhiTotal1, LiftForce1, disp1, FOS1] = calculation(20.4124, 0.0046, 0.0046, 0, 10, 1, 0.4679);
% 
% % 2cm shelled hollow wing
% [X2, T2, Phi2, PhiTotal2, LiftForce2, disp2, FOS2] = calculation(20.4124, 0.0028, 0.0028, 0, 10, 1, 0.4679); 
% 
% % 1cm shelled hollow wing
% [X3, T3, Phi3, PhiTotal3, LiftForce3, disp3, FOS3] = calculation(20.4123, 0.0016, 0.0016, 0, 10, 1, 0.4679);

% 6cm shelled 3m cord
[X4, T4, Phi4, PhiTotal4, LiftForce4, disp4, FOS4] = calculation(61.2347, 0.2275, 0.2275, 0.2234, 10, 3, 1.4598);

% 3m cord, solid wing
[X5, T5, Phi5, PhiTotal5, LiftForce5, disp5, FOS5] = calculation(61.2347, 0.3731, 0.3731,  0.3676, 10, 3, 1.4598);

% 5cm shelled 3m cord
[X6, T6, Phi6, PhiTotal6, LiftForce6, disp6, FOS6] = calculation(61.2347, 0.1987,  0.1987, 0.1950, 10, 3, 1.4598);

% 4cm shelled 3m cord
[X7, T7, Phi7, PhiTotal7, LiftForce7, disp7, FOS7] = calculation(61.2347, 0.1665,  0.1665, 0.1634, 10, 3, 1.4598);

% 3cm shelled 3m cord
[X8, T8, Phi8, PhiTotal8, LiftForce8, disp8, FOS8] = calculation(61.2347, 0.1304,  0.1304, 0.1279, 10, 3, 1.4598);

% 2cm shelled 3m cord
[X9, T9, Phi9, PhiTotal9, LiftForce9, disp9, FOS9] = calculation(61.2347, 0.0912,  0.0912, 0.0894, 10, 3, 1.4598);

% figure
% hold all
% plot(X1, disp1)
% plot(X2, disp2)
% plot(X3, disp3)
% title('Displacement from Bending of Wings')
% legend('Solid Wing', '2cm Shelled Wing', '1cm Shelled Wing', 'Location', 'Best')
% xlabel('Position over length of wing (m)')
% ylabel('Linear Displacement of the Wing (m)')
% 
% figure
% hold all
% plot(X1, PhiTotal1*180/pi)
% plot(X2, PhiTotal2*180/pi)
% plot(X3, PhiTotal3*180/pi)
% title('Twist of Wings, 1m Cord Length') 
% legend('Solid Wing','2cm Shelled Wing','1cm Shelled Wing', 'Location', 'Best')
% xlabel('Position over length of wing (m)')
% ylabel('Angular Displacement of the Wing (degrees)')
% 
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

figure
hold all
plot(X7, FOS7)
plot(X8, FOS8)
plot(X9, FOS9)
title('FOS for Different Wings')
legend('4cm Shelled Wing', '3cm Shelled Wing', '2cm Shelled Wing')
xlabel('Position over length of wing (m)')
ylabel('FOS')



    function [X, T, Phi, PhiTotal, LiftForce, disp, FOS] = calculation(A, J, Izz, Ixx, L, cord, centroid)

            LiftConst = (1/2 * rho * V^2 * A); %lift bullshit Lift = 2pi(1/2 * rho * V^2 * A); 
            f = 2*pi; % lift fudge factor
            C = LiftConst * f; %coefficient of lift
            
            d = abs(centroid - cord/4);

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
            
            % Find shear strain over wing due to twisting
            ShearStrain = zeros(1, length(X)); 
            
            for i = 1 : length(X)
                ShearStrain(i) = calculateStrain(cord, centroid, Phi(i), X(i));
            end
         
            % Find shear stress over wing due to twisting
            ShearStress = zeros(1, length(X)); 
            
            for i = 1 : length(X)
                ShearStress(i) = ShearStrain(i) * G;
            end
                ShearStress(1) = 0; 
            
            % Find normal stress due to bending
            y = 0.1795; 
            Sigma = zeros(1, length(X));
            for i = 1 : length(X)
                Mx = LiftForce(i) * X(i);
                Sigma(i) = Mx*y/Ixx;
                % Add up shear stresses to the things
                ShearStress(i) = ShearStress(i) + Sigma(i)/2; 
            end
            
            % Find normal stresses/VonMises/FOS
            VonMises = zeros(1, length(X));
            FOS = zeros(1, length(X));
            for i = 1 : length(X)
                Mat = [Sigma(i), ShearStress(i); ShearStress(i), 0];
                [Vectors, Values] = eig(Mat); 
                sigma1 = Values(1,1); 
                sigma2 = Values(2,2);
                VonMises(i) = sqrt(sigma1^2 + sigma2^2 - sigma1*sigma2); 
                
                yield = 276e6; 
                FOS(i) = yield/VonMises(i); 
            end
    end

    function ShearStrain = calculateStrain(cord, centroid, Phi, L)
           %twisting shearstrain
        x = cord-centroid; 
        
        if x > centroid
            c = x; 
        else
            c = centroid; 
        end
        
        ShearStrain = c*Phi/L; 
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





