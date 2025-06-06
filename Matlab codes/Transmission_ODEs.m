% The provided code incorporates exact solutions for the fields A and B, representing their propagation 
% in a Distributed Feedback laser (DFB). 
% To validate the accuracy of the methods, these exact solutions are inserted into the coupled system of
% ordinary differential equations (ODEs), and the resulting transmissions of the fields are visualized.

%% THE CODE REQUIRES ABOUT 30 SECONDS TO COMPILE

close all
clear all
format long

% Definition of the parameters
L = 300e-6;                % Length of the active medium
Nx = 50000;                % Spatial nodes
dx = L/Nx;                 % Spatial step mesh
x = [0:dx:L]';             % Spatial mesh
c = 3e8;                   % Speed of light
g = 0;                     % Gain of the laser
dn=2e-3;                   % Amplitude of the modulation of the refractive index
n_0=3;                     % Average value of the refractive index
lambda=0.1e-6;             % Period of the modulation of the refractive index

beta_0 = pi/lambda;        % Propagation constant at Bragg's condition
q0= (dn/(2*n_0))*beta_0;   % Coupling parameter in the case of uniform grating
u0 = c/n_0;                % Speed of the wave in the medium

n_k=200;                   % Number of nodes in delta_beta_L space
db=20/n_k;                 % Discretization of the delta_beta_L space
delta_beta_L=[-10:db:10];  
delta_beta=delta_beta_L/L; % Detuning Parameter

% Initialization of the vectors for the transmission in the different cases
Transmission=zeros(1,n_k+1);       % Transmission of the exact solutions
Transmission_eul=zeros(1,n_k+1);   % Transmission of the solutions computed with Euler method
Transmission_rk=zeros(1,n_k+1);    % Transmission of the solutions computed with Runge-Kutta method

% This first for loop is employed to iterate through all potential values of beta, 
% facilitating the analysis of the propagation of various modes. 
for p=1:n_k+1
    
    % Definition of the parameters necessary for the exact solutions
    delta_k = delta_beta(p) - 1j*g/2;
    gamma = sqrt(q0^2 - delta_k^2);
    Omega = delta_beta(p) * c/n_0;
    a = -1j*Omega;

    % The spatial derivatives needed for the rk-method must be evaluated for
    % each value of a, and so they change upon changing the index p
    F = @(x, A, B) (g/2 - n_0/c * a)*A + 1j*q0 *B;
    G = @(x, A, B) (-g/2 + n_0/c * a)*B - 1j*conj(q0) *A;

    % Definition of the initial conditions, which depend on gamma, and so on the index p
    A0 = 1;
    B0 = 1j*A0*(conj(q0)/gamma)*sinh(gamma*L)/(cosh(gamma*L)-1j*(delta_k/gamma)*sinh(gamma*L));
    
    % Definition of the exact solutions as functions only depending on x 
    Atilde = @(x) A0 *(cosh(gamma*x)+ 1j*(delta_k/gamma)*sinh(gamma*x))+ 1j*(q0/gamma) * sinh(gamma*x)*B0;
    Btilde = @(x) -1j*conj(q0)/gamma * sinh(gamma*x)*A0 + (cosh(gamma*x) - 1j*(delta_k/gamma)*sinh(gamma*x))*B0;
    Aex = Atilde(x);
    Bex = Btilde(x);

    % The transmission is a vector (in order to investigate all the possible frequencies
    % of different modes) in which each element is obtained by the following formula
    Transmission(p)=abs(Aex(Nx)/Aex(1))^2;
    
    % Implementation of the Euler method to calculate the transmission
    A_euler = zeros(Nx+1,1);
    A_euler(1)=A0;
    B_euler = zeros(Nx+1,1);
    B_euler(1)=B0;
    
    for k=1:Nx
        A_euler(k+1) = A_euler(k)+dx*((g/2 - n_0/c * a)*A_euler(k) + 1j*q0 *B_euler(k));
        B_euler(k+1) = B_euler(k)+dx*((-g/2 + n_0/c * a)*B_euler(k) - 1j*conj(q0) *A_euler(k));
    end
    
    % Transmission obtained with Euler method
    Transmission_eul(p)=abs(A_euler(Nx)/A_euler(1))^2;
    
    % Implementation of the three stage Runge-Kutta method to calculate the transmission
    A_rk3= zeros(Nx+1,1);
    A_rk3(1)=A0;
    B_rk3 = zeros(Nx+1,1);
    B_rk3(1)=B0;
    
    for w=1:Nx
        F1_3 = F(x(w), A_rk3(w), B_rk3(w));
        F2_3 = F(x(w)+dx, A_rk3(w) + dx*F1_3, B_rk3(w)+ dx*F1_3);
        F3_3 = F(x(w)+dx/2, A_rk3(w)+dx*(F1_3+F2_3)/4, B_rk3(w)+dx*(F1_3+F2_3)/4);
        G1_3 = G(x(w), A_rk3(w), B_rk3(w));
        G2_3 = G(x(w)+dx, A_rk3(w) + dx*G1_3, B_rk3(w)+ dx*G1_3);
        G3_3 = G(x(w)+dx/2, A_rk3(w)+dx*(G1_3+G2_3)/4, B_rk3(w)+dx*(G1_3+G2_3)/4);
        A_rk3(w+1) = A_rk3(w) + dx*(F1_3+F2_3+4*F3_3)/6;
        B_rk3(w+1) = B_rk3(w) + dx*(G1_3+G2_3+4*G3_3)/6;
    end
    
    % Transimmision obtained with three stage Runge-Kutta method
    Transmission_rk(p)=abs(A_rk3(Nx)/A_rk3(1))^2;

end


figure(1)
% Plot of the Exact Transmission over delta_beta_L
plot(delta_beta_L,Transmission,'b',LineWidth=2) 
hold on
grid on
% Plot of the Transmission computed with Euler method
plot(delta_beta_L,Transmission_eul,'rx','MarkerSize',8,'MarkerFaceColor', 'red', LineWidth=2)
% Plot of the Transmission computed with three stage Runge-Kutta method
plot(delta_beta_L, Transmission_rk, 'o','MarkerSize',4,'MarkerFaceColor', 'yellow', LineWidth=2)
xlabel('Δβ*L')
ylabel('Transmission')
legend('Transimission','Transmission Euler', 'Transmission RK')
title(sprintf('gain = %d , L = 300 µm', g),"FontSize", 16)

% Definition of the error in the transmission in the two methods
errT_eul = abs(Transmission - Transmission_eul);    % absolute error of Euler method
errT_rk = abs(Transmission - Transmission_rk);      % absolute error of Runge-Kutta method
% Plot of both errors
figure(2)
plot(delta_beta_L, errT_eul, 'r', delta_beta_L, errT_rk, 'yellow', LineWidth=2)
grid on
xlabel('Δβ*L')
ylabel('Error')
legend('Error Euler', 'Error RK')
title('Errors','FontSize',16)


