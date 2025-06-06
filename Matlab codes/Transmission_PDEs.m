% The provided code incorporates exact solutions for the fields A and B, representing their propagation in 
% a Distributed Feedback (DFB) laser. 
% To validate the accuracy of the code, these exact solutions are inserted into the coupled system of
% partial differential equations (PDEs), and the resulting transmissions of the fields are visualized. 
% Additionally, a comparison is made between the transmission obtained using the exact solutions and the 
% one computed through the code verification process. 
% The second figure depicts the error between these two transmissions.

%% THE CODE TAKES ABOUT 10 MINUTES FOR COMPILATION

close all
clear all
format long


L = 300e-6;                 % Lenght of the active medium
T = 1e-12;                  % Time of observation
c = 3e8;                    % Speed of light
g = 0;                      % Gain of the laser
dn=2e-3;                    % Amplitude of the modulation of the refractive index
n_0=3;                      % Average value of the refractive index
lambda=0.1e-6;              % Period of the modulation of the refractive index
beta_0 = pi/lambda;         % Propagation costant at Bragg's condition
q0= (dn/(2*n_0))*beta_0;    % Coupling parameter in case of uniform grating
u0 = c/n_0;                 % Speed of the wave in the medium

n_k=200;                    % Number of nodes in delta_beta_L space
db=20/n_k;                  % Discretization of the delta_beta_L space
delta_beta_L=[-10:db:10];
delta_beta=delta_beta_L/L;  % Detuning parameter


Nx = 2000;                  % Spatial nodes
Nt = 5000;                  % Temporal nodes
dt = T / Nt;                % Temporal step mesh
dx = L / Nx;                % Spatial step mesh
x = [0:dx:L]';              % Spatial mesh
xint = [dx:dx:L-dx]';  

% Verification of the stability
s=u0*(dt/dx)

% Initialization of the vectors for the transmission in the different cases
err_transm=zeros(1,n_k+1);      % Transmission's error 
Transmission=zeros(1,n_k+1);    % Transmission of the exact solutions
Transmission_nm=zeros(1,n_k+1); % Approximated transmission


% The for loop is employed to iterate through all potential values of beta, 
% facilitating the analysis of the propagation of various modes.
for k=1:n_k+1
    
    % Definition of parameters 
    delta_k = delta_beta(k) - 1j*g/2;
    gamma = sqrt(q0^2 - delta_k^2);
    Omega = delta_beta(k) * c/n_0;
    
    % A0 e BO initial conditions at t=0;
    A0 = 1;
    B0 = 1j*A0*(conj(q0)/gamma)*sinh(gamma*L)/(cosh(gamma*L)-1j*(delta_k/gamma)*sinh(gamma*L));

    % Exact solutions depending only on x
    Atilde = @(x) A0 *(cosh(gamma*x)+ 1j*(delta_k/gamma)*sinh(gamma*x))+ 1j*q0/gamma * sinh(gamma*x)*B0;
    Btilde = @(x) -1j*conj(q0)/gamma * sinh(gamma*x)*A0 + (cosh(gamma*x) - 1j*(delta_k/gamma)*sinh(gamma*x))*B0;
    
    % Total exact solutions: they depend both on space and time
    A=@(x,t) Atilde(x) * exp(- 1j*Omega*t);
    B=@(x,t) Btilde(x) * exp(- 1j*Omega*t);
    
    % The transmittion is a vector (in order to investigate all the possible frequencies of 
    % different modes) in which each element is obtained by the following formula
    Transmission(k)=abs(A(L,T)/A(0,T))^2;

    % gA0 e gBL are the injected signals
    gA0 = @(t) 0*t + Atilde(0)*exp(-1j*Omega*t); 
    gBL = @(t) 0*t; 

    % Discretization of time
    tempo=linspace(0,T,Nt);

    % The exact solutions are assigned to the approximated ones
    A_approx = A(xint,tempo);
    B_approx = B(xint,tempo);
    
    % Time initialization
    t = 0;
    
    % The subsequent pair of for loops implement the finite difference 
    % method to solve a pair of coupled advection equations.
    for n=1:Nt-1
    A_approx(1,n+1) = A_approx(1,n) - u0 * dt/dx * (A_approx(1,n)-gA0(t)) + 1j * u0 * dt * q0 * B_approx(1,n) + u0 * dt * g/2 * A_approx(1,n);
    B_approx(1,n+1) = B_approx(1,n) + u0 * dt/dx * (B_approx(2,n)-B_approx(1,n)) + 1j * u0 * dt * conj(q0) * A_approx(1,n) + u0 * dt * g/2 * B_approx(1,n);
        for m=2:Nx-2
            A_approx(m,n+1)=  A_approx(m,n) - u0 * dt/dx * (A_approx(m,n)-A_approx(m-1,n)) + 1j * u0 * dt * q0 * B_approx(m,n) + u0 * dt * g/2 * A_approx(m,n);
            B_approx(m,n+1)=  B_approx(m,n) + u0 * dt/dx * (B_approx(m+1,n)-B_approx(m,n)) + 1j * u0 * dt * conj(q0) * A_approx(m,n) + u0 * dt * g/2 * B_approx(m,n);
        end
    A_approx(Nx-1,n+1) = A_approx(Nx-1,n) - u0 * dt/dx * (A_approx(Nx-1,n)-A_approx(Nx-2,n)) + 1j * u0 * dt * q0 * B_approx(Nx-1,n) + u0 * dt * g/2 * B_approx(Nx-1,n);
    B_approx(Nx-1,n+1) = B_approx(Nx-1,n) + u0 * dt/dx * (gBL(t)-B_approx(Nx-1,n)) + 1j * u0 * dt * conj(q0) * A_approx(Nx-1,n) + u0 * dt * g/2 * B_approx(Nx-1,n);

    % Time updating
    t = t+dt;
    end

    % Transmission obtained from the finite difference method
    Transmission_nm(k)=abs(A_approx(Nx-1,Nt)/A_approx(1,Nt))^2;

end

% Absolute error between the two transmissions
err_transm=abs(Transmission_nm-Transmission);

figure(1)
% Plot of the Exact Transmission over delta_beta_L
plot(delta_beta_L,Transmission,'b',LineWidth=2) 
hold on
% Plot of the Approximated Transmission over delta_beta_L
plot(delta_beta_L,Transmission_nm,'rx', LineWidth=2)
xlabel('Δβ*L')
ylabel('Transmission')
legend('Transimission','Transmission nm')
title(sprintf('gain = %d , L = 300 µm', g),"FontSize", 16)
grid on

figure(2)
% Plot of the error
plot(delta_beta_L,err_transm,'r',LineWidth=2) 
xlabel('Δβ*L')
ylabel('Error')
title('Error', 'Fontsize', 16)
grid on

