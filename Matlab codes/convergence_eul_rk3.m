% The following code aims to verify the convergence of both Euler and Runge-Kutta method.
% In particular, we'll calculate the order of convergence using the definition:
% p = -log2[ err(dx/2) / err(dx)]
% where the errors are calculated using the maximum norm.
% Since the exact values are needed to calculate the error, we'll use the initial conditions that lead to the exact solutions.

close all
clear all
format long

% Definition of the parameters of the problem, needed in both methods
L = 300e-6;                 % Length of the active medium
c = 3e8;                    % Speed of light
g = 0;                      % Gain of the laser
dn=2e-3;                    % Amplitude of the modulation of the refractive index
n_0=3;                      % Average valure of the refractive index
lambda=0.1e-6;              % Period of the modulation of the refractive index

beta_0 = pi/lambda;         % Propagation constant at Bragg's condition
q0= (dn/(2*n_0))*beta_0;    % Coupling parameter in the case of uniform grating
u0 = c/n_0;                 % Speed of the wave in the medium

% Parameters needed for the exact solutions
delta_beta = -1.8e5; 
delta_k = delta_beta - 1j*g/2;
gamma = sqrt(q0^2 - delta_k^2);
Omega = delta_beta * c/n_0;
a = -1j*Omega;

 % Definition of the initial conditions
A0 = 1;
B0 = 1j*A0*(conj(q0)/gamma)*sinh(gamma*L)/(cosh(gamma*L)-1j*(delta_k/gamma)*sinh(gamma*L));
%The initial conditions only depend on the specifics of the problem, so they will remain the same for every method and every Nx

%% Calculation of the order of convergence of the Euler method

% Implementation of the Euler method, in the case of Nx=50000;
Nx = 50000;                 % Spatial nodes
dx = L/Nx;                  % Spatial step mesh
x = [0:dx:L]';              % Spatial mesh

  % Definition of the exact solutions
Atilde = @(x) A0 *(cosh(gamma*x)+ 1j*(delta_k/gamma)*sinh(gamma*x))+ 1j*(q0/gamma) * sinh(gamma*x)*B0;
Btilde = @(x) -1j*conj(q0)/gamma * sinh(gamma*x)*A0 + (cosh(gamma*x) - 1j*(delta_k/gamma)*sinh(gamma*x))*B0;
Aex = Atilde(x);
Bex = Btilde(x);

  % Initialization of the numerical solutions
A_euler = zeros(Nx+1,1);
A_euler(1)=A0;
B_euler = zeros(Nx+1,1);
B_euler(1)=B0;

  % Computation of the method
for k=1:Nx
    A_euler(k+1) = A_euler(k)+dx*((g/2 - n_0/c * a)*A_euler(k) + 1j*q0 *B_euler(k));
    B_euler(k+1) = B_euler(k)+dx*((-g/2 + n_0/c * a)*B_euler(k) - 1j*conj(q0) *A_euler(k));
end

err_eul_1 = abs(A_euler-Aex);                            % Absolute error
abs_err_eul_inf1=norm(err_eul_1,inf)/norm(Aex,inf);      % Maximum norm error, in the case of Nx=50000

% Implementation of the Euler method in the case of Nx=100000 so that the step mesh is halved
Nx=2*Nx;
dx = L/Nx;
x = [0:dx:L]';

Atilde = @(x) A0 *(cosh(gamma*x)+ 1j*(delta_k/gamma)*sinh(gamma*x))+ 1j*(q0/gamma) * sinh(gamma*x)*B0;
Btilde = @(x) -1j*conj(q0)/gamma * sinh(gamma*x)*A0 + (cosh(gamma*x) - 1j*(delta_k/gamma)*sinh(gamma*x))*B0;
Aex = Atilde(x);
Bex = Btilde(x);

A_euler = zeros(Nx+1,1);
A_euler(1)=A0;
B_euler = zeros(Nx+1,1);
B_euler(1)=B0;

for k=1:Nx
    A_euler(k+1) = A_euler(k)+dx*((g/2 - n_0/c * a)*A_euler(k) + 1j*q0 *B_euler(k));
    B_euler(k+1) = B_euler(k)+dx*((-g/2 + n_0/c * a)*B_euler(k) - 1j*conj(q0) *A_euler(k));
end

err_eul_2 = abs(A_euler-Aex);                         
abs_err_eul_inf2=norm(err_eul_2,inf)/norm(Aex,inf);   % Maximum norm error in the case of Nx=100000

% Order of convergence of the Euler method
p_euler = -log2(abs_err_eul_inf2/abs_err_eul_inf1)


%% Calculation of the order of convergence of the four stage Runge-Kutta method

% Implementation of the four stage Runge-Kutta method with Nx=50000
Nx = 50000;
dx = L/Nx;
x = [0:dx:L]';

  % Definition of the spatial derivatives
F = @(x, A, B) (g/2 - n_0/c * a)*A + 1j*q0 *B;
G = @(x, A, B) (-g/2 + n_0/c * a)*B - 1j*conj(q0) *A;

  % Definition of the exact solutions
Atilde = @(x) A0 *(cosh(gamma*x)+ 1j*(delta_k/gamma)*sinh(gamma*x))+ 1j*(q0/gamma) * sinh(gamma*x)*B0;
Btilde = @(x) -1j*conj(q0)/gamma * sinh(gamma*x)*A0 + (cosh(gamma*x) - 1j*(delta_k/gamma)*sinh(gamma*x))*B0;
Aex = Atilde(x);
Bex = Btilde(x);

  % Initialization of the numerical solutions
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


err_rk3_1 = abs(A_rk3-Aex);                           
abs_err_rk3_inf1=norm(err_rk3_1,inf)/norm(Aex,inf);   % Maximum norm error with Nx=50000

% Implementation of the four stage Runge-Kutta method with Nx=100000
Nx = 2*Nx;
dx = L/Nx;
x = [0:dx:L]';

F = @(x, A, B) (g/2 - n_0/c * a)*A + 1j*q0 *B;
G = @(x, A, B) (-g/2 + n_0/c * a)*B - 1j*conj(q0) *A;

Atilde = @(x) A0 *(cosh(gamma*x)+ 1j*(delta_k/gamma)*sinh(gamma*x))+ 1j*(q0/gamma) * sinh(gamma*x)*B0;
Btilde = @(x) -1j*conj(q0)/gamma * sinh(gamma*x)*A0 + (cosh(gamma*x) - 1j*(delta_k/gamma)*sinh(gamma*x))*B0;
Aex = Atilde(x);
Bex = Btilde(x);

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


err_rk3_2 = abs(A_rk3-Aex);                            
abs_err_rk3_inf2=norm(err_rk3_2,inf)/norm(Aex,inf);   % Maximum norm error with Nx=100000

  % Order of convergence of the four stage Runge-Kutta method 
p_rk3 = -log2(abs_err_rk3_inf2/abs_err_rk3_inf1)
