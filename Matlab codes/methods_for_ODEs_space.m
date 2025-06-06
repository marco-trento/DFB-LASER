% The provided code addresses the solution of coupled mode equations with respect to spatial  
% coordinates, assuming the fields are independent of time. 
% The Euler and Runge-Kutta methods, are employed to solve these equations. 
% Additionally, the code includes the plotting of solutions obtained through the ode45 function.

close all
clear all
format long

global L c g dn n_0 lambda beta_0 q0 u0 a 

L = 300e-6;                     % Lenght of the active medium
c = 3e8;                        % Speed of the waves in the medium
g = 10000;                      % Gain of the laser
dn=2e-3;                        % Amplitude of the modulation of the refractive index
n_0=3;                          % Average value of the refractive index
lambda=0.1e-6;                  % Period of the modulation of the refractive index
beta_0 = pi/lambda;             % Propagation costant at Bragg's condition
q0= (dn/(2*n_0))*beta_0;        % Coupling parameter in case of uniform grating
u0 = c/n_0;                     % Speed of the wave in the medium

Nx = 5000;                      % Spatial nodes
dx = L/Nx;                      % Spatial step mesh
x = [0:dx:L]';                  % Spatial mesh

% Parameters needed for the exact solutions
delta_beta=-1.8e4;              % Detuning parameter
delta_k=delta_beta-1j*g/2;
gamma=sqrt(q0^2-delta_k^2);
Omega = delta_beta*c/n_0;
a=-1j*Omega;

% Definition of the spatial derivatives
F = @(x, A, B) (g/2 - n_0/c * a)*A + 1j*q0 *B;          
G = @(x, A, B) (-g/2 + n_0/c * a)*B - 1j*conj(q0) *A;

% Initial Conditions
A0=1;
B0=1j*[sinh(gamma*L/4)+cosh(gamma*L/4)];

% Initial Conditions for the exact solution
% A0 = whatever;
% B0 = 1j*A0*(conj(q0)/gamma)*sinh(gamma*L)/(cosh(gamma*L)-1j*(delta_k/gamma)*sinh(gamma*L));
% Definition of the exact solutions
% Atilde = @(x) A0 *(cosh(gamma*x)+ 1j*(delta_k/gamma)*sinh(gamma*x))+ 1j*(q0/gamma) * sinh(gamma*x)*B0;
% Btilde = @(x) -1j*conj(q0)/gamma * sinh(gamma*x)*A0 + (cosh(gamma*x) - 1j*(delta_k/gamma)*sinh(gamma*x))*B0;
% Aex = Atilde(x);
% Bex = Btilde(x);

%% EULER'S METHOD

% The field are vectors in which the first position is occupied by the initial conditions
A_euler= zeros(Nx+1,1);
A_euler(1)=A0;
B_euler = zeros(Nx+1,1);
B_euler(1)=B0;

% The for loop perform the Euler's method that is a one step method
for k=1:Nx
    A_euler(k+1) = A_euler(k) + dx*F(x(k),A_euler(k), B_euler(k));
    B_euler(k+1) = B_euler(k) + dx*G(x(k),A_euler(k), B_euler(k));
end


%% 3th ORDER RUNGA KUTTA METHOD

% The field are vectors in which the first position is occupied by the initial conditions
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


%% 4th ORDER RUNGA KUTTA METHOD 

% The field are vectors in which the first position is occupied by the initial conditions
A_rk4= zeros(Nx+1,1);
A_rk4(1)=A0;
B_rk4 = zeros(Nx+1,1);
B_rk4(1)=B0;

for w=1:Nx
    F1_4 = F(x(w), A_rk4(w), B_rk4(w));
    F2_4 = F(x(w)+dx/2, A_rk4(w) + dx*F1_4/2, B_rk4(w)+ dx*F1_4/2);
    F3_4 = F(x(w)+dx/2, A_rk4(w)+dx*F2_4/2, B_rk4(w)+dx*F2_4/2);
    F4_4 = F(x(w+1), A_rk4(w)+ dx*F3_4, B_rk4(w) +dx*F3_4);
    G1_4 = G(x(w), A_rk4(w), B_rk4(w));
    G2_4 = G(x(w)+dx/2, A_rk4(w) + dx*G1_4/2, B_rk4(w)+ dx*G1_4/2);
    G3_4 = G(x(w)+dx/2, A_rk4(w)+dx*G2_4/2, B_rk4(w)+dx*G2_4/2);
    G4_4 = G(x(w+1), A_rk4(w)+ dx*G3_4, B_rk4(w) +dx*G3_4);
    A_rk4(w+1) = A_rk4(w) + dx*(F1_4+2*F2_4+2*F3_4 + F4_4)/6;
    B_rk4(w+1) = B_rk4(w) + dx*(G1_4+2*G2_4+2*G3_4 + G4_4)/6;
end


%% ODE45

Y0 = [A0; B0]; % Initial condition vector
xspan = [0 L]; % Integration range from x=0 to x=10

% Call the ODE solver
[x45, Y] = ode45(@odefun, xspan, Y0);

% Extract solutions for A and B
A45 = Y(:, 1);
B45 = Y(:, 2);

% Plot of the real part of the fields A and B calculated with the different methods
figure;
subplot(2, 2, 1);
plot(x, real(A_euler),'b-',x,real(B_euler), 'r-',LineWidth=2)
xlabel('x')
ylabel('Amplitude')
grid on
legend('real(A euler)','real(B euler)')
title('Euler s method','FontSize',16)

subplot(2, 2, 2);
plot(x,real(A_rk3), 'b',x, real(B_rk3), 'r', LineWidth=2)
xlabel('x')
ylabel('Amplitude')
grid on
legend('real(A rk3)','real(B rk3)')
title('RK3 method','FontSize',16)

subplot(2, 2, 3);
plot(x,real(A_rk4), 'b',x, real(B_rk4), 'r', LineWidth=2)
xlabel('x')
ylabel('Amplitude')
grid on
legend('real(A rk4)','real(B rk4)')
title('RK4 method','FontSize',16)

subplot(2, 2, 4);
plot(x45, real(A45), 'b-',x45, real(B45), 'r-',LineWidth=2)
xlabel('x')
ylabel('Amplitude')
grid on
legend('real(A45)','real(B45)')
title('ODE45','FontSize',16)

% The errors can be calculated with respect to the exact solution when the right initial conditions
% are defined. See line 41 
% err_eul = abs(A_euler - Aex);
% err_rk3 = abs(A_rk3 - Aex);
% err_rk4 = abs(A_rk4 - Aex);
% figure
% plot(x, err_eul, 'yellow', x, err_rk3, 'bo', x, err_rk4, 'r', LineWidth=2)
% xlabel('x')
% ylabel('Errors')
% grid on
% legend('err eul','err rk3', 'err rk4')
% title('Errors','FontSize',16)

% Definition of the function required for the computation of the ode45 method
function dYdx = odefun(x45, Y)
    global L c g dn n_0 lambda beta_0 q0 u0 a
    % Unpack variables
    A = Y(1);
    B = Y(2);
    
    % Compute derivatives
    dA_dx = (g/2 - n_0/c * a) * A + 1j * q0 * B;
    dB_dx = (-g/2 + n_0/c * a) * B - 1j * conj(q0) * A;
    
    % Pack the derivatives into a column vector
    dYdx = [dA_dx; dB_dx];
end
