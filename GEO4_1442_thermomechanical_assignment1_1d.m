clearvars; close all; clc;
%% Parameters
% code parameters
n = 40; % number of nodes, the bigger the more precise but more computing time

% physics parameters
v_bc = 0; % boundary condition (start and end) [m.s⁻¹]
xsize = 100; % width (x-axis) of the dyke [m]
eta = 1e16; % viscosity [Pa.s]
rho = 2800; % density [kg.m⁻³]
dPdy = 31000; % pressure over height (y-axis) [Pa.m⁻¹]
g = 10; % gravity acceleration [m.s⁻¹]

%% Functions
% physics
dx = xsize/(n-1); % set dx size relative to the node number, -1 (because slots between)
RHS = (dPdy-rho*g)/eta; % Righ hand side of the equation, constants
conversion_year = 100*60*60*24*365.25; % conversion from m.s⁻¹ to cm.y⁻¹

% code initialisation
x = 0:dx:xsize; % create x vector with dx interval from 0 to L
R = zeros(n,1); % initiate response vector
L = zeros(n,n); % initiate stiffness matrix

% boundary condition
R(1,1) = v_bc; % boundary condition left
R(n,1) = v_bc; % boundary condition right
L(1,1) = 1; % boundary condition for L top left
L(n,n) = 1; % boundary condition for L bottom right

% compute solution
for i=2:1:n-1 % loop to fill matrix
    R(i,1) = RHS; % create RHS vector of length 2 to n-1
    L(i,i-1) = dx^-2; % compute first diagonal L
    L(i,i) = -2*dx^-2; % compute second diagonal L
    L(i,i+1) = dx^-2; % compute third diagonal L
end
S = L\R; % compute linear solution of L*S=R for S

% compute analytical for veryfing the results
x_ana = 0:0.001:xsize; % analytical x, for comparison later
V_ana = (1/(2*eta))*(dPdy-rho*g)*(x_ana.^2 - xsize*x_ana); % anatycal solution for comparison, same answer but not able to modify as we could the other solution

%% Plots
figure(1)
hold on
plot(x,-S.*conversion_year,'r.-','MarkerSize',20,'LineWidth',1)
plot(x_ana,-V_ana.*conversion_year,'g','LineWidth',1)
title(['1D magma flow with ',num2str(n),' nodes'])
xlabel('Length [m]')
ylabel('Velocity [cm.yr⁻¹]')
legend('Linear solution result','Analytical solution result','Location','best')

