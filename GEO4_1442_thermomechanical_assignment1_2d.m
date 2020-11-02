clearvars; close all; clc;
%% Parameters
% code parameters
nx = 50; % number of nodes on x-axis, the bigger the more precise but more computing time
ny = 50; % number of nodes on y-axis

% physics parameters
v_bc = 0; % boundary condition (start and end) [m.s⁻¹]
xsize = 100; % length dyke on x-axis [m]
ysize = 100; % length dyke on y-axis [m]
eta = 1e16; % viscosity [Pa.s]
rho = 2800; % density [kg.m⁻³]
dPdy = 31000; % pressure over height (y-axis) [Pa.m⁻¹]
g = 10; % gravity acceleration [m.s⁻¹]

%% Functions
% physics
dx = xsize/(nx-1); % set dx size relative to the node number, -1 (because slots between)
dy = ysize/(ny-1);
RHS = (dPdy-rho*g)/eta; % Righ hand side of the equation, constants
conversion_year = 100*60*60*24*365.25; % conversion from m.s⁻¹ to cm.y⁻¹

% code initialisation
x = 0:dx:xsize; % create x vector with dx interval from 0 to L
y = 0:dy:ysize; % create y vector with dy interval
R = zeros(nx*ny,1); % initiate response vector, square because grid n*n
L = zeros(nx*ny,nx*ny); % initiate stiffness matrix
T = zeros(nx,ny);

% boundary condition
% R(1,1) = v_bc; % boundary condition left
% R(nx*ny,1) = v_bc; % boundary condition right
% L(1,1) = 1; % boundary condition for L top left
% L(nx*ny,nx*ny) = 1; % boundary condition for L bottom right

% compute solution
for i=1:1:ny
    for j=1:1:nx % loop to fill matrix
        k = (j-1)*ny+i;
%         R(k,1) = RHS; % compute RHS vector
        if(i==1 || i==nx || j==1 || j==ny) % boundary conditions, add to take into account both start and end of each direction
            R(k,1)=v_bc;
            L(k,k)=1;
        else % filling matrix
            R(k,1) = RHS;
            L(k,k-ny) = dx^-2;
            L(k,k-1) = dy^-2; 
            L(k,k) = -2*dx^-2-2*dy^-2; 
            L(k,k+1) = dy^-2; 
            L(k,k+ny) = dx^-2;
        end
    end
end
S = L\R; % compute linear solution of L*S=R for S
S = reshape(S,[nx,ny]); % reshape S into a matrix of size nx*xy

%% Plots
figure(1)
mesh(x,y,-S.*conversion_year) % 2D plots
title('2D magma flow')
xlabel('[m]')
xlabel('[m]')
zlabel('Velocity [cm.yr⁻¹]')
colorbar