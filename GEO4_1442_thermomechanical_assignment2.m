clearvars; close all; clc;
%% Parameters  
% code parameters
n = 500; % number of nodes, the bigger the more precise but more computing time
CFL = 0.99; % Courant–Friedrichs–Lewy condition making dt smaller
StopCond = 3000*365.25*24*60*60; % set at what time it should stop, 3000 years
plot_speed = 4; % set number of iteration before plot (higher = faster)

% physics parameters
T_bc = 1000; % boundary condition [deg K]
T_intr = 1300; % intrusion temperature [deg K]
L0_intr = 1e4; % top of inculsion depth [m]
L1_intr = 11e3; % bottom of inclusion depth [m]

Ly = 15000; % depth [m]
rho = 3000; % density [kg.m⁻³]
eta = 1e16; % viscosity [Pa.s]
k = 3; % thermal coeff [W/m/K]
Cp = 1000; % thermal capacity [J/kg/K]
dPdy = 31000; % pressure over height (y-ayis)
g = 10; % gravity acceleration

%% Functions
% physics
dy = Ly/(n-1); % set dy size relative to the node number, -1 because slot between [m]
kappa = k/rho/Cp; % compute kappa constant
dt = dy^2/2/kappa*CFL; % set dt according to dy [s]

% code intialisation
y = 0:dy:Ly; % create y vector with dy interval from 0 to L
T = zeros(n,1); % initiate temperature vector
L = zeros(n,n); % initiate stiffness matriy
ii=0; % initialisation of while itreation counter
time = 0; % initialize time counter for "for" loop

% boundary condition
L(1,1) = 1; % boundary condition for L top left
L(n,n) = 1; % boundary condition for L bottom right

for i=1:1:n % loop initiating temperature profile
    if y(1,i)>=L0_intr && y(1,i)<=L1_intr % create T vector with the temperature evolution
        T(i,1) = T_intr; % set intrusion temperature between L0 and L1 intrusion
    else
        T(i,1) = T_bc;
    end
end
R = T; % copy T array to response arrray = same size and same boundary conditions since boundary aren't touch by loops
T0 = T; % store initial tempearture vector

% plot
figure(1) % index figure
while round(time) < StopCond % loop over each row between 2 to end-1
    ii = ii+1; % update while iteration counter [-]
    time = time+dt; % update time counter [s]
    
    for i=2:1:n-1 % fill up solution array
        R(i,1) = T(i,1)/dt; % Right Hand side is T/dt
        L(i,i-1) = -kappa/(dy^2); % compute first diagonal L
        L(i,i) = 1/dt+2*kappa/(dy^2); % compute second diagonal L
        L(i,i+1) = -kappa*dy^-2; % compute third diagonal L
    end
    T = L\R; % compute linear solution of L*S=R for S with S being the new Temperature and R beeing old T/dt
   
    % Plots
    if mod(ii,plot_speed)==0 % allow to plot every few iteration, set by plot_speed
        plot(T0,-y/1000,'cyan',T,-y/1000,'r')
        xlim([950 1350])
        ylim([-15 0])
        title(['Elapsed time: ',num2str(round(time/(365.25*24*60*60))),' years'])
        xlabel(['Temperature [',char(176),'K]'])
        ylabel('Depth [km]')
        legend('Initial temperature profile','Updated temperature profile')
        grid on
        drawnow % update plot each time this function is executed
    end
end
disp(['Question 2C: The peak temperature after ',num2str(round(time/(365.25*24*60*60))),' years is ',num2str(round(max(T))),char(176),'K']) % answer question 2C
