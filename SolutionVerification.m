%% Solution Verification
% Author: Sean McInerney; Last Update: 08/08/2018
%
% Used to generate Figure 3.2

%% Set-up

clear
clc
close all

% Parameters
D1 = 0.09;
D2 = 0.009;
L1 = 1.6;
L2 = 4;

% Spatial discretisation
N = 201; % Select N such that K is an integer
K = (N-1)*(L1/L2)+1;

% Find appropriate point to stop temporal sampling
thresh = 0.99; % Once this temperature is reached, temporal sampling stops
fun = @(t) Model_Solve(L2,t,D1,D2,L1,L2) - thresh;
t_c = fzero(fun,[0,10^10*L2^2/D2]); % Initial guess is scaled
% Designate plot times
t_plot = [0.01,0.05,0.1,0.2,0.5,1]*t_c; % Plot times, dependent on t_c
nt = length(t_plot); % Number of plots

% Set up background colours
colour = [0.8275,0.5529,0.7176; 0.8745,0.8039,0.8784];

% Plot details
figure
hold on
rectangle('Position',[0,0,L1,1],'FaceColor',colour(1,:),'EdgeColor','none')
rectangle('Position',[L1,0,L2-L1,1],'FaceColor',colour(2,:),'EdgeColor','none')
xlabel('$x$ (mm)','Interpreter','Latex')
ylabel('Temperature','Interpreter','Latex')
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','Latex')
axis([0,L2,0,1])
annotation('arrow',[0.2,0.7],[0.25,0.9],'LineWidth',1.5)
box on

% Determine which nodes are plotted for numerical solution
nx = 11; % Number of bullet points in numerical visualisation
index = round(linspace(1,N,nx),0); % Assists in numerical visualisation

%% Semi-analytical Approach

x = linspace(0,L2,N); % Spatial points where solution is evaluated

soln = Model_Solve(x,t_plot,D1,D2,L1,L2); % Semi-analytical solution
for i = 1:nt
    plot(x,soln,'k','LineWidth',1.5);
end

%% Numerical Approach

% Discretisation
Delta_x = x(2) - x(1);
Delta_t = 10^(-6)*t_c;
t_current = Delta_t; % Keeps track of the time
% To ensure that the numerical solution is evaluated at the times
% designated in t_plot, Delta_t_current is introduced. It takes on the time
% until the next plot time, if t_current+Delta_t were to exceed the next
% plot time
Delta_t_current = Delta_t;
plot_count = 0; % How many of the plots have been plotted 

% Initial time (Refer to Equation (3.47)-(3.48))
T = zeros(N,1); % T is the numerical approximation to the temperature
T(1) = 1; % Dirichlet BC (Refer to Equation (3.49))

% Construct A (Refer to Equations (3.57)-(3.59))
A_L = [D1*ones(1,K-1),D2*ones(1,N-K-1),2*D2];
A_D = [0,-2*D1*ones(1,K-2),-(D1+D2),-2*D2*ones(1,N-K)];
A_U = [0,D1*ones(1,K-2),D2*ones(1,N-K)];
A = sparse(diag(A_L,-1)+diag(A_D)+diag(A_U,1));

while plot_count < nt
    
    % Iterate through time (Refer to Equation (3.56))
    T = T + (Delta_t_current/Delta_x^2)*A*T;
    
    % Check if deltat needs to be changed to obtain the plot time
    Delta_t_current = t_plot(plot_count+1) - t_current; % Time until next plot
    if Delta_t_current > Delta_t
        Delta_t_current = Delta_t; % Use standard time step if it would occur before
    end
    
    % Plot if necessary
    if t_current == t_plot(plot_count+1)
        plot_count = plot_count+1;
        plot(x(index),T(index),'k.','MarkerSize',18)
    end
    
    % Increment time
    t_current = t_current + Delta_t_current;
end