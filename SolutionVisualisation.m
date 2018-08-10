%% Solution Visualisation
% Author: Sean McInerney; Last Update: 14/06/2018
%
% Used to generate Figure 2.2

%% Set-up

clear
clc
close all

% Target parameters
D1_hat = 0.09;
D2_hat = 0.009;

% Spatial structure
L1 = 1.6;
L2 = 4;
N = 101;
x = linspace(0,L2,N);

% Probe location
p = L2;

% Set up background colours for Figure 2.2(a) and (c)
colour = [0.8275,0.5529,0.7176; 0.8745,0.8039,0.8784];

% Plot times
t_plot_a = [10,20,50,100,200,500,1000]; % For Figure 2.2(a)
t_plot_c = [10,200,500]; % For Figure 2.2(c)

% Parameters used for comparison in Figure 2.2(c)-(d)
D1 = [D1_hat,0.45];
D2 = [D2_hat,0.0077];

% Line colours in Figure 2.2(c)-(d)
% This would need to be adapted if more than two parameter pairs are compared
linecolours =    [0,0,0; 0.15,0.65,0.15];

% Find appropriate point to stop temporal sampling
threshold = 0.99; % Once this temperature is reached at probe, sampling stops
fun = @(t) Model_Solve(p,t,D1_hat,D2_hat,L1,L2) - threshold;
t_c = fzero(fun,[0,10^10*L2^2/D2_hat]); % Initial guess is scaled

% Data points
M = 100; % Number of temporal samples
delta_t = t_c/M;
t_data = delta_t:delta_t:t_c;

% Generate synthetic temperature data
T_hat = Model_Solve(p,t_data,D1_hat,D2_hat,L1,L2);

%% Figure 2.2(a)

% Plot details
figure
hold on
rectangle('Position',[0,0,L1,1],'FaceColor',colour(1,:),'EdgeColor','none')
rectangle('Position',[L1,0,L2-L1,1],'FaceColor',colour(2,:),'EdgeColor','none')
xlabel('$x$ (mm)','Interpreter','Latex')
ylabel('Temperature','Interpreter','Latex')
axis([0,L2,0,1])
set(gca,'FontSize',16)
box on
set(gca,'TickLabelInterpreter','Latex')
annotation('arrow',[0.2,0.7],[0.25,0.9],'LineWidth',1.5)

% Plot model solution
T = Model_Solve(x,t_plot_a,D1_hat,D2_hat,L1,L2);
plot(x,T,'k','LineWidth',1.5);

%% Figure 2.2(b)

% Plot details
figure
hold on
xlabel('$t$ (sec)','Interpreter','Latex')
ylabel('Temperature','Interpreter','Latex')
set(gca,'FontSize',16)
box on
set(gca,'TickLabelInterpreter','Latex')
axis([0,1.2*t_c,0,1])
plot([t_c,t_c],[0,1],'k--')

% Plot data
plot(t_data,T_hat,'k.')

%% Figure 2.2(c)

% Plot details
figure
hold on
rectangle('Position',[0,0,L1,1],'FaceColor',colour(1,:),'EdgeColor','none')
rectangle('Position',[L1,0,L2-L1,1],'FaceColor',colour(2,:),'EdgeColor','none')
axis([0,4,0,1])
xlabel('$x$ (mm)','Interpreter','Latex')
ylabel('Temperature','Interpreter','Latex')
set(gca,'FontSize',16)
box on
set(gca,'TickLabelInterpreter','Latex')
annotation('arrow',[0.2,0.7],[0.25,0.9],'LineWidth',1.5)

% Loop over the different parameter pairs
for ii = 1:length(D1)
    T_ii = Model_Solve(x,t_plot_c,D1(ii),D2(ii),L1,L2);
    % Loop over plot times
    for jj = 1:length(t_plot_c)
        % Plot solution using different parameter pairs
        plot(x,T_ii(:,jj),'Color',linecolours(ii,:),'LineWidth',1.5);
    end
end

%% Figure 2.2(d)

% Plot details
figure
hold on
xlabel('$t$ (sec)','Interpreter','Latex')
ylabel('Temperature','Interpreter','Latex')
set(gca,'FontSize',16)
box on
set(gca,'TickLabelInterpreter','Latex')
plot([t_c,t_c],[0,1],'k--')
axis([0,1.2*t_c,0,1])

% Loop over different parameter pairs
for ii = 1:length(D1)
    % Plot solution at probe location using different parameter pairs
    T_ii = Model_Solve(p,t_data,D1(ii),D2(ii),L1,L2);
    plot(t_data,T_ii,'Color',linecolours(ii,:),'LineWidth',1.5)
end
