%% Optimal Probe Location
% Author: Sean McInerney; Last Update: 21/06/2018
%
% Used to generate Figures 2.5 and 2.7
%
% Uses TestParameterSpace.m to generate data
% Ensure that parameters in TestParameterSpace.m are appropriate and
% saveData = 'yes' is selected
%
% Is computational costly for high values of Npts

%% Set-up

clear
clc
close all

% Spatial Structure
L1 = 1.6;
L2 = 4;

% Different probe locations
PROBES = cell(40,1);
for ii = 1:20
    PROBES{ii} = 0.2*ii; % Single probe set-ups
    PROBES{ii+20} = [4,0.2*ii]; % Two probe set-ups
end

% Thresholds for indicator function.
epsilon = [1.5,1.0,0.5];

% Colour scheme of lines. Currently set up for three thresholds
map = [
    0,102/255,204/255 % Blue
    1,0,0; % Red
    0,204/255,102/255; % Green
    ];

% Set-up background colours
colour = [0.8275,0.5529,0.7176; 0.8745,0.8039,0.8784];

%% Generate Data using TestParameterSpace.m

% Loop over different probe set-ups 
for ii = 1:length(PROBES)
    p = PROBES{ii}(1);
    if length(PROBES{ii})==2 % If there is a second probe
        q = PROBES{ii}(2);
    else
        q = 'none';
    end
    TestParameterSpace(p,q) % Generate data
end

close all

%% Work with Data

nP = length(PROBES); % Number of differet probe set-ups
Neps = length(epsilon); % Number of different thresholds considered
% mathcalA will contain the proportion of area where I = 1 for each
% probe set-up and each choice of epsilon.
mathcalA = zeros(nP,Neps);
nSP = 0; % Number of single probe set-ups investigated
n2P = 0; % Number of two probe set-ups investigated
% To determine nSp and n2P:
for ii = 1:nP
    x_data = PROBES{ii}; % Probe locations in this set-up
    nx = length(x_data);
    if nx == 1
        nSP = nSP+1;
    elseif nx == 2
        n2P = n2P+1;
    end
end

singleProbeIndex = zeros(nSP,1); % Keeps track of single probes set-ups
twoProbeIndex = zeros(n2P,1); % Keeps track of two probe set-ups

singleProbes = zeros(nSP,1); % Locations of single probes
twoProbes = zeros(n2P,1); % Location of second probe in two probe set-ups

for ii = 1:nP % Loop over probe set-ups
    x_data = PROBES{ii};
    nx = length(x_data);
    p = x_data(1);
    if nx == 1 % If it is a single probe set up
        load(['ResultsForSingleProbeAt',num2str(10*p)])
        for jj = 1:Neps % Loop over different choices of epsilon
            % Find the parameter pairs where I_1(D1,D2|p,eps) \leq epsilon
            mathcalA(ii,jj) = length(find(Results<=epsilon(jj)));
        end
        index = find(singleProbeIndex==0,1); % Find next index to adjust
        singleProbeIndex(index) = ii;
        singleProbes(index) = p; % Store the probe location
    elseif nx == 2 % If it is a two probe set-up
        q = x_data(2);
        load(['ResultsForTwoProbesAt',num2str(10*p),'and',num2str(10*q)])
        for jj = 1:Neps % Loop over different choices of epsilon
            % Find the parameter pairs where I_1(D1,D2|p,eps) \leq epsilon
            mathcalA(ii,jj) = length(find(Results<=epsilon(jj)));
        end
        index = find(twoProbeIndex==0,1); % Find next index to adjust
        twoProbeIndex(index) = ii;
        twoProbes(index) = q; % Store the second probe location
    end
    mathcalA(ii,:) = mathcalA(ii,:)/numel(Results); % Recalculate area as a proportion
end

%% Produce plot

% If any single probe set-ups have been considered 
if nSP > 0
    
    % Plot details
    figure5 = figure;
    axes5 = axes('parent',figure5);
    hold on
    
    maxSingle = max(mathcalA(singleProbeIndex)); % Find upper bound for plot
    rectangle(axes5,'Position',[0,0,L1,1.2*maxSingle],'FaceColor',colour(1,:),'EdgeColor','none')
    rectangle(axes5,'Position',[L1,0,L2-L1,1.2*maxSingle],'FaceColor',colour(2,:),'EdgeColor','none')
    axis(axes5,[0,L2,0,1.2*maxSingle])
    box on
    set(axes5,'FontSize',12)
    xlabel(axes5,'$p$ (mm)','Interpreter','Latex')
    ylabel(axes5,'$\mathcal{A}_1$','Interpreter','Latex')
    set(axes5,'TickLabelInterpreter','Latex')
    
    % Pre-plot sorting of indices
    [singleProbes,index] = sort(singleProbes);
    singleProbeIndex = singleProbeIndex(index);
    
    % Plot
    for jj = 1:Neps
        plot(axes5,singleProbes,mathcalA(singleProbeIndex,jj),'.','Color',map(jj,:),'MarkerSize',18)
        plot(axes5,singleProbes,mathcalA(singleProbeIndex,jj),'-','Color',map(jj,:),'LineWidth',1.5)
    end
end

% If any two probe set-ups have been considered
if n2P > 0
    
    % Plot details
    figure7 = figure;
    axes7 = axes('parent',figure7);
    hold on
    
    maxTwo = max(mathcalA(twoProbeIndex)); % Find upper bound for plot
    rectangle(axes7,'Position',[0,0,L1,1.2*maxTwo],'FaceColor',colour(1,:),'EdgeColor','none')
    rectangle(axes7,'Position',[L1,0,L2-L1,1.2*maxTwo],'FaceColor',colour(2,:),'EdgeColor','none')
    axis(axes7,[0,L2,0,1.2*maxTwo])
    box on
    set(axes7,'FontSize',12)
    xlabel(axes7,'$q$ (mm)','Interpreter','Latex')
    ylabel(axes7,'$\mathcal{A}_2$','Interpreter','Latex')
    set(axes7,'TickLabelInterpreter','Latex')
    
    % Pre-plot sorting of indices
    [twoProbes,index] = sort(twoProbes);
    twoProbeIndex = twoProbeIndex(index);
    
    % Plot
    for jj = 1:Neps
        plot(axes7,twoProbes,mathcalA(twoProbeIndex,jj),'.','Color',map(jj,:),'MarkerSize',18)
        plot(axes7,twoProbes,mathcalA(twoProbeIndex,jj),'-','Color',map(jj,:),'LineWidth',1.5)
    end
end