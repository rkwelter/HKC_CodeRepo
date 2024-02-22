% Simulator and heat transport calculator for the rotating Lorenz system
% by Roland Welter

clear all, close all, clc
addpath('../ModelConstruction');
addpath('../FluidSolver');
date = '2024_01_19_';

%% generate Data 

% generate variables
M=21; maxShell = floor((sqrt(8*M+1)-1)/2);
shear = 'yes';
[uPlusModes,uMinusModes,thetaModes] = VariableConstructor(M,shear);

% select parameters

Pr = 10;
k1 = 1/sqrt(2);
V = pi/sqrt(2*k1);
%Ra = [1 (1:1:4)*10^3/4 (1:1:4)*10^4/4 (1:1:4)*10^5/4 (3:1:20)*10^5/2 ]; % Large, course Rayleigh grid
Ra = [1 , 5:5:2000]; % [1 (1:1:10)*50 (6:1:10)*10^2 (2:1:5)*10^3 ]; % small, fine Rayleigh grid
smallRay = 1;
%Ro = (0:1:10)*10^2; % Large, course rotation grid
Ro = 0; % 0:50:300; % Small, fine rotation grid

% model name definition
if(strcmp(shear,'no'))
    diffEq = strcat('HKC',num2str(M),'NoShear');
else
    diffEq = strcat('HKC',num2str(M));
end

numUp = size(uPlusModes,1);
numUm = size(uMinusModes,1);
numTheta = size(thetaModes,1);
n=numUp+numUm+numTheta;

% build grid

RaGrid = ones(length(Ra),length(Ro));
RoGrid = ones(length(Ra),length(Ro));

for i =1:length(Ra)
    for j =1:length(Ro)
        RaGrid(i,j) = Ra(i);
        RoGrid(i,j) = Ro(j);
    end
end

heatTransportGrid = zeros(length(Ra),length(Ro));
limitStdDevPercentGrid = ones(length(Ra),length(Ro));

% Integrate

for i=1:length(Ra)
    for j =1:length(Ro)

        name = strcat('../SavedData/Trajectories_Prandtl_10/',diffEq,'_Trajectory_Ra_',num2str(Ra(i)),'_Ro_',num2str(Ro(j)),'.mat');
            
        if(isfile(name))
            allData = load(name);
            heatTransportGrid(i,j) = allData.heatTransport(end);
            limitStdDevPercentGrid(i,j) = allData.limitStdDevPercent;
        else
            message = strcat('Error: Ra ',num2str(Ra(i)),' Ro ',num2str(Ro(j)),' not found.');
            error(message);
        end 
    end
end

%surf(RoGrid,RaGrid,heatTransportGrid,'FaceAlpha',0.75);
if(smallRay)
    name = strcat('../SavedData/HeatTransport/',diffEq,'_HeatTransport_FineRes.mat'); % SmallRay.mat');
else
    name = strcat('../SavedData/HeatTransport/',diffEq,'_HeatTransport.mat');
end
save(name,'M','RaGrid','Pr','k1','RoGrid','heatTransportGrid','limitStdDevPercentGrid');
