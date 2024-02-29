%% HeatTransport_Assembler by Roland Welter
% 
% Combines the heat transport data from a large number of trajectories into
% a single heatTransportGrid variable

clear all, close all, clc
addpath('../ModelConstruction');
addpath('../FluidSolver');
date = '2024_01_19_';

%% generate Data 

% generate variables
M=21; maxShell = floor((sqrt(8*M+1)-1)/2);
shear = 'yes';
[uPlusModes,uMinusModes,thetaModes] = VariableConstructor(M);

% select parameters

Pr = 10;
k1 = 1/sqrt(2); V = pi/sqrt(2*k1);
%Ra = [1 (1:1:4)*10^3/4 (1:1:4)*10^4/4 (1:1:4)*10^5/4 (3:1:20)*10^5/2 ]; % Large, course Rayleigh grid
%Ra = [1 (1:1:10)*50 (6:1:10)*10^2 (2:1:5)*10^3 ]; % small, fine Rayleigh grid
Ra = [1 , 5:5:2000]; % hyperfine Rayleigh grid
smallRay = 2; % determines name for file saving, 0 if large grid, 1 if small, 2 if hyperfine
%Ro = (0:1:10)*10^2; % Large, course rotation grid
%Ro = 0:50:300; % Small, fine rotation grid
Ro = 0; % non-rotating

% model name definition
diffEq = strcat('HKC',num2str(M));

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

        name = strcat('../SavedData/Trajectories/Prandtl_',num2str(Pr),'/',diffEq,'_Trajectory_Ra_',num2str(Ra(i)),'_Ro_',num2str(Ro(j)),'.mat');
            
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
if(smallRay==2)
    name = strcat('../SavedData/HeatTransport/',diffEq,'_HeatTransport_FineRes.mat');
elseif(smallRay==1)
    name = strcat('../SavedData/HeatTransport/',diffEq,'_HeatTransport_SmallRay.mat');
else
    name = strcat('../SavedData/HeatTransport/',diffEq,'_HeatTransport.mat');
end
save(name,'M','RaGrid','Pr','k1','RoGrid','heatTransportGrid','limitStdDevPercentGrid');
