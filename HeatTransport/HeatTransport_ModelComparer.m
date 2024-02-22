% Compute the empirical Nusselt dimension and compare to the dimension of
% the minimal HKC model containing the unstable manifold 
% by Roland Welter

clear all, close all, clc
addpath('../ModelConstruction');
addpath('../FluidSolver');
date = '2024_01_19_';

% select parameters

Pr = 10;
k1 = 1/sqrt(2);
V = pi/sqrt(2*k1);
Ra = [1 (1:1:10)*50 (6:1:10)*10^2 (2:1:5)*10^3 ]; % small, fine Rayleigh grid
Ro = 0:50:300;

% build grid

RaGrid = ones(length(Ra),length(Ro));
RoGrid = ones(length(Ra),length(Ro));
nusseltDim = zeros(size(RaGrid));

for i =1:length(Ra)
    for j =1:length(Ro)
        RaGrid(i,j) = Ra(i);
        RoGrid(i,j) = Ro(j);
    end
end

% load heat transport data
heatTransport = cell(1,14);

for maxShell = 1:14
    M = maxShell*(maxShell+1)/2;
    diffEq = strcat('HKC',num2str(M));
    name = strcat('../SavedData/HeatTransport/',diffEq,'_HeatTransport_SmallRay.mat');
    allData = load(name);
    heatTransport{maxShell} = allData.heatTransportGrid;
end

for i=1:length(Ra)
    for j =1:length(Ro)
        for maxShell = 1:13
            M = maxShell*(maxShell+1)/2;
            if(nusseltDim(i,j) == 0 && (heatTransport{maxShell}(i,j) >= .9*heatTransport{14}(i,j) ) && (heatTransport{maxShell}(i,j) <= 1.1*heatTransport{14}(i,j) ) )
                nusseltDim(i,j) = (3*M+4*floor((sqrt(8*M+1)-1)/2)-1);
            end
        end
    end
end

surf(RoGrid,RaGrid,nusseltDim,'FaceAlpha',.75);
xlim([0 300]); ylim([0 5000]);
hold on;

name = strcat('../SavedData/HeatTransport/UnstableManifoldAmbientDimension.mat');
unstableManifoldDim = load(name);

% surf(RoGrid,RaGrid,320*pi^3/(Pr*(1+Pr)*min(1,k1^2))*(1+RaGrid)); % Hausdorff dimension bound
%surf(RoGrid,RaGrid,floor(sqrt(2/27)*RaGrid.^(1/2)./(k1*(1+RoGrid.^2./RaGrid))),'FaceColor','r'); % Dimension of the unstable manifold of the origin
surf(RoGrid,RaGrid,unstableManifoldDim.unstableDim,'FaceColor','r','FaceAlpha',.5); 

name = '../SavedData/HeatTransport/HeatTransportModelComparison.mat';
save(name,'M','RaGrid','Pr','k1','RoGrid','nusseltDim');
