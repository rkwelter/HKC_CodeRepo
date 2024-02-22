% Simulator and heat transport calculator for the rotating Lorenz system
% by Roland Welter

%clear all, close all, clc
addpath('../ModelConstruction');
addpath('../FluidSolver');
date = '2024_01_26_';

%% generate Data 

% generate variables
M=105;
shear = 'yes';
[uPlusModes,uMinusModes,thetaModes] = VariableConstructor(M,shear);

% select parameters

Pr = 10;
k1 = 1/sqrt(2);
V = pi/sqrt(2*k1);
%RaVec = [1 (1:1:4)*10^3/4 (1:1:4)*10^4/4 (1:1:4)*10^5/4 (3:1:20)*10^5/2 ]; % Large, course Rayleigh grid
%RoVec = (0:1:10)*10^2; % Large, course rotation grid
RaVec = [1 (1:1:10)*50 (6:1:10)*10^2 (2:1:5)*10^3 ]; % small, fine Rayleigh grid
RoVec = 0:50:300; % Small, fine rotation grid
numRuns = 1;
extendLimit = 100;

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

% near uniform IC

numVert=1;

while(M > (numVert)*(numVert+1)/2)
    numVert=numVert+1;
end

The0nearUni = zeros(numTheta,1);

for i=1:numVert
    The0nearUni(i) = The0nearUni(i) - 2*pi/(sqrt(k1)*thetaModes(i,2));
end

% Integrate

options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
extensions = 0;
progress = 0;
totaljobs = length(RaVec)*length(RoVec)*numRuns;

for i=1:length(RaVec)
    for j =1:length(RoVec)
        for k=1:numRuns

            Ra = RaVec(i);
            Ro = RoVec(j);
            name = strcat('../SavedData/Trajectories_Prandtl_',num2str(Pr),'/',diffEq,'_Trajectory_Ra_',num2str(Ra),'_Ro_',num2str(Ro),'.mat');
            
            if(isfile(name))
                allData = load(name);
                X = allData.X;
                tInc = allData.tInc;
                Tf = allData.Tf;
                runTime = allData.runTime;
                heatTransport = allData.heatTransport;
                limitStdDevPercent = allData.limitStdDevPercent;
                addTime = 10^3*tInc;
            else
                Up0 = 0.001*randn(numUp,1);
                Um0 = 0.001*randn(numUm,1);
                The0 = The0nearUni + 0.001*randn(numTheta,1);
                
                X=[Up0;Um0;The0]';
                tInc=min(0.0001,10/(RaVec(i)));
                Tf = 0;
                runTime = 0;
                heatTransport = [];
                limitStdDevPercent = 1;
                addTime = 10^4*tInc;
            end 

            while( limitStdDevPercent > .01 && extensions < extendLimit ) % && size(X,1) < 5*10^4 )
                %addTime = 10^4*(extendLimit-extensions)*tInc;
                [t,X,runTime] = FluidSolver([Pr Ro k1 Ra tInc Tf M runTime],X,addTime,0);

                % Nusselt integrand
            
                NusseltIntegrand = -sqrt(k1)*thetaModes(1,2)*X(:,numUp+numUm+1)/(pi^2);
                if(numVert > 1)
                    for ell=2:numVert
                        NusseltIntegrand = NusseltIntegrand -sqrt(k1)*thetaModes(ell,2)*X(:,numUp+numUm+ell)/(pi^2);
                    end
                end

                heatTransport = 1+tInc*cumsum(NusseltIntegrand)./t;
                limitStdDevPercent = std(heatTransport(floor(end/2):end))/heatTransport(end);
                extensions = extensions + 1;
                Tf = Tf + addTime;
                addTime = 10^3*tInc;
            end
            
            if(extensions > 0)
                save(name,'Ra','Pr','k1','Ro','Tf','tInc','M','runTime','X','heatTransport','limitStdDevPercent');
            end
            extensions = 0;
        end
        progress=progress+1.0;
        progress/totaljobs
    end
end
