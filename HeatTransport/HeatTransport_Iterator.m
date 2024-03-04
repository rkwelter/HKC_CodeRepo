%% HeatTransport_Iterator by Roland Welter
% 
% Computes the heat transport for a large number of Rayleigh and rotation
% numbers by creating trajectories for each.

%clear all, close all, clc
addpath('../ModelConstruction');
addpath('../FluidSolver'); addpath('../FluidSolver/TimeSteppers/');

%% generate Data 

% select parameters

M=190;
Pr = 10;
k1 = 1/sqrt(2); V = pi/sqrt(2*k1);
RaVec = 7500;
%RaVec = [1 (1:1:4)*10^3/4 (1:1:4)*10^4/4 (1:1:4)*10^5/4 (3:1:20)*10^5/2 ]; % Large, course Rayleigh grid
%RaVec = [1 (1:1:10)*50 (6:1:10)*10^2 (2:1:5)*10^3 ]; % small, fine Rayleigh grid
%RaVec = [1 , 5:5:2000]; % hyperfine Rayleigh grid
%RoVec = (0:1:10)*10^2; % Large, course rotation grid
%RoVec = 0:50:300; % Small, fine rotation grid
RoVec = 0; % no rotation
errTol = 10^(-7);

numRuns = 1; % numRuns controls how many different trajectories are computed (ie different init. cond.) for each Rayleigh/rotation number
extendLimit = 20; % extendLimit sets a maximum number of iterations that the solver continues a trajectory if the variation in the heat transport remains large

% generate variables

[uPlusModes,uMinusModes,thetaModes] = VariableConstructor(M);
numUp = size(uPlusModes,1); 
numUm = size(uMinusModes,1); 
numTheta = size(thetaModes,1);
n=numUp+numUm+numTheta;
diffEq = strcat('HKC',num2str(M));

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

extensions = 0;
progress = 0;
totaljobs = length(RaVec)*length(RoVec)*numRuns;

for i=1:length(RaVec)
    for j =1:length(RoVec)
        for k=1:numRuns

            Ra = RaVec(i);
            Ro = RoVec(j);
            name = strcat('../SavedData/Trajectories/Prandtl_',num2str(Pr),'/',diffEq,'_Trajectory_Ra_',num2str(Ra),'_Ro_',num2str(Ro),'.mat');
            
            if(isfile(name))
                allData = load(name);
                X = allData.X;
                parameters = allData.parameters;
                heatTransport = allData.heatTransport;
                limitStdDevPercent = allData.limitStdDevPercent;
                addTime = 500*tInc;
            else
                
                Up0 = 0.05*randn(numUp,1)./(uPlusModes(:,1).^2+uPlusModes(:,2).^2);
                Um0 = 0.05*randn(numUm,1)./(uMinusModes(:,1).^2+uMinusModes(:,2).^2);
                The0 = The0nearUni + 0.05*randn(numTheta,1)./(thetaModes(:,1).^2+thetaModes(:,2).^2);
                
                X=[Up0;Um0;The0]';
                tInc=10^(-5); %min(0.0001,10/(RaVec(i)));
                Tf = 0;
                runTime = 0;
                heatTransport = [];
                limitStdDevPercent = 1;
                addTime = 500*tInc;
                parameters = [Pr Ro k1 Ra tInc Tf M runTime errTol];
            end 


            while( extensions < extendLimit ) % limitStdDevPercent > .01 &&
                %addTime = 10^4*(extendLimit-extensions)*tInc;
                [t,X,parameters] = FluidSolver(parameters,X,addTime,0);

                % Compute Nusselt integrand
            
                NusseltIntegrand = -sqrt(k1)*thetaModes(1,2)*X(:,numUp+numUm+1)/(pi^2);
                if(numVert > 1)
                    for ell=2:numVert
                        NusseltIntegrand = NusseltIntegrand -sqrt(k1)*thetaModes(ell,2)*X(:,numUp+numUm+ell)/(pi^2);
                    end
                end

                heatTransport = 1+tInc*cumsum(NusseltIntegrand)./t;
                limitStdDevPercent = std(heatTransport(floor(end/2):end))/heatTransport(end);
                extensions = extensions + 1
                %Tf = Tf + addTime;
                addTime = 500*tInc;
            end
                            
            if(extensions > 0)
                save(name,'parameters','X','heatTransport','limitStdDevPercent');
            end
            extensions = 0;
        end
        progress=progress+1.0;
        progress/totaljobs
    end
end
