function [t,X,runTime] = FluidSolver(varargin)
%FLUIDSOLVER This function integrates the HKC models using ode45
%   Note Tf MUST be a multiple of tInc


addpath('../ModelConstruction');

parameters = varargin{1};
Xprev = varargin{2};
addTime = varargin{3};
saveData = varargin{4};

% parameters

Pr = parameters(1);
Ro = parameters(2);
k1 = parameters(3); V = pi/sqrt(2*k1);
Ra = parameters(4); 
tInc=parameters(5);
Tf = parameters(6);
M = parameters(7);
runTime = parameters(8);

% generate variables
shear = 'yes';

% model name definition
if(strcmp(shear,'no'))
    diffEq = strcat('HKC',num2str(M),'NoShear');
else
    diffEq = strcat('HKC',num2str(M));
end

% initial conditions

if(Tf == 0)
    X0 = InitialConditions(M,shear,k1);
else
    X0=Xprev(end,:);
end

n=size(X0,2);

% integrate
tspan=0:tInc:addTime;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));

if(runTime >= 0)
    tic;
end

odeFunc = @(t,X) feval(diffEq,X,Pr,Ra,Ro,k1,V);
[t,X]=ode45(odeFunc,tspan,X0,options);

if(runTime >= 0)
    runTime = runTime + toc
end


if(Tf > 0)
    tOld = tInc*(0:1:round(Tf/tInc))'; % The colon operator may fail to include the last index if the quantities are not integers.
    t = [ tOld ; (t(2:end)+Tf) ];
    X = [ Xprev ; X(2:end,:)];
end

Tf = Tf + addTime;
name = strcat('../SavedData/Trajectories_Prandtl_1/',diffEq,'_Trajectory_Ra_',num2str(Ra),'_Ro_',num2str(Ro),'.mat');

if(nargin == 4 && saveData)
    save(name,'Ra','Pr','k1','Ro','Tf','tInc','M','runTime','X');
elseif(nargin == 6 && saveData)
    heatTransport = varargin{5};
    limitStdDevPercent = varargin{6};
    save(name,'Ra','Pr','k1','Ro','Tf','tInc','M','runTime','X','heatTransport','limitStdDevPercent');
elseif(nargin ~= 6 && nargin ~= 4)
    error('Wrong number of input arguments!');
end


end

