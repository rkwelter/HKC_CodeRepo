function [t,X,parameters] = FluidSolver(varargin)
%FLUIDSOLVER by Roland Welter
% 
% Solves and saves a trajectory from a user-selected HKC model
% by calling ode45 with its corresponding time stepper.  Note that this
% function can either initialize a new trajectory with a user-selected
% initial condition or continue an existing trajectory
% 
% Required input arguments:
%     - parameters: 8x1 vector, specifying the Prandtl number, rotation
%                   number, shape parameter k1, Rayleigh number,
%                   integration step size, final time, HKC model number,
%                   and cumulative run time.  Note Tf MUST be a multiple of
%                   tInc
%     - Xprev: a matrix storing the previously computed values of the
%              trajectory.  Should be an empty matrix if the trajectory has
%              not been previously initialized.  Otherwise the different
%              columns represent the different variables in the HKC model,
%              and the rows represent the solution at different times.
%     - addTime: a scalar specifying how much time to add onto the
%                trajectory. Note addTime MUST be a multiple of tInc
%     - saveData: 1 if you would like to automatically save the trajectory
%                 to the "SavedData" folder, 0 if not.
% 
% Optional input arguments:
%     - heatTransport: a vector specifying the running average heat
%                      transport
%     - limitStdDevPercent: a vector specifying the running standard
%                           deviation for the heat transport
%
% Output arguments:
%     - t: a vector storing the time values for the trajectory
%     - X: a matrix storing the trajectory, in which the different columns
%          represent the different variables in the HKC model, and the rows
%          represent the solution at different times.
%     - parameters: same as above, but with updated values for Tf, runTime

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

% model name definition
diffEq = strcat('HKC',num2str(M));

% initial conditions

if(Tf == 0)
    X0 = InitialConditions(M,k1);
else
    X0=Xprev(end,:);
end

n=size(X0,2);

% integrate
tspan=0:tInc:addTime;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));

if(runTime >= 0) % uncomment here and below to test the computational time
    tic;
end

odeFunc = @(t,X) feval(diffEq,X,Pr,Ra,Ro,k1,V);
[t,X]=ode45(odeFunc,tspan,X0,options);

if(runTime >= 0)
    runTime = runTime + toc;
end


if(Tf > 0)
    tOld = tInc*(0:1:round(Tf/tInc))'; % The colon operator may fail to include the last index if the quantities are not integers.
    t = [ tOld ; (t(2:end)+Tf) ];
    X = [ Xprev ; X(2:end,:)];
end

Tf = Tf + addTime;
parameters(6) = Tf;
parameters(8) = runTime;
folderName = strcat('../SavedData/Trajectories/Prandtl_',num2str(Pr));
fileName = strcat(folderName,'/',diffEq,'_Trajectory_Ra_',num2str(Ra),'_Ro_',num2str(Ro),'.mat');

if(nargin == 4 && saveData)
    if(isfolder(folderName))
        save(fileName,'Ra','Pr','k1','Ro','Tf','tInc','M','runTime','X');
    elseif(mkdir(folderName))
        save(fileName,'Ra','Pr','k1','Ro','Tf','tInc','M','runTime','X');
    else
        error('Folder creation failed.');
    end
elseif(nargin == 6 && saveData)
    heatTransport = varargin{5};
    limitStdDevPercent = varargin{6};
    if(isfolder(folderName))
        save(fileName,'Ra','Pr','k1','Ro','Tf','tInc','M','runTime','X','heatTransport','limitStdDevPercent');
    elseif(mkdir(folderName))
        save(fileName,'Ra','Pr','k1','Ro','Tf','tInc','M','runTime','X','heatTransport','limitStdDevPercent');
    else
        error('Folder creation failed.');
    end
elseif(nargin ~= 6 && nargin ~= 4)
    error('Wrong number of input arguments!');
end


end

