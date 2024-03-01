function [] = FluidVis_PhaseSpace(parameters,X,numPhotos,plotInc)
%FLUIDVIS_PHASESPACE by Roland Welter
% 
% Given a solution X from a HKC model, this function plots the solution as
% a trajectory in phase space.  Specifically, the solution is projected 
% onto its Lorenz triple u^+_{(1,1)},\theta_{(1,1)},\theta_{(0,2)}.  The
% entire trajectory is plotted using color_line3, and then snapshots are
% taken with the current state of the trajectory depicted as a point
% travelling along the trajectory.
%
% Required input arguments:
%     - parameters: 8x1 vector, specifying the Prandtl number, rotation
%                   number, shape parameter k1, Rayleigh number,
%                   integration step size, final time, HKC model number,
%                   and cumulative run time.  Note Tf MUST be a multiple of
%                   tInc
%     - X: a matrix storing the trajectory, in which the different columns
%          represent the different variables in the HKC model, and the rows
%          represent the solution at different times.
%     - numPhotos: a positive integer specifying how many snapshots should
%                  be generated.
%     - plotInc: a positive integer specifying the increment between
%                snapshots.  Namely, the time increment between the values 
%                of the solution may be small, hence the user may wish to 
%                view only every 5th time increment.  This depends on the 
%                Rayleigh number, but I typically found that for low
%                Rayleigh numbers and a time increment of 0.001 a good
%                value of plotInc would be 4 or 5.
% 

%clear all, close all, clc
addpath('../ModelConstruction');
addpath('../FluidSolver');

date = '2024_01_24_';

% parameters

Ro = parameters(2);
k1 = parameters(3);
Ra = parameters(4); 
tInc=parameters(5);
Tf = parameters(6); t = tInc:tInc:Tf;
M = parameters(7);

%fracSlow = .65;
%plotInc = Tf/(timeInc*numPhotos);

% model name definition
diffEq = strcat('HKC',num2str(M));

[uPlusModes,uMinusModes,~] = VariableConstructor(M);
numUp = size(uPlusModes,1);
numUm = size(uMinusModes,1);

% Find Lorenz triple

numVert=1;

while(M > (numVert)*(numVert+1)/2)
    numVert=numVert+1;
end

% Fixed points in the Lorenz-Stenflo model.

%if(Ra >= ((1+k1^2)^3+Ro^2)/k1^2)
%    Xf = [2*pi*(1+k1^2)*sqrt(k1^2*Ra/((1+k1^2)^3+Ro^2)-1)/k1^(3/2);
%          2*pi*((1+k1^2)^3+Ro^2)*sqrt(k1^2*Ra/((1+k1^2)^3+Ro^2)-1)/(k1^(5/2)*sqrt(1+k1^2)*Ra);
%          pi*(k1^2*Ra-(1+k1^2)^3-Ro^2)/(k1^(5/2)*Ra)];
%else
%    Xf=[0;0;0];
%end

% Plot fluid flow

for j=1:numPhotos
    
    current = j*plotInc + 7000;

    f=figure;
    f.Position(3:4) = [900 600]; %[900 413];

    % plot3([0,Xf(1),-Xf(1)],[0,Xf(2),-Xf(2)],[0,Xf(3),Xf(3)],'o','Color','b','MarkerSize',7,'MarkerFaceColor','#D9FFFF')
    % hold on;
    color_line3(X(:,numVert+1),-1*X(:,numUp+numUm+numVert+1),-X(:,numUp+numUm+1),t,'LineWidth',1.5); hold on;
    plot3([X(current,numVert+1)],-X(current,numUp+numUm+numVert+1),-X(current,numUp+numUm+1),'o','Color','b','MarkerSize',7,'MarkerFaceColor','#D9FFFF')
    view(27,16)
    grid on
    xlabel('x','FontSize',13)
    ylabel('y','FontSize',13)
    zlabel('z','FontSize',13)
    set(gca,'FontSize',13)
    grid on;

    name = strcat('../SavedImages/pics/',date,diffEq,'_Phase_Ra_',num2str(Ra),'_Ro_',num2str(Ro),'_',num2str(j),'.jpg');
    saveas(gcf,name);
    close;

end


end