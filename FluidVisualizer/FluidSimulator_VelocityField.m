function [] = FluidSimulator_VelocityField(parameters,X,numPhotos,N,Ny,plotInc)
%FLUIDSIMULATOR_VELOCITYFIELD by Roland Welter
% 
% Given a solution X from a HKC model, this function generates a sequence
% of snapshots of the corresponding velocity field which show its evolution
% as a function of time.  
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
%     - N: a positive integer specifying how many points the velocity field 
%          should be evaluated at.  Specifically, this function evaluates 
%          the velocity on an NxN grid in the (x1,x3) plane. Typically 
%          15-20 is a good choice, since the image becomes very dense for 
%          higher resolution.
%     - Ny: a positive integer specifying how many points in the x2
%           direction that the velocity field should be evaluated at.
%           Typically 5-8 is a good choice, again because the image becomes
%           very dense for higher resolution.
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

date = '2024_01_18_';

% parameters

Ro = parameters(2);
k1 = parameters(3); V = pi/sqrt(2*k1);
Ra = parameters(4); 
tInc=parameters(5);
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

% Build grid

[Xgrid,Ygrid,Zgrid] = ndgrid((1:N)/N,(1:Ny)/Ny,(N:-1:1)/N);
Zgrid = pi*Zgrid;
Xgrid = 2*pi*Xgrid/k1;
Ygrid = pi*Ygrid/k1;

% Plot fluid flow

for j=1:numPhotos
    
    current = j*plotInc;

    f=figure;
    f.Position(3:4) = [1350 600];

    Velocity1 = zeros(N,Ny,N);
    Velocity2 = zeros(N,Ny,N);
    Velocity3 = zeros(N,Ny,N);

    for i=1:numUp
        m = uPlusModes(i,:);
        if( m(1)~=0 && m(2)~=0)
            cm=1;
        else
            cm=1/sqrt(2);
        end

        if(mod(m(1)+m(2),2)==0)
            Velocity1 = Velocity1 + cm*X(current,i)*m(2)*sin(k1*m(1)*Xgrid).*cos(m(2)*Zgrid)/(V*(k1*m(1)+m(2)));
            Velocity3 = Velocity3 - cm*X(current,i)*k1*m(1)*cos(k1*m(1)*Xgrid).*sin(m(2)*Zgrid)/(V*(k1*m(1)+m(2)));
        else
            Velocity1 = Velocity1 + cm*X(current,i)*m(2)*cos(k1*m(1)*Xgrid).*cos(m(2)*Zgrid)/(V*(k1*m(1)+m(2)));
            Velocity3 = Velocity3 + cm*X(current,i)*k1*m(1)*sin(k1*m(1)*Xgrid).*sin(m(2)*Zgrid)/(V*(k1*m(1)+m(2)));
        end
    end

    for i=1:numUm
        m = uMinusModes(i,:);
        if( m(1)~=0 && m(2)~=0)
            cm=1;
        else
            cm=1/sqrt(2);
        end

        if(mod(m(1)+m(2),2)==0)
            Velocity2 = Velocity2 + cm*X(current,numUp+i)*sin(k1*m(1)*Xgrid).*cos(m(2)*Zgrid)/(V);
        else
            Velocity2 = Velocity2 + cm*X(current,numUp+i)*cos(k1*m(1)*Xgrid).*cos(m(2)*Zgrid)/(V);
        end
    end

    quiver3(Xgrid,Ygrid,Zgrid,Velocity1,Velocity2,Velocity3);
    axis equal
    xlabel('X','FontSize',13);
    ylabel('Y','FontSize',13);
    zlabel('Z','FontSize',13);
    view(-10,10);
    title(strcat('Ra = ',num2str(Ra),' , Ro = ',num2str(Ro),' , t = ',num2str(tInc*current)));
    %pause(.5)
    name = strcat('../SavedImages/pics/',date,diffEq,'_Velocity_Ra_',num2str(Ra),'_Ro_',num2str(Ro),'_',num2str(j),'.jpg');
    saveas(gcf,name);
    close;

end


end
