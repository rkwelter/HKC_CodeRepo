function [] = FluidVis_TemperatureField(parameters,X,numPhotos,N,reso,plotInc)
%FLUIDVIS_TEMPERATUREFIELD by Roland Welter
% 
% Given a solution X from a HKC model, this function generates a sequence
% of snapshots of the corresponding temperature field which show its 
% evolution as a function of time.  
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
%     - N: a positive integer specifying how many points the temperature
%          field should be evaluated at.  Specifically, this function 
%          evaluates the temperature on an NxN grid. Typically 100 is a 
%          good choice, but for high Rayleigh number flows the user may 
%          want to include more to capture the small scale structure.
%     - reso: a positive integer specifying the number of contours used in 
%             the contour plot, ie the resolution of the image.  Typically 
%             40 is a good choice, but again for high Rayleigh number flows
%             the user may want to include more.
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

date = '2024_01_29_';

% parameters

Ro = parameters(2);
k1 = parameters(3); V = pi/sqrt(2*k1);
Ra = parameters(4); 
tInc=parameters(5);
%Tf = parameters(6);
M = parameters(7);

%fracSlow = .65;
%plotInc = Tf/(timeInc*numPhotos);

% model name definition
diffEq = strcat('HKC',num2str(M));

[uPlusModes,uMinusModes,thetaModes] = VariableConstructor(M);
numUp = size(uPlusModes,1);
numUm = size(uMinusModes,1);
numTheta = size(thetaModes,1);

% Find Lorenz triple

numVert=1;

while(M > (numVert)*(numVert+1)/2)
    numVert=numVert+1;
end

% Build grid

[Xgrid,Zgrid] = meshgrid((1:N)/N,(N:-1:1)/N);
Zgrid = pi*Zgrid;
Xgrid = 2*pi*Xgrid/k1;

% Plot fluid flow

for j=1:numPhotos
    
    current = j*plotInc;

    f=figure;
    f.Position(3:4) = [900 600]; %[900 413];

    Temp = 1 - Zgrid/pi;
    
    for i=1:numTheta
        m = thetaModes(i,:);
        if( m(1)~=0 && m(2)~=0)
            cm=1;
        else
            cm=1/sqrt(2);
        end

        if(mod(m(1)+m(2),2)==0)
            Temp = Temp + cm*X(current,numUp+numUm+i)*cos(k1*m(1)*Xgrid).*sin(m(2)*Zgrid)/(V*pi);
        else
            Temp = Temp + cm*X(current,numUp+numUm+i)*sin(k1*m(1)*Xgrid).*sin(m(2)*Zgrid)/(V*pi);
        end
    end
    
    maxTemp = max(max(max(Temp)),1);
    minTemp = min(min(min(Temp)),0);

    contourf(Xgrid,Zgrid,Temp,'LevelList',(maxTemp-minTemp)*(((-reso:1:reso)/reso).^3/2+(-reso:1:reso)/(2*reso)+1)/2+minTemp,'LineColor','none');
    colormap jet;
    colorbar;
    set(gca,'FontSize',16)
    title(strcat('Ra = ',num2str(Ra),' , Ro = ',num2str(Ro),' , t = ',num2str(tInc*current)));
    %pause(.5)
    name = strcat('../SavedImages/pics/',date,diffEq,'_Temp_Ra_',num2str(Ra),'_Ro_',num2str(Ro),'_',num2str(j),'.jpg');
    saveas(gcf,name);
    close;

end


end
