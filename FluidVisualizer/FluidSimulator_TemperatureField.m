function [] = FluidSimulator_TemperatureField(parameters,X,numPhotos,N,reso,plotInc)
% This function takes a trajectories from the HKC model and plots the
% associated temperature field as a function of time.

%clear all, close all, clc
addpath('../ModelConstruction');
addpath('../FluidSolver');

date = '2024_01_29_';

% parameters

Ro = parameters(2);
k1 = parameters(3); V = pi/sqrt(2*k1);
Ra = parameters(4); 
tInc=parameters(5);
Tf = parameters(6); t = tInc:tInc:Tf;
M = parameters(7);
shear = 'yes';

%fracSlow = .65;
%plotInc = Tf/(timeInc*numPhotos);

% model name definition
if(strcmp(shear,'no'))
    diffEq = strcat('HKC',num2str(M),'NoShear');
else
    diffEq = strcat('HKC',num2str(M));
end

[uPlusModes,uMinusModes,thetaModes] = VariableConstructor(M,shear);
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

if(Ra >= ((1+k1^2)^3+Ro^2)/k1^2)
    Xf = [2*pi*(1+k1^2)*sqrt(k1^2*Ra/((1+k1^2)^3+Ro^2)-1)/k1^(3/2);
          2*pi*((1+k1^2)^3+Ro^2)*sqrt(k1^2*Ra/((1+k1^2)^3+Ro^2)-1)/(k1^(5/2)*sqrt(1+k1^2)*Ra);
          pi*(k1^2*Ra-(1+k1^2)^3-Ro^2)/(k1^(5/2)*Ra)];
else
    Xf=[0;0;0];
end

% Plot fluid flow

for j=1:numPhotos
    
    current = j*plotInc;

    f=figure;
    f.Position(3:4) = [900 600]; %[900 413];

    %subplot(1,2,1)
    %plot3([0,Xf(1),-Xf(1)],[0,Xf(2),-Xf(2)],[0,Xf(3),Xf(3)],'o','Color','b','MarkerSize',7,'MarkerFaceColor','#D9FFFF')
    %hold on;
    %color_line3(X(:,numVert+1),-1*X(:,numUp+numUm+numVert+1),-X(:,numUp+numUm+1),t,'LineWidth',1.5); hold on;
    %plot3([X(current,numVert+1)],[-X(current,numUp+numUm+numVert+1)],[-X(current,numUp+numUm+1)],'o','Color','b','MarkerSize',7,'MarkerFaceColor','#D9FFFF')
    %view(27,16)
    %grid on
    %xlabel('x','FontSize',13)
    %ylabel('y','FontSize',13)
    %zlabel('z','FontSize',13)
    %set(gca,'FontSize',13)
    %grid on

    %subplot(1,2,2)

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
    %title(strcat('Ra = ',num2str(Ra),' , Ro = ',num2str(Ro),' , t = ',num2str(tInc*current)));
    %pause(.5)
    name = strcat('../SavedImages/pics/',date,diffEq,'_Temp_Ra_',num2str(Ra),'_Ro_',num2str(Ro),'_',num2str(j),'.jpg');
    saveas(gcf,name);
    close;

end


end
