function [] = FluidSimulator_All(parameters,X,numPhotos,Nt,reso,Nv,plotInc)
% This function takes a trajectories from the HKC model and plots the
% associated temperature and velocity fields as a function of time, as well
% as the plotting the trajectory in phase space.

%clear all, close all, clc
addpath('../ModelConstruction');
addpath('../FluidSolver');

date = '2024_01_24_';

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

[Xtgrid,Ztgrid] = meshgrid((1:Nt)/Nt,(Nt:-1:1)/Nt);
Ztgrid = pi*Ztgrid;
Xtgrid = 2*pi*Xtgrid/k1;

[Xvgrid,Zvgrid] = meshgrid((1:Nv)/Nv,(Nv:-1:1)/Nv);
Zvgrid = pi*Zvgrid;
Xvgrid = 2*pi*Xvgrid/k1;

if(Ra >= ((1+k1^2)^3+Ro^2)/k1^2)
    Xf = [2*pi*(1+k1^2)*sqrt(k1^2*Ra/((1+k1^2)^3+Ro^2)-1)/k1^(3/2);
          2*pi*((1+k1^2)^3+Ro^2)*sqrt(k1^2*Ra/((1+k1^2)^3+Ro^2)-1)/(k1^(5/2)*sqrt(1+k1^2)*Ra);
          pi*(k1^2*Ra-(1+k1^2)^3-Ro^2)/(k1^(5/2)*Ra)];
else
    Xf=[0;0;0];
end

% Plot fluid flow

for j=1:numPhotos
    
    current = j*plotInc + 7000;

    f=figure;
    f.Position(3:4) = [1500 600]; %[900 413];

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

    %sfh1 = subplot(1,2,1);
    %sfh1.Position(3:4) = [900 600];

    subplot(1,2,1)
    Temp = 1 - Ztgrid/pi;
    
    for i=1:numTheta
        m = thetaModes(i,:);
        if( m(1)~=0 && m(2)~=0)
            cm=1;
        else
            cm=1/sqrt(2);
        end

        if(mod(m(1)+m(2),2)==0)
            Temp = Temp + cm*X(current,numUp+numUm+i)*cos(k1*m(1)*Xtgrid).*sin(m(2)*Ztgrid)/(V*pi);
        else
            Temp = Temp + cm*X(current,numUp+numUm+i)*sin(k1*m(1)*Xtgrid).*sin(m(2)*Ztgrid)/(V*pi);
        end
    end
    
    maxTemp = max(max(max(Temp)),1);
    minTemp = min(min(min(Temp)),0);

    contourf(Xtgrid,Ztgrid,Temp,'LevelList',(maxTemp-minTemp)*(((-reso:1:reso)/reso).^3/2+(-reso:1:reso)/(2*reso)+1)/2+minTemp,'LineColor','none');
    colormap jet;
    colorbar;
    set(gca,'FontSize',16)
    %title(strcat('Ra = ',num2str(Ra),' , Ro = ',num2str(Ro),' , t = ',num2str(tInc*current)));
    %pause(.5)

    %sfh2 = subplot(1,2,2);
    %sfh2.Position(3:4) = [900 600];

    subplot(1,2,2)

    Velocity1 = zeros(Nv,Nv);
    %Velocity2 = zeros(N,Ny,N);
    Velocity3 = zeros(Nv,Nv);
    
    for i=1:numUp
        m = uPlusModes(i,:);
        if( m(1)~=0 && m(2)~=0)
            cm=1;
        else
            cm=1/sqrt(2);
        end

        if(mod(m(1)+m(2),2)==0)
            Velocity1 = Velocity1 + cm*X(current,i)*m(2)*sin(k1*m(1)*Xvgrid).*cos(m(2)*Zvgrid)/(V*(k1*m(1)+m(2)));
            Velocity3 = Velocity3 - cm*X(current,i)*k1*m(1)*cos(k1*m(1)*Xvgrid).*sin(m(2)*Zvgrid)/(V*(k1*m(1)+m(2)));
        else
            Velocity1 = Velocity1 + cm*X(current,i)*m(2)*cos(k1*m(1)*Xvgrid).*cos(m(2)*Zvgrid)/(V*(k1*m(1)+m(2)));
            Velocity3 = Velocity3 + cm*X(current,i)*k1*m(1)*sin(k1*m(1)*Xvgrid).*sin(m(2)*Zvgrid)/(V*(k1*m(1)+m(2)));
        end
    end
    
    %maxTemp = max(max(max(Temp)),1);
    %minTemp = min(min(min(Temp)),0);

    quiver(Xvgrid,Zvgrid,Velocity1,Velocity3);
    %axis equal
    xlabel('X','FontSize',16);
    zlabel('Z','FontSize',16);
    set(gca,'FontSize',16)
    %view(-10,10);

    name = strcat('../SavedImages/pics/',date,diffEq,'_Temp_Ra_',num2str(Ra),'_Ro_',num2str(Ro),'_',num2str(j),'.jpg');
    saveas(gcf,name);
    %close;

end


end