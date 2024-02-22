function [] = HeatTransport_Trajectory(parameters,X)
% Given a trajectory this function computes its corresponding heat
% transport

addpath('../ModelConstruction');
date = '2024_01_18_';

% parameters

Pr = parameters(1);
Ro = parameters(2);
k1 = parameters(3); V = pi/sqrt(2*k1);
Ra = parameters(4); 
tInc=parameters(5);
Tf = parameters(6); t = (tInc:tInc:Tf)';
M = parameters(7);
runTime = parameters(8);
%heatTransport = parameters(9);
%limitStdDevPercent = parameters(10);

% generate variables
M=78;
shear = 'yes';
[uPlusModes,uMinusModes,thetaModes] = VariableConstructor(M,shear);
numUp = size(uPlusModes,1);
numUm = size(uMinusModes,1);
numTheta = size(thetaModes,1);

% Find Lorenz triple

numVert=1;

while(M > (numVert)*(numVert+1)/2)
    numVert=numVert+1;
end

%% Heat transport

% Nusselt integrand
            
NusseltIntegrand = -sqrt(k1)*thetaModes(1,2)*X(:,numUp+numUm+1)/(pi^2);
      if(numVert > 1)
         for ell=2:numVert
             NusseltIntegrand = NusseltIntegrand - sqrt(k1)*thetaModes(ell,2)*X(:,numUp+numUm+ell)/(pi^2);
         end
      end

heatTransport = tInc*cumsum(NusseltIntegrand)./t;
limitStdDevPercent = std(heatTransport(floor(end/2):end))/heatTransport(end);

plot(t,heatTransport,'LineWidth',2);
picTitle = strcat('Ra = ', num2str(Ra),' , Ro = ', num2str(Ro));
title(picTitle);
xlim([0 Tf]); ylim([0 ceil(max(heatTransport))]);
ylabel(strcat('Nu^{',num2str(M),'}(t)')); xlabel('t'); 

%name = strcat(date,diffEq,'_HeatTransport_',num2str((i-1)*length(Ro)+j),'.jpg');
%saveas(gcf,name);
%close;

end

