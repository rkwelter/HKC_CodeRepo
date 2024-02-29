function [X0] = InitialConditions(M,k1)
%INITIALCONDITIONS Defines an initial condition for a new trajectory

[uPlusModes,uMinusModes,thetaModes] = VariableConstructor(M);

numUp = size(uPlusModes,1);
numUm = size(uMinusModes,1);
numTheta = size(thetaModes,1);

% random IC
Up0 = 0.05*randn(numUp,1);
Um0 = 0.05*randn(numUm,1);
The0 = 0.05*randn(numTheta,1);

% near uniform IC

numVert=1;

while(M > (numVert)*(numVert+1)/2)
    numVert=numVert+1;
end

for i=1:numVert
    The0(i) = The0(i) - 2*pi/(sqrt(k1)*thetaModes(i,2));
end

% homoclinic IC
%Up0 = .01;
%Um0 = 0;
%The0 = [0; .01];

X0=[Up0;Um0;The0];

end

