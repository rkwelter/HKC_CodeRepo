function [uPlusModes,uMinusModes,thetaModes] = VariableConstructor(M,shear)
%VARIABLECONSTRUCTOR Summary of this function goes here
%   Detailed explanation goes here

%M = model number > 0, full shell if M = n(n+1)/2 for some n

n=0;

while(M > (n+1)*(n+2)/2)
    n=n+1;
end

%% Adjoin common modes

commonModes = ones(M,2);
m1 = 1; m3 = 1;
M=M-1;

while(M > 0)
    if(m3 == 1)
        m3 = m1+1;
        m1 = 1;
        commonModes(M,:) = [m1,m3];
        M = M-1;
    else
        m1 = m1 + 1;
        m3 = m3-1;
        commonModes(M,:) = [m1,m3];
        M = M-1;
    end
end

commonModes=sortrows(commonModes);
uPlusModes = commonModes;
uMinusModes = commonModes;
thetaModes = commonModes;

%% Adjoin shear and stratified modes

strat = [0,2];
shearUp = [0,1];
shearUm = [0,1];

if(n>0)
    evens = (2*(1:n)+2)';
    strat = [strat; zeros(n,1), evens];
end

thetaModes= [strat; thetaModes];

if(strcmp(shear,'yes'))
    if(n>0)
        odds = 2*(1:n)'+1;
        uPlusModes= [zeros(n,1), odds;uPlusModes];
        uMinusModes= [zeros(n,1), odds; (1:n)', zeros(n,1); uMinusModes];
    end
    uPlusModes= [shearUp; uPlusModes];
    uMinusModes= sortrows([shearUm; uMinusModes]);
end

end

