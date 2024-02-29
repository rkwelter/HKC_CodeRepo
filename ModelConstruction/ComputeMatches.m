function [uPlusMatches,uMinusMatches,thetaMatches] = ComputeMatches(uPlusModes,uMinusModes,thetaModes)
%COMPUTEMATCHES This function creates the libraries uPlusMatches, uMinusMatches, thetaMatches
% which enable one to easily determine across the uPlus, uMinus, theta
% variables which element in the combined vector X corresponds to the same
% Fourier index. 

uPlusMatches = -1*ones(size(uPlusModes,1),2);
uMinusMatches = -1*ones(size(uMinusModes,1),1);
thetaMatches = -1*ones(size(thetaModes,1),1);

for i=1:size(uPlusModes,1)
    m=uPlusModes(i,:);
    for j=1:size(uMinusModes,1)
        mt = uMinusModes(j,:);
        if(m(1)==mt(1)&&m(2)==mt(2))
            uPlusMatches(i,1) = j;
            uMinusMatches(j,1) = i;
        end
    end

    for j=1:size(thetaModes,1)
        mt = thetaModes(j,:);
        if(m(1)==mt(1)&&m(2)==mt(2))
            uPlusMatches(i,2) = j;
            thetaMatches(j,1) = i;
        end
    end
end

end

