function [nonlinear] = Mathematica_NonlinearUm(i,triples,uPlusModes,uMinusModes)
%NONLINEARUM Summary of this function goes here
%   Detailed explanation goes here

nonlinear = '';
numTrips = size(triples,2);
m = uMinusModes(i,:);

for i=1:numTrips

  p = uPlusModes(triples(1,i),:);
  q = uMinusModes(triples(2,i),:);
  indexCond = triples(3,i);

  zeroCount = ZeroCounter(m(1))+ZeroCounter(m(2))+ZeroCounter(p(1))+ZeroCounter(p(2))+ZeroCounter(q(1))+ZeroCounter(q(2));
  [S1,S2,S3,S4] = SignCoefficients(m,q,indexCond);
  coef = (q(1)*p(2)*S1-p(1)*q(2)*S3);

  %Prefactor
  if(zeroCount==0)
    prefactor = '';
  elseif(mod(zeroCount,2)==0)
    prefactor = strcat(num2str(2^(zeroCount/2)),'*');
  else
    prefactor = strcat(num2str(2^((zeroCount-1)/2)),'*Sqrt[2]*');
  end

  %Coefficient * Nonlinearity
  if(coef ~= 0)
    KpSq = AnisoLaplacian(p);
    nonlinear = strcat(nonlinear,prefactor,'(',num2str(coef),')*X[[',num2str(triples(1,i)),']]*X[[',num2str(size(uPlusModes,1)+triples(2,i)),']]');

    if( p(1) >0 || p(2) > 1 )
        nonlinear = strcat(nonlinear,'/Sqrt[',KpSq,']');
    else
    end

    %End factor
    if(i < numTrips)
      nonlinear = strcat(nonlinear,' + '); %if the last terms are zero then this might end with a plus +
    end
  end

end

end

