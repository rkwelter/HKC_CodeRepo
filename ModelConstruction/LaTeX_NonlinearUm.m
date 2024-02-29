function [nonlinear] = LaTeX_NonlinearUm(i,triples,uPlusModes,uMinusModes)
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
  else
    prefactor = strcat('2^{\frac{',num2str(zeroCount),'}{2}}');
  end

  %Coefficient * Nonlinearity
  if(coef ~= 0)
    KpSq = AnisoLaplacian_LaTeX(p);
    if(p(1) == 0)
        nonlinear = strcat(nonlinear,' \frac{',prefactor,'(',num2str(coef),')}{',num2str(p(2)),'} u^+_{(',num2str(p(1)),',',num2str(p(2)),')} u^-_{(',num2str(q(1)),',',num2str(q(2)),')} ');
    else
         nonlinear = strcat(nonlinear,' \frac{',prefactor,'(',num2str(coef),')}{\sqrt{',KpSq,'}} u^+_{(',num2str(p(1)),',',num2str(p(2)),')} u^-_{(',num2str(q(1)),',',num2str(q(2)),')} ');
    end

    %End factor
    if(i < numTrips)
      nonlinear = strcat(nonlinear,' + '); %if the last terms are zero then this might end with a plus +
    end
  end

end

end

