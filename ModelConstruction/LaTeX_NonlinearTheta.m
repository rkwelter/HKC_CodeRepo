function [nonlinear] = LaTeX_NonlinearTheta(i,triples,uPlusModes,thetaModes,numUm)
%NONLINEARTHETA Summary of this function goes here
%   Detailed explanation goes here

nonlinear = '';
numTrips = size(triples,2);
m = thetaModes(i,:);

for i=1:numTrips

  p = uPlusModes(triples(1,i),:);
  q = thetaModes(triples(2,i),:);
  indexCond = triples(3,i);

  zeroCount = ZeroCounter(m(1))+ZeroCounter(m(2))+ZeroCounter(p(1))+ZeroCounter(p(2))+ZeroCounter(q(1))+ZeroCounter(q(2));
  [S1,S2,S3,S4] = SignCoefficients(m,q,indexCond);
  coef = (q(1)*p(2)*S2-p(1)*q(2)*S4);

  %Prefactor
  if(mod(q(1)+q(2)+m(1)+m(2),2)==1)
    prefactor = strcat('(-1)');
  else
    prefactor = '';
  end

  if(zeroCount~=0)
    prefactor = strcat(prefactor,'2^{\frac{',num2str(zeroCount),'}{2}}');
  end

  %Coefficient * Nonlinearity
  if(coef ~= 0)
    KpSq = AnisoLaplacian_LaTeX(p);
    if(p(1) == 0)
        nonlinear = strcat(nonlinear,' \frac{',prefactor,'(',num2str(coef),')}{',num2str(p(2)),'} u^+_{(',num2str(p(1)),',',num2str(p(2)),')} \theta_{(',num2str(q(1)),',',num2str(q(2)),')} ');
    else
         nonlinear = strcat(nonlinear,' \frac{',prefactor,'(',num2str(coef),')}{\sqrt{',KpSq,'}} u^+_{(',num2str(p(1)),',',num2str(p(2)),')} \theta_{(',num2str(q(1)),',',num2str(q(2)),')} ');
    end

    %End factor
    if(i < numTrips)
      nonlinear = strcat(nonlinear,' + '); %if the last terms are zero then this might end with a plus +
    end
  end

end

end

