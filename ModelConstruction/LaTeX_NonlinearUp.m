function [nonlinear] = LaTeX_NonlinearUp(i,triples,uPlusModes)
%NONLINEARUP Summary of this function goes here
%   Detailed explanation goes here

nonlinear = '';
numTrips = size(triples,2);
m = uPlusModes(i,:);

for i=1:numTrips

  p = uPlusModes(triples(1,i),:);
  q = uPlusModes(triples(2,i),:);
  KpSq = AnisoLaplacian_LaTeX(p);
  KqSq = AnisoLaplacian_LaTeX(q);
  indexCond = triples(3,i);

  zeroCount = ZeroCounter(m(1))+ZeroCounter(m(2))+ZeroCounter(p(1))+ZeroCounter(p(2))+ZeroCounter(q(1))+ZeroCounter(q(2));
  [S1,S2,S3,S4] = SignCoefficients(m,q,indexCond);
  coef0 = (q(1)*p(2)*S1-p(1)*q(2)*S3)*q(2)*m(2);
  coef1 = (q(1)*p(2)*S2-p(1)*q(2)*S4)*q(1)*m(1);

  %Prefactor
  if(zeroCount==0)
    prefactor = '';
  else
    prefactor = strcat('2^{\frac{',num2str(zeroCount),'}{2}}');
  end

  %Coefficient
  if(coef0 ~= 0)
    if(coef1 > 0)
      nonlinear = strcat(nonlinear,' \frac{',prefactor,'(',num2str(coef0),'+',num2str(coef1),'k_1^2)}{\sqrt{(',KpSq,')(',KqSq,')}}');
    elseif(coef1<0)
      nonlinear = strcat(nonlinear,' \frac{',prefactor,'(',num2str(coef0),'-',num2str(-1*coef1),'k_1^2)}{\sqrt{(',KpSq,')(',KqSq,')}}');
    else
      nonlinear = strcat(nonlinear,' \frac{',prefactor,'(',num2str(coef0),')}{\sqrt{(',KpSq,')(',KqSq,')}}');
    end
  else
    if(coef1 ~= 0)
      nonlinear = strcat(nonlinear,' \frac{',prefactor,'(',num2str(coef1),')k_1^2}{\sqrt{(',KpSq,')(',KqSq,')}}');
    end
  end

  %Nonlinearity
  if(coef0 ~= 0 || coef1 ~= 0)
    nonlinear = strcat(nonlinear,'u^+_{(',num2str(p(1)),',',num2str(p(2)),')} u^+_{(',num2str(q(1)),',',num2str(q(2)),')}');

    %End factor
    if(i < numTrips)
      nonlinear = strcat(nonlinear,' + ');
    end
  end


end

end

