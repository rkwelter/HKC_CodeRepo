function [nonlinear] = Mathematica_NonlinearUp(i,triples,uPlusModes)
%NONLINEARUP Summary of this function goes here
%   Detailed explanation goes here

nonlinear = '';
numTrips = size(triples,2);
m = uPlusModes(i,:);

for i=1:numTrips

  p = uPlusModes(triples(1,i),:);
  q = uPlusModes(triples(2,i),:);
  indexCond = triples(3,i);

  zeroCount = ZeroCounter(m(1))+ZeroCounter(m(2))+ZeroCounter(p(1))+ZeroCounter(p(2))+ZeroCounter(q(1))+ZeroCounter(q(2));
  [S1,S2,S3,S4] = SignCoefficients(m,q,indexCond);
  coef0 = (q(1)*p(2)*S1-p(1)*q(2)*S3)*q(2)*m(2);
  coef1 = (q(1)*p(2)*S2-p(1)*q(2)*S4)*q(1)*m(1);

  %Prefactor
  if(zeroCount==0)
    prefactor = '';
  elseif(mod(zeroCount,2)==0)
    prefactor = strcat(num2str(2^(zeroCount/2)),'*');
  else
    prefactor = strcat(num2str(2^((zeroCount-1)/2)),'*Sqrt[2]*');
  end

  %Coefficient
  if(coef0 ~= 0)
    if(coef1 > 0)
      nonlinear = strcat(nonlinear,prefactor,'(',num2str(coef0),'+',num2str(coef1),'*k1^2)*');
    elseif(coef1<0)
      nonlinear = strcat(nonlinear,prefactor,'(',num2str(coef0),'-',num2str(-1*coef1),'*k1^2)*');
    else
      nonlinear = strcat(nonlinear,prefactor,'(',num2str(coef0),')*');
    end
  else
    if(coef1 ~= 0)
      nonlinear = strcat(nonlinear,prefactor,'(',num2str(coef1),'*k1^2)*');
    end
  end

  %Nonlinearity
  if(coef0 ~= 0 || coef1 ~= 0)
    KpSq = AnisoLaplacian(p);
    KqSq = AnisoLaplacian(q);
    
    nonlinear = strcat(nonlinear,'X[[',num2str(triples(1,i)),']]*X[[',num2str(triples(2,i)),']]/(Sqrt[');

    if( p(1) >0 || p(2) > 1 )
        nonlinear = strcat(nonlinear,'(',KpSq,')');
        if(q(1) > 0 || q(2) > 1)
            nonlinear = strcat(nonlinear,'*(',KqSq,')])');
        else
            nonlinear = strcat(nonlinear,'])');
        end
    else
        if(q(1) > 0 || q(2) > 1)
            nonlinear = strcat(nonlinear,'(',KqSq,')])');
        else
            nonlinear = strcat(nonlinear,'])');
        end
    end


    %End factor
    if(i < numTrips)
      nonlinear = strcat(nonlinear,' + ');
    end
  end


end

end

