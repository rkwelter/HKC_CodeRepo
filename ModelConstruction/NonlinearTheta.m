function [nonlinear] = NonlinearTheta(i,triples,uPlusModes,thetaModes,numUm)
%NONLINEARTHETA by Roland Welter
%   writes the nonlinear terms for the theta variables

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
  if(zeroCount==0)
    prefactor = '';
  elseif(mod(zeroCount,2)==0)
    prefactor = strcat(num2str(2^(zeroCount/2)),'*');
  else
    prefactor = strcat(num2str(2^((zeroCount-1)/2)),'*sqrt(2)*');
  end

  if(mod(q(1)+q(2)+m(1)+m(2),2)==1)
    prefactor = strcat(prefactor,'(-1)*');
  end

  %Coefficient * Nonlinearity
  if(coef ~= 0)
    KpSq = AnisoLaplacian(p);
    nonlinear = strcat(nonlinear,prefactor,'(',num2str(coef),')*X(',num2str(triples(1,i)),')*X(',num2str(size(uPlusModes,1)+numUm+triples(2,i)),')');

    if( p(1) >0 || p(2) > 1 )
        nonlinear = strcat(nonlinear,'/sqrt(',KpSq,')');
    else
    end

    %End factor
    if(i < numTrips)
      nonlinear = strcat(nonlinear,' + '); %if the last terms are zero then this might end with a plus +
    end

  else % uncomment for verbose output in order to check for zero coefficients
      %strcat('Zero coef for theta, m = (',num2str(p(1)),',',num2str(p(2)),') , p = (',num2str(q(1)),',',num2str(q(2)),') , q = (',num2str(m(1)),',',num2str(m(2)),')')
  end

end

end

