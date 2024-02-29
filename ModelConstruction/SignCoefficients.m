function [S1,S2,S3,S4] = SignCoefficients(m,q,indexCond)
%SignCoefficients by Roland Welter 
% Computes the sign coefficients in the nonlinear term

if(indexCond==1)
  S1 = (-1)^((m(1)+m(2)+1)*(q(1)+q(2)));
  S2 = S1;
  S3 = S1;
  S4 = S1;
elseif(indexCond==2)
  S1 = (-1)^((m(1)+m(2)+1)*(q(1)+q(2)));
  S2 = -S1;
  S3 = -S1;
  S4 = S1;
elseif(indexCond==3)
  S1 = (-1)^((m(1)+m(2)+1)*(q(1)+q(2)));
  S2 = S1;
  S3 = -S1;
  S4 = -S1;
elseif(indexCond==4)
  S1 = (-1)^((m(1)+m(2))*(q(1)+q(2)));
  S2 = -S1;
  S3 = -S1;
  S4 = S1;
elseif(indexCond==5)
  S1 = (-1)^((m(1)+m(2))*(q(1)+q(2)));
  S2 = S1;
  S3 = S1;
  S4 = S1;
elseif(indexCond==6)
  S1 = (-1)^((m(1)+m(2))*(q(1)+q(2)));
  S2 = -S1;
  S3 = S1;
  S4 = -S1;
elseif(indexCond==7)
  S1 = (-1)^((m(1)+m(2))*(q(1)+q(2)+1)+1);
  S2 = S1;
  S3 = -S1;
  S4 = -S1;
elseif(indexCond==8)
  S1 = (-1)^((m(1)+m(2))*(q(1)+q(2)+1)+1);
  S2 = -S1;
  S3 = S1;
  S4 = -S1;
elseif(indexCond==9)
  S1 = (-1)^((m(1)+m(2))*(q(1)+q(2)+1)+1);
  S2 = S1;
  S3 = S1;
  S4 = S1;
else
  error('Improper index condition');
end

end
