function [dX] =HKC1(X,Pr,Ra,Ro,k1,V)

dX = [ -Pr*(1)*X(1) + Pr*Ro*1*X(3)/sqrt(1);
-Pr*(k1^2*1 +1)*X(2) - Pr*Ra*k1*1*X(6)/sqrt(k1^2*1 +1) + Pr*Ro*1*X(4)/sqrt(k1^2*1 +1);
-Pr*(1)*X(3) - Pr*Ro*1*X(1)/sqrt(1);
-Pr*(k1^2*1 +1)*X(4) - Pr*Ro*1*X(2)/sqrt(k1^2*1 +1);
-(4)*X(5) - k1*(1*sqrt(2)*(-2)*X(2)*X(6)/(sqrt(k1^2*1 +1)))/(4*V);
-(k1^2*1 +1)*X(6) - k1*1*X(2)/sqrt(k1^2*1 +1) - k1*(1*sqrt(2)*(2)*X(2)*X(5)/(sqrt(k1^2*1 +1)))/(4*V)];
end