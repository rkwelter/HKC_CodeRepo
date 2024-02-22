function [coef] = CoefGenTheta(m,p,q,cond)
%COEFGENTHETA Summary of this function goes here
%   Detailed explanation goes here

syms k1 V;
zeroCount = ZeroCounter(m(1))+ZeroCounter(m(2))+ZeroCounter(p(1))+ZeroCounter(p(2))+ZeroCounter(q(1))+ZeroCounter(q(2));
[S1,S2,S3,S4] = SignCoefficients(m,q,cond);
coef1 = (q(1)*p(2)*S2-p(1)*q(2)*S4);
prefactor = 2^(zeroCount/2-2);

coef = (-1)^(q(1)+q(2)+m(1)+m(2)+1)*k1*prefactor*coef1/(V*sqrt(p(1)^2*k1^2+p(2)^2));

end

