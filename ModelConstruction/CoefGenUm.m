function [coef] = CoefGenUm(m,p,q,cond)
%COEFGENUM Summary of this function goes here
%   Detailed explanation goes here

syms k1 V;
zeroCount = ZeroCounter(m(1))+ZeroCounter(m(2))+ZeroCounter(p(1))+ZeroCounter(p(2))+ZeroCounter(q(1))+ZeroCounter(q(2));
[S1,S2,S3,S4] = SignCoefficients(m,q,cond);
coef0 = (q(1)*p(2)*S1-p(1)*q(2)*S3);
prefactor = 2^(zeroCount/2-2);

coef = -k1*prefactor*coef0/(V*sqrt(p(1)^2*k1^2+p(2)^2));

end

