function [coef] = CoefGenUp(m,p,q,cond)
%COEFGENUP Summary of this function goes here
%   Detailed explanation goes here

syms k1 V;
zeroCount = ZeroCounter(m(1))+ZeroCounter(m(2))+ZeroCounter(p(1))+ZeroCounter(p(2))+ZeroCounter(q(1))+ZeroCounter(q(2));
[S1,S2,S3,S4] = SignCoefficients(m,q,cond);
coef0 = (q(1)*p(2)*S1-p(1)*q(2)*S3)*q(2)*m(2);
coef1 = (q(1)*p(2)*S2-p(1)*q(2)*S4)*q(1)*m(1);
prefactor = 2^(zeroCount/2-2);

coef = -k1*prefactor*(coef0+coef1*k1^2)/(V*sqrt(p(1)^2*k1^2+p(2)^2)*sqrt(q(1)^2*k1^2+q(2)^2)*sqrt(m(1)^2*k1^2+m(2)^2));

end

