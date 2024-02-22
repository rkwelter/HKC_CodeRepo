function indexCond = CompatibilityTest(m,p,q)
%COMPATIBILITYTEST Summary of this function goes here
%   Detailed explanation goes here

if(m(1)== p(1)+q(1) && m(2)==p(2)+q(2))
    indexCond = 1;
elseif(m(1)== p(1)+q(1) && m(2)==p(2)-q(2))
    indexCond = 2;
elseif(m(1)== p(1)+q(1) && m(2)==-p(2)+q(2))
    indexCond = 3;
elseif(m(1)== p(1)-q(1) && m(2)==p(2)+q(2))
    indexCond = 4;
elseif(m(1)== p(1)-q(1) && m(2)==p(2)-q(2))
    indexCond = 5;
elseif(m(1)== p(1)-q(1) && m(2)==-p(2)+q(2))
    indexCond = 6;
elseif(m(1)== -p(1)+q(1) && m(2)==p(2)+q(2))
    indexCond = 7;
elseif(m(1)== -p(1)+q(1) && m(2)==p(2)-q(2))
    indexCond = 8;
elseif(m(1)== -p(1)+q(1) && m(2)==-p(2)+q(2))
    indexCond = 9;
else
    indexCond = 0;
end

end

