function indexCond = CompatibilityTest(m,p,q)
%COMPATIBILITYTEST tests whether wave numbers p,q are compatible with a
%wave number m

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

