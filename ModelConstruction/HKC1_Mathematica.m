function [dX] =HKC1_Mathematica(X,Pr,Ra,Ro,k1,V)

X = { uP0v1 ,
uP1v1 ,
uM0v1 ,
uM1v1 ,
th0v2 ,
th1v1 ,
}

dX = { -Pr*(1)*X[[1]] + Pr*Ro*1*X[[3]]/Sqrt[1],
-Pr*(k1^2*1 +1)*X[[2]] - Pr*Ra*k1*1*X[[6]]/Sqrt[k1^2*1 +1] + Pr*Ro*1*X[[4]]/Sqrt[k1^2*1 +1],
-Pr*(1)*X[[3]] - Pr*Ro*1*X[[1]]/Sqrt[1],
-Pr*(k1^2*1 +1)*X[[4]] - Pr*Ro*1*X[[2]]/Sqrt[k1^2*1 +1],
-(4)*X[[5]] - k1*(1*Sqrt[2]*(-2)*X[[2]]*X[[6]]/(Sqrt[k1^2*1 +1]))/(4*V),
-(k1^2*1 +1)*X[[6]] - k1*1*X[[2]]/Sqrt[k1^2*1 +1] - k1*(1*Sqrt[2]*(2)*X[[2]]*X[[5]]/(Sqrt[k1^2*1 +1]))/(4*V)}
end