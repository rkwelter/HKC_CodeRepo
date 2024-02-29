function [KmSq] = AnisoLaplacian(Modes)
%ANISOLAPLACIAN Summary of this function goes here
%   Detailed explanation goes here

    if(Modes(1) == 0 && Modes(2)==1)
        KmSq = '';
    elseif(Modes(1) == 0 && Modes(2)>1)
        KmSq = num2str(Modes(2)^2);
    elseif(Modes(2) == 0 && Modes(1)==1)
        KmSq = 'k1^2';
    elseif(Modes(2) == 0 && Modes(1)>1)
        KmSq = strcat(num2str(Modes(1)^2),'*k1^2');
    elseif(Modes(2) > 0 && Modes(1)==1)
        KmSq = strcat('k1^2 + ',num2str(Modes(2)^2));
    else
        KmSq = strcat(num2str(Modes(1)^2),'*k1^2 + ',num2str(Modes(2)^2));
    end

end

