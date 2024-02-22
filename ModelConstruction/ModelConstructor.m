%% ModelConstructor
% ModelConstructor builds the time stepper for a desired HKC spectral model

M=45;
shear = 'yes';
[uPlusModes,uMinusModes,thetaModes] = VariableConstructor(M,shear);

numUp = size(uPlusModes,1);
numUm = size(uMinusModes,1);
numTheta = size(thetaModes,1);

%calculate compatible triples
[uPlusTriples,uMinusTriples,thetaTriples] = CompatibleTriples(uPlusModes,uMinusModes,thetaModes);
%calculate match indices
[uPlusMatches,uMinusMatches,thetaMatches] = ComputeMatches(uPlusModes,uMinusModes,thetaModes);

%% write function

if(strcmp(shear,'no'))
    title = strcat('HKC',num2str(M),'NoShear');
else
    title = strcat('HKC',num2str(M));
end

FID=fopen(strcat(title,'.m'),'w');
introString = [strcat('function [dX] = ',title,'(X,Pr,Ra,Ro,k1,V)'),newline,newline,'dX = [ '];
fprintf(FID, '%s', introString);

%% uPlus terms

for i=1:numUp

    %viscous term
    KmSq = AnisoLaplacian(uPlusModes(i,:));
    if(uPlusModes(i,1) == 0 && uPlusModes(i,2)==1)
        ODEstring = strcat('-Pr*X(',num2str(i),')');
    elseif((uPlusModes(i,1) >= 1 && uPlusModes(i,2)==0) || (uPlusModes(i,1) == 0 && uPlusModes(i,2)>1))
        ODEstring = strcat('-Pr*',KmSq,'*X(',num2str(i),')');
    else
        ODEstring = strcat('-Pr*(',KmSq,')*X(',num2str(i),')');
    end

    %buoyancy term
    if(uPlusModes(i,1) > 0 && uPlusMatches(i,2)>0)
        if(mod(uPlusModes(i,1)+uPlusModes(i,2),2)==0)
            ODEstring = strcat(ODEstring,' - Pr*Ra*k1*',num2str(uPlusModes(i,1)),'*X(',num2str(numUp+numUm+uPlusMatches(i,2)),')/sqrt(',KmSq,')');
        else
            ODEstring = strcat(ODEstring,' + Pr*Ra*k1*',num2str(uPlusModes(i,1)),'*X(',num2str(numUp+numUm+uPlusMatches(i,2)),')/sqrt(',KmSq,')');
        end
    end

    %Coriolis term
    if(uPlusModes(i,2) > 0 && uPlusMatches(i,1)>0)
        if(uPlusModes(i,1) == 0)
            ODEstring = strcat(ODEstring,' + Pr*Ro*X(',num2str(numUp+uPlusMatches(i,1)),')');
        else
            ODEstring = strcat(ODEstring,' + Pr*Ro*',num2str(uPlusModes(i,2)),'*X(',num2str(numUp+uPlusMatches(i,1)),')/sqrt(',KmSq,')');
        end
    end

    %Nonlinear term
    triples = cell2mat(uPlusTriples(i));
    if(isempty(triples))
        ODEstring = [strcat(ODEstring,';'),newline];
    else
        nonlinear = NonlinearUp(i,triples,uPlusModes);
        if( uPlusModes(i,1) > 0 || uPlusModes(i,2) > 1)
            ODEstring = [strcat(ODEstring,' - k1*(',nonlinear,')/(4*V*sqrt(',KmSq,'));'),newline];
        else
            ODEstring = [strcat(ODEstring,' - k1*(',nonlinear,')/(4*V);'),newline];
        end
    end

    fprintf(FID, '%s', ODEstring);
end

%% uMinus terms

for i=1:numUm

    %viscous term
    KmSq = AnisoLaplacian(uMinusModes(i,:));
    if(uMinusModes(i,1) == 0 && uMinusModes(i,2)==1)
        ODEstring = strcat('-Pr*X(',num2str(numUp+i),')');
    elseif((uMinusModes(i,1) >= 1 && uMinusModes(i,2)==0) || (uMinusModes(i,1) == 0 && uMinusModes(i,2)>1))
        ODEstring = strcat('-Pr*',KmSq,'*X(',num2str(numUp+i),')');
    else
        ODEstring = strcat('-Pr*(',KmSq,')*X(',num2str(numUp+i),')');
    end

    %Coriolis term
    if(uMinusModes(i,2) > 0 && uMinusMatches(i,1)>0)
        if(uMinusModes(i,1) == 0)
            ODEstring = strcat(ODEstring,' - Pr*Ro*X(',num2str(uMinusMatches(i,1)),')');
        else
            ODEstring = strcat(ODEstring,' - Pr*Ro*',num2str(uMinusModes(i,2)),'*X(',num2str(uMinusMatches(i,1)),')/sqrt(',KmSq,')');
        end
    end

    %Nonlinear term
    triples = cell2mat(uMinusTriples(i));
    if(isempty(triples))
        ODEstring = [strcat(ODEstring,';'),newline];
    else
        nonlinear = NonlinearUm(i,triples,uPlusModes,uMinusModes);
        ODEstring = [strcat(ODEstring,' - k1*(',nonlinear,')/(4*V);'),newline];
    end

    fprintf(FID, '%s', ODEstring);
end

%% theta terms

for i=1:numTheta

    %viscous term
    KmSq = AnisoLaplacian(thetaModes(i,:));
    if(thetaModes(i,1) == 0 && thetaModes(i,2)==1)
        ODEstring = strcat('-X(',num2str(numUp+numUm+i),')');
    elseif((thetaModes(i,1) >= 1 && thetaModes(i,2)==0) || (thetaModes(i,1) == 0 && thetaModes(i,2)>1))
        ODEstring = strcat('-',KmSq,'*X(',num2str(numUp+numUm+i),')');
    else
        ODEstring = strcat('-(',KmSq,')*X(',num2str(numUp+numUm+i),')');
    end

    %buoyancy term
    if(thetaModes(i,1) > 0 && thetaMatches(i,1)>0)
        if(mod(thetaModes(i,1)+thetaModes(i,2),2)==0)
            ODEstring = strcat(ODEstring,' - k1*',num2str(thetaModes(i,1)),'*X(',num2str(thetaMatches(i,1)),')/sqrt(',KmSq,')');
        else
            ODEstring = strcat(ODEstring,' + k1*',num2str(thetaModes(i,1)),'*X(',num2str(thetaMatches(i,1)),')/sqrt(',KmSq,')');
        end
    end

    %Nonlinear term
    triples = cell2mat(thetaTriples(i));
    if(isempty(triples))
    else
        nonlinear = NonlinearTheta(i,triples,uPlusModes,thetaModes,numUm);
        ODEstring = [strcat(ODEstring,' - k1*(',nonlinear,')/(4*V)')];
    end
    if(i < numTheta)
        ODEstring = [strcat(ODEstring,';'),newline];
    else

    end

    fprintf(FID, '%s', ODEstring);
end

closeString = ['];',newline,'end'];
fprintf(FID, '%s', closeString);
fclose(FID);
