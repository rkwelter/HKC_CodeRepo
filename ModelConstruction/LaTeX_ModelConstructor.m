%% ModelConstructor
% ModelConstructor builds the time stepper for a desired HKC spectral model

M=2;
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
    title = strcat('HKC',num2str(M),'NoShear_LaTeX');
else
    title = strcat('HKC',num2str(M),'_LaTeX');
end

FID=fopen(strcat(title,'.m'),'w');
introString = [strcat('function [dX] = ',title,'(X,Pr,Ra,Ro,k1,V)'),newline,newline,'\begin{equation} X = \begin{pmatrix} '];
for i =1:numUp
  introString = [introString,strcat('u^+_{(',num2str(uPlusModes(i,1)),',',num2str(uPlusModes(i,2)),')} \\ '),newline];
end
for i =1:numUm
  introString = [introString,strcat('u^-_{(',num2str(uMinusModes(i,1)),',',num2str(uMinusModes(i,2)),')} \\ '),newline];
end
for i =1:numTheta
  introString = [introString,strcat('\theta_{(',num2str(thetaModes(i,1)),',',num2str(thetaModes(i,2)),')} \\ '),newline];
end
introString = [introString,'\end{pmatrix} \end{equation}',newline,newline,'\begin{equation} \begin{split} '];
fprintf(FID, '%s', introString);


%% uPlus terms

for i=1:numUp


    ODEstring = strcat('\frac{d}{dt} u^+_{(',num2str(uPlusModes(i,1)),',',num2str(uPlusModes(i,2)),')} & = ');

    %viscous term
    KmSq = AnisoLaplacian_LaTeX(uPlusModes(i,:));
    if(uPlusModes(i,1) > 0 && uPlusModes(i,2)>0)
      ODEstring = strcat(ODEstring,'-(',KmSq,')\Pra u^+_{(',num2str(uPlusModes(i,1)),',',num2str(uPlusModes(i,2)),')}');
    else
      ODEstring = strcat(ODEstring,'-',KmSq,'\Pra u^+_{(',num2str(uPlusModes(i,1)),',',num2str(uPlusModes(i,2)),')}');
    end

    %buoyancy term
    if(uPlusModes(i,1) > 0 && uPlusMatches(i,2)>0)
        if(mod(uPlusModes(i,1)+uPlusModes(i,2),2)==0)
            ODEstring = strcat(ODEstring,' - \frac{',num2str(uPlusModes(i,1)),'\Pra \Ray k_1}{\sqrt{',KmSq,'}} \theta_{(',num2str(uPlusModes(i,1)),',',num2str(uPlusModes(i,2)),')}');
        else
            ODEstring = strcat(ODEstring,' + \frac{',num2str(uPlusModes(i,1)),'\Pra \Ray k_1}{\sqrt{',KmSq,'}} \theta_{(',num2str(uPlusModes(i,1)),',',num2str(uPlusModes(i,2)),')}');
        end
    end

    %Coriolis term
    if(uPlusModes(i,2) > 0 && uPlusMatches(i,1)>0)
        if(uPlusModes(i,1)==0)
            ODEstring = strcat(ODEstring,' + \Pra \Rot u^-_{(',num2str(uPlusModes(i,1)),',',num2str(uPlusModes(i,2)),')}');
        else
            ODEstring = strcat(ODEstring,' + \frac{',num2str(uPlusModes(i,2)),'\Pra \Rot }{\sqrt{',KmSq,'}} u^-_{(',num2str(uPlusModes(i,1)),',',num2str(uPlusModes(i,2)),')}');
        end
    end

    %Nonlinear term
    triples = cell2mat(uPlusTriples(i));
    if(isempty(triples))
        ODEstring = [strcat(ODEstring,' \\ '),newline];
    else
        nonlinear = LaTeX_NonlinearUp(i,triples,uPlusModes);
        if(uPlusModes(i,1)==0)
          ODEstring = [strcat(ODEstring,' - \frac{k_1}{',num2str(4*uPlusModes(i,2)),'V} \Big [',nonlinear,' \Big ] \\'),newline];
        else
          ODEstring = [strcat(ODEstring,' - \frac{k_1}{4V \sqrt{',KmSq,'}} \Big [',nonlinear,' \Big ] \\'),newline];
        end
    end

    fprintf(FID, '%s', ODEstring);
end

%% uMinus terms

for i=1:numUm

    ODEstring = strcat('\frac{d}{dt} u^-_{(',num2str(uMinusModes(i,1)),',',num2str(uMinusModes(i,2)),')} & = ');

    %viscous term
    KmSq = AnisoLaplacian_LaTeX(uMinusModes(i,:));
    if(uMinusModes(i,1) > 0 && uMinusModes(i,2)>0)
      ODEstring = strcat(ODEstring,'-(',KmSq,')\Pra u^-_{(',num2str(uMinusModes(i,1)),',',num2str(uMinusModes(i,2)),')}');
    else
      ODEstring = strcat(ODEstring,'-',KmSq,'\Pra u^-_{(',num2str(uMinusModes(i,1)),',',num2str(uMinusModes(i,2)),')}');
    end


    %Coriolis term
    if(uMinusModes(i,2) > 0 && uMinusMatches(i,1)>0)
        if(uMinusModes(i,1)==0)
            ODEstring = strcat(ODEstring,' - \Pra \Rot u^+_{(',num2str(uMinusModes(i,1)),',',num2str(uMinusModes(i,2)),')}');
        else
            ODEstring = strcat(ODEstring,' - \frac{',num2str(uMinusModes(i,2)),'\Pra \Rot }{\sqrt{',KmSq,'}} u^-_{(',num2str(uMinusModes(i,1)),',',num2str(uMinusModes(i,2)),')}');
        end
    end

    %Nonlinear term
    triples = cell2mat(uMinusTriples(i));
    if(isempty(triples))
        ODEstring = [strcat(ODEstring,' \\ '),newline];
    else
        nonlinear = LaTeX_NonlinearUm(i,triples,uPlusModes,uMinusModes);
        ODEstring = [strcat(ODEstring,' - \frac{k_1}{4V} \Big [',nonlinear,' \Big ] \\ '),newline];
    end

    fprintf(FID, '%s', ODEstring);
end

%% theta terms

for i=1:numTheta

    ODEstring = strcat('\frac{d}{dt} \theta_{(',num2str(thetaModes(i,1)),',',num2str(thetaModes(i,2)),')} & = ');

    %viscous term
    KmSq = AnisoLaplacian_LaTeX(thetaModes(i,:));
    if(thetaModes(i,1) > 0 && thetaModes(i,2)>0)
      ODEstring = strcat(ODEstring,'-(',KmSq,') \theta_{(',num2str(thetaModes(i,1)),',',num2str(thetaModes(i,2)),')}');
    else
      ODEstring = strcat(ODEstring,'-',KmSq,' \theta_{(',num2str(thetaModes(i,1)),',',num2str(thetaModes(i,2)),')}');
    end

    %buoyancy term
    if(thetaModes(i,1) > 0 && thetaMatches(i,1)>0)
        if(mod(thetaModes(i,1)+thetaModes(i,2),2)==0)
            ODEstring = strcat(ODEstring,' - \frac{',num2str(thetaModes(i,1)),' k_1}{\sqrt{',KmSq,'}} u^+_{(',num2str(thetaModes(i,1)),',',num2str(thetaModes(i,2)),')}');
        else
            ODEstring = strcat(ODEstring,' + \frac{',num2str(thetaModes(i,1)),' k_1}{\sqrt{',KmSq,'}} u^+_{(',num2str(thetaModes(i,1)),',',num2str(thetaModes(i,2)),')}');
        end
    end

    %Nonlinear term
    triples = cell2mat(thetaTriples(i));
    if(isempty(triples))
    else
        nonlinear = LaTeX_NonlinearTheta(i,triples,uPlusModes,thetaModes,numUm);
        ODEstring = [strcat(ODEstring,' - \frac{k_1}{4V} \Big [',nonlinear,' \Big ] ')];
    end
    if(i < numTheta)
        ODEstring = [strcat(ODEstring,' \\ '),newline];
    else

    end

    fprintf(FID, '%s', ODEstring);
end

closeString = ['\end{split} \end{equation}',newline,'end'];
fprintf(FID, '%s', closeString);
fclose(FID);
