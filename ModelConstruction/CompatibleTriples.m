function [uPlusTriples,uMinusTriples,thetaTriples] = CompatibleTriples(uPlusModes,uMinusModes,thetaModes)
%COMPATIBLETRIPLES For each of the variables in a given HKC model, this
%function determines which quadratic terms belong in the variable's
%nonlinear term

uPlusTriples = cell(size(uPlusModes,1),1);
uMinusTriples = cell(size(uMinusModes,1),1);
thetaTriples = cell(size(thetaModes,1),1);

for i=1:size(uPlusModes,1)
    triples = [];
    for j=1:size(uPlusModes,1)
        for k=1:size(uPlusModes,1)
            indexCond = CompatibilityTest(uPlusModes(i,:),uPlusModes(j,:),uPlusModes(k,:));
            if(indexCond>0)
                triple = [j;k;indexCond];
                triples = [triples triple];
            end
        end
    end
    uPlusTriples(i,1) = {triples};
end

for i=1:size(uMinusModes,1)
    triples = [];
    for j=1:size(uPlusModes,1)
        for k=1:size(uMinusModes,1)
            indexCond = CompatibilityTest(uMinusModes(i,:),uPlusModes(j,:),uMinusModes(k,:));
            if(indexCond>0)
                triple = [j;k;indexCond];
                triples = [triples triple];
            end
        end
    end
    uMinusTriples(i,1) = {triples};
end

for i=1:size(thetaModes,1)
    triples = [];
    for j=1:size(uPlusModes,1)
        for k=1:size(thetaModes,1)
            indexCond = CompatibilityTest(thetaModes(i,:),uPlusModes(j,:),thetaModes(k,:));
            if(indexCond>0)
                triple = [j;k;indexCond];
                triples = [triples triple];
            end
        end
    end
    thetaTriples(i,1) = {triples};
end

end

