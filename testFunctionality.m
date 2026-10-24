function [totTest, numTests, mets] = testFunctionality(model)
    [~,bmMets,~] = xlsread('biomass_rxn.xls');%'union.xls');
    bmMets = setdiff(bmMets,{'atp[c]','adp[c]','pi[c]','h2o[c]','h[c]'});
    numTests = numel(bmMets) + 1;
    comp = {'[c]','[e]','[g]','[l]','[m]','[n]','[r]','[x]'};
    totTest = 0;
    for i = 1:numel(bmMets)

        cMet = bmMets{i}(1:end-3);
        cTest = 0;
        for j = 1:numel(comp)
            met = [cMet,comp{j}];
            if ~isempty(intersect(model.mets,met))
                
                nm = addReaction(model,['BMS_',met],{met},-0.5,false,0,1000,0,'Biomass Metabolite Sink');
                nm = changeObjective(nm, ['BMS_',met]);
                sol = optimizeCbModel(nm);

                if sol.f > 1e-1;
                    cTest = 1;
                    mets{i}=met;

                    break;
                end
            end
        end
        totTest = cTest + totTest;
        a(i)=sol.f;
    end
    
    warning off all
    nm = addReaction(model,'BMS_atp(c)',{'atp[c]','h2o[c]','adp[c]','pi[c]','h[c]'},[-0.5,-0.5,0.5,0.5,0.5],false,0,1000,0,'Biomass Metabolite Sink'); %Make a copy of ATP demand
    nm = changeRxnBounds(nm, 'DM_atp_c_',0,'l');
    warning on all
    nm = changeObjective(nm, 'BMS_atp(c)');
    sol = optimizeCbModel(nm);
    if sol.f > 1e-1;
        totTest = totTest + 1;
    end
    a(i+1)=sol.f
    disp(['Functionality test: ',num2str(totTest),'/',num2str(numel(bmMets)+1)])
end
