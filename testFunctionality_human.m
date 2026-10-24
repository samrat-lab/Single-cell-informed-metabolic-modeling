function [totTest, numTests, mets] = testFunctionality_human(model)
[~,bmMets,~] = xlsread('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\biomass_rxn_human.xlsx');
%     bmMets = setdiff(bmMets,{'atp[c]','adp[c]','pi[c]','h2o[c]','h[c]'});
bmMets = setdiff(bmMets,{'MAM01371c','MAM01285c','MAM02751c','MAM02040c','MAM02039c'});

numTests = numel(bmMets) + 1;
% comp = {'[c]','[e]','[g]','[l]','[m]','[n]','[r]','[x]'};%,'[i]'
% comp = {'c','e','g','l','m','n','r','x','i'};
comp = {'e','x','m','c','l','r','g','n','i'};


totTest = 0;
for i = 1:numel(bmMets)
    
%     cMet = bmMets{i}(1:end-3);
    cMet = bmMets{i}(1:end-1); % since we have compartments written as single alphabet at the end rather than in [].

    cTest = 0;
    for j = 1:numel(comp)
        met = [cMet,comp{j}];
%         met
       
        if ~isempty(intersect(model.mets,met))
            
            nm = addReaction(model,['BMS_',met],{met},-0.5,false,0,1000,0,'Biomass Metabolite Sink');
            nm = changeObjective(nm, ['BMS_',met]);
            sol = optimizeCbModel(nm);
            
            if sol.f > 1e-1
                cTest = 1;
                mets{i}=met;
                
                break;
            end
        end
    end
    if isempty(intersect(model.mets,met))
            a(i)=0;
            continue
    else
            totTest = cTest + totTest;

       a(i)=sol.f; 
    end
    
end

warning off all
%     nm = addReaction(model,'BMS_atp(c)',{'atp[c]','h2o[c]','adp[c]','pi[c]','h[c]'},[-0.5,-0.5,0.5,0.5,0.5],false,0,1000,0,'Biomass Metabolite Sink'); %Make a copy of ATP demand
nm = addReaction(model,'BMS_MAM01371(c)',{'MAM01371c','MAM02040c','MAM01285c','MAM02751c','MAM02039c'},[-0.5,-0.5,0.5,0.5,0.5],false,0,1000,0,'Biomass Metabolite Sink'); %Make a copy of ATP demand
%     nm = changeRxnBounds(nm, 'DM_atp_c_',0,'l');
nm = changeRxnBounds(nm, 'DM_MAM01371_c_',0,'l');
warning on all
% nm = changeObjective(nm, 'BMS_atp(c)');
nm = changeObjective(nm, 'BMS_MAM01371(c)');

sol = optimizeCbModel(nm);
if sol.f > 1e-1;
    totTest = totTest + 1;
end
a(i+1)=sol.f
disp(['Functionality test: ',num2str(totTest),'/',num2str(numel(bmMets)+1)])
end
