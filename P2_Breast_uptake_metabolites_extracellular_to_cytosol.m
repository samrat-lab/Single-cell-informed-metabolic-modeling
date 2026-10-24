clear
clc

%% Load the models of P2_Breast
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P2_Breast_tINIT.mat')

%%%%%%%%%%% Cluster 1 %%%%%%%%%%%

%%%% Uptake Metabolites 

model=P2_Breast_tINIT.clst1.model;
model_flux=P2_Breast_tINIT.clst1.final_flux_values;

% Threshold to avoid numerical noise
fluxThreshold = 1e-6;

P2_Breast_clst1_uptake_rxns={};
n=0;
for i = 1:length(model.rxns)
    mets = model.mets(model.S(:,i) ~= 0);
    coeffs = full(model.S(model.S(:,i) ~= 0,i));
    
    % Get compartments (last character of metabolite ID)
    comps = cellfun(@(x) x(end), mets, 'UniformOutput', false);
    
    if  any(contains(comps,'c')) && any(contains(comps,'e'))
        if  any(coeffs(contains(comps,'c'))>0) && model_flux(i)> fluxThreshold
            n=n+1;
            comps
            P2_Breast_clst1_uptake_rxns{n,1}=model.rxns{i};
            P2_Breast_clst1_uptake_rxns{n,2}=mets;
            P2_Breast_clst1_uptake_rxns{n,3}=coeffs;
            P2_Breast_clst1_uptake_rxns{n,4}=model_flux(i);  
        end 
    end  
end

P2_Breast_clst1_uptake_mets=strings();
for i=1:size(P2_Breast_clst1_uptake_rxns,1)
    tmp=P2_Breast_clst1_uptake_rxns{i,2};
    for j=1:length(tmp)
        P2_Breast_clst1_uptake_mets=[P2_Breast_clst1_uptake_mets; tmp{j}(1:end-1)];
    end
    
end
P2_Breast_clst1_uptake_mets=unique(P2_Breast_clst1_uptake_mets(2:end)); % First cell is empty

save P2_Breast_clst1_uptake_rxns P2_Breast_clst1_uptake_rxns
save P2_Breast_clst1_uptake_mets P2_Breast_clst1_uptake_mets

clearvars -except P2_Breast_tINIT fluxThreshold P2_Breast_clst1_uptake_rxns P2_Breast_clst1_uptake_mets


%%%%%%%%%%% Cluster 2 %%%%%%%%%%%

%%%% Uptake Metabolites 

model=P2_Breast_tINIT.clst2.model;
model_flux=P2_Breast_tINIT.clst2.final_flux_values;

% Threshold to avoid numerical noise
fluxThreshold = 1e-6;

P2_Breast_clst2_uptake_rxns={};
n=0;
for i = 1:length(model.rxns)
    mets = model.mets(model.S(:,i) ~= 0);
    coeffs = full(model.S(model.S(:,i) ~= 0,i));
    
    % Get compartments (last character of metabolite ID)
    comps = cellfun(@(x) x(end), mets, 'UniformOutput', false);
    
    if  any(contains(comps,'c')) && any(contains(comps,'e'))
        if  any(coeffs(contains(comps,'c'))>0) && model_flux(i)> fluxThreshold
            n=n+1;
            comps
            P2_Breast_clst2_uptake_rxns{n,1}=model.rxns{i};
            P2_Breast_clst2_uptake_rxns{n,2}=mets;
            P2_Breast_clst2_uptake_rxns{n,3}=coeffs;
            P2_Breast_clst2_uptake_rxns{n,4}=model_flux(i);  
        end 
    end  
end

P2_Breast_clst2_uptake_mets=strings();
for i=1:size(P2_Breast_clst2_uptake_rxns,1)
    tmp=P2_Breast_clst2_uptake_rxns{i,2};
    for j=1:length(tmp)
        P2_Breast_clst2_uptake_mets=[P2_Breast_clst2_uptake_mets; tmp{j}(1:end-1)];
    end
    
end
P2_Breast_clst2_uptake_mets=unique(P2_Breast_clst2_uptake_mets(2:end)); % First cell is empty

save P2_Breast_clst2_uptake_rxns P2_Breast_clst2_uptake_rxns
save P2_Breast_clst2_uptake_mets P2_Breast_clst2_uptake_mets

clearvars -except P2_Breast_clst1_uptake_rxns P2_Breast_clst1_uptake_mets ...
    P2_Breast_clst2_uptake_rxns P2_Breast_clst2_uptake_mets

P2_Breast_uptake_mets=intersect(P2_Breast_clst1_uptake_mets, P2_Breast_clst2_uptake_mets); % 37 uptake metabolites

union(P2_Breast_clst1_uptake_mets, P2_Breast_clst2_uptake_mets); % 48 uptake metabolites

save P2_Breast_uptake_mets P2_Breast_uptake_mets

% P2_Breast_clst1_uptake_mets - 42 uptake metabolites
% P2_Breast_clst2_uptake_mets - 43 uptake metabolites
% Common between them - 37 uptake metabolites

