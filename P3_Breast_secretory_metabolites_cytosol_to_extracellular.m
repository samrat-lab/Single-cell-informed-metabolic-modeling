clear
clc

%% Load the models of P3_Breast
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P3_Breast_tINIT.mat')

%%%%%%%%%%% Cluster 1 %%%%%%%%%%%

%%%% Secretory Metabolites 

model=P3_Breast_tINIT.clst1.model;
model_flux=P3_Breast_tINIT.clst1.final_flux_values;

% Threshold to avoid numerical noise
fluxThreshold = 1e-6;

P3_Breast_clst1_secretory_rxns={};
n=0;
for i = 1:length(model.rxns)
    mets = model.mets(model.S(:,i) ~= 0);
    coeffs = full(model.S(model.S(:,i) ~= 0,i));
    
    % Get compartments (last character of metabolite ID)
    comps = cellfun(@(x) x(end), mets, 'UniformOutput', false);
    
    if  any(contains(comps,'c')) && any(contains(comps,'e'))
        if  any(coeffs(contains(comps,'e'))>0) && model_flux(i)> fluxThreshold
            n=n+1;
            comps
            P3_Breast_clst1_secretory_rxns{n,1}=model.rxns{i};
            P3_Breast_clst1_secretory_rxns{n,2}=mets;
            P3_Breast_clst1_secretory_rxns{n,3}=coeffs;
            P3_Breast_clst1_secretory_rxns{n,4}=model_flux(i);  
        end 
    end  
end

P3_Breast_clst1_secretory_mets=strings();
for i=1:size(P3_Breast_clst1_secretory_rxns,1)
    tmp=P3_Breast_clst1_secretory_rxns{i,2};
    for j=1:length(tmp)
        P3_Breast_clst1_secretory_mets=[P3_Breast_clst1_secretory_mets; tmp{j}(1:end-1)];
    end
    
end
P3_Breast_clst1_secretory_mets=unique(P3_Breast_clst1_secretory_mets(2:end)); % First cell is empty

save P3_Breast_clst1_secretory_rxns P3_Breast_clst1_secretory_rxns
save P3_Breast_clst1_secretory_mets P3_Breast_clst1_secretory_mets

clearvars -except P3_Breast_tINIT model model_flux fluxThreshold P3_Breast_clst1_secretory_rxns P3_Breast_clst1_secretory_mets


%%%%%%%%%%% Cluster 2 %%%%%%%%%%%

%%%% Secretory Metabolites 

model=P3_Breast_tINIT.clst2.model;
model_flux=P3_Breast_tINIT.clst2.final_flux_values;

% Threshold to avoid numerical noise
fluxThreshold = 1e-6;

P3_Breast_clst2_secretory_rxns={};
n=0;
for i = 1:length(model.rxns)
    mets = model.mets(model.S(:,i) ~= 0);
    coeffs = full(model.S(model.S(:,i) ~= 0,i));
    
    % Get compartments (last character of metabolite ID)
    comps = cellfun(@(x) x(end), mets, 'UniformOutput', false);
    
    if  any(contains(comps,'c')) && any(contains(comps,'e'))
        if  any(coeffs(contains(comps,'e'))>0) && model_flux(i)> fluxThreshold
            n=n+1;
            comps
            P3_Breast_clst2_secretory_rxns{n,1}=model.rxns{i};
            P3_Breast_clst2_secretory_rxns{n,2}=mets;
            P3_Breast_clst2_secretory_rxns{n,3}=coeffs;
            P3_Breast_clst2_secretory_rxns{n,4}=model_flux(i);  
        end 
    end  
end

P3_Breast_clst2_secretory_mets=strings();
for i=1:size(P3_Breast_clst2_secretory_rxns,1)
    tmp=P3_Breast_clst2_secretory_rxns{i,2};
    for j=1:length(tmp)
        P3_Breast_clst2_secretory_mets=[P3_Breast_clst2_secretory_mets; tmp{j}(1:end-1)];
    end
    
end
P3_Breast_clst2_secretory_mets=unique(P3_Breast_clst2_secretory_mets(2:end)); % First cell is empty

save P3_Breast_clst2_secretory_rxns P3_Breast_clst2_secretory_rxns
save P3_Breast_clst2_secretory_mets P3_Breast_clst2_secretory_mets

clearvars -except P3_Breast_clst1_secretory_rxns P3_Breast_clst1_secretory_mets ...
    P3_Breast_clst2_secretory_rxns P3_Breast_clst2_secretory_mets


P3_Breast_secretory_mets=intersect(P3_Breast_clst1_secretory_mets, P3_Breast_clst2_secretory_mets); % 39 secretory metabolites

union(P3_Breast_clst1_secretory_mets, P3_Breast_clst2_secretory_mets); % 41 secretory metabolites

save P3_Breast_secretory_mets P3_Breast_secretory_mets

% P3_Breast_clst1_secretory_mets - 41 secretory metabolites
% P3_Breast_clst2_secretory_mets - 39 secretory metabolites
% Common between them - 39 secretory metabolites

