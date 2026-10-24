%% Finding the flux values of reactions corresponding to each uptake metabolites

clear
clc

% P1 Breast
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Uptake_metabolites_and_reactions\P1_Breast_clst1_uptake_mets.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Uptake_metabolites_and_reactions\P1_Breast_clst1_uptake_rxns.mat')

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Uptake_metabolites_and_reactions\P1_Breast_clst2_uptake_mets.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Uptake_metabolites_and_reactions\P1_Breast_clst2_uptake_rxns.mat')


% P2 Breast
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Uptake_metabolites_and_reactions\P2_Breast_clst1_uptake_mets.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Uptake_metabolites_and_reactions\P2_Breast_clst1_uptake_rxns.mat')

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Uptake_metabolites_and_reactions\P2_Breast_clst2_uptake_mets.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Uptake_metabolites_and_reactions\P2_Breast_clst2_uptake_rxns.mat')


% P3 Breast
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Uptake_metabolites_and_reactions\P3_Breast_clst1_uptake_mets.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Uptake_metabolites_and_reactions\P3_Breast_clst1_uptake_rxns.mat')

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Uptake_metabolites_and_reactions\P3_Breast_clst2_uptake_mets.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Uptake_metabolites_and_reactions\P3_Breast_clst2_uptake_rxns.mat')


% Total uptake metabolites
Total_uptake_mets=union(union(union(union(union(P1_Breast_clst1_uptake_mets, P1_Breast_clst2_uptake_mets), P2_Breast_clst1_uptake_mets), ...
    P2_Breast_clst2_uptake_mets), P3_Breast_clst1_uptake_mets), P3_Breast_clst2_uptake_mets);


% P1_Breast_Clst1
temp1=zeros(size(Total_uptake_mets,1),1);
for i=1:size(Total_uptake_mets,1)
    target=Total_uptake_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P1_Breast_clst1_uptake_rxns(:,2));
    if any(idx)
    ind=find(idx);
    temp1(i,1)=length(P1_Breast_clst1_uptake_rxns(ind,1));
    end
    clear target idx ind
end

% P1_Breast_Clst2
temp2=zeros(size(Total_uptake_mets,1),1);
for i=1:size(Total_uptake_mets,1)
    target=Total_uptake_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P1_Breast_clst2_uptake_rxns(:,2));
    if any(idx)
    ind=find(idx);
    temp2(i,1)=length(P1_Breast_clst2_uptake_rxns(ind,1));
    end
    clear target idx ind
end


% P2_Breast_Clst1
temp3=zeros(size(Total_uptake_mets,1),1);
for i=1:size(Total_uptake_mets,1)
    target=Total_uptake_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P2_Breast_clst1_uptake_rxns(:,2));
    if any(idx)
    ind=find(idx);
    temp3(i,1)=length(P2_Breast_clst1_uptake_rxns(ind,1));
    end
    clear target idx ind
end

% P2_Breast_Clst2
temp4=zeros(size(Total_uptake_mets,1),1);
for i=1:size(Total_uptake_mets,1)
    target=Total_uptake_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P2_Breast_clst2_uptake_rxns(:,2));
    if any(idx)
    ind=find(idx);
    temp4(i,1)=length(P2_Breast_clst2_uptake_rxns(ind,1));
    end
    clear target idx ind
end


% P3_Breast_Clst1
temp5=zeros(size(Total_uptake_mets,1),1);
for i=1:size(Total_uptake_mets,1)
    target=Total_uptake_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P3_Breast_clst1_uptake_rxns(:,2));
    if any(idx)
    ind=find(idx);
    temp5(i,1)=length(P3_Breast_clst1_uptake_rxns(ind,1));
    end
    clear target idx ind
end


% P3_Breast_Clst2
temp6=zeros(size(Total_uptake_mets,1),1);
for i=1:size(Total_uptake_mets,1)
    target=Total_uptake_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P3_Breast_clst2_uptake_rxns(:,2));
    if any(idx)
    ind=find(idx);
    temp6(i,1)=length(P3_Breast_clst2_uptake_rxns(ind,1));
    end
    clear target idx ind
end


%% Combining all of them
g=[max([temp1(:),temp2(:)],[],2), max([temp3(:),temp4(:)],[],2), max([temp5(:),temp6(:)],[],2)];

% Finding the metabolite names of uptake metabolites
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Human-GEM-main\model\Human-GEM.mat')

ihuman_met_ID= string(cellfun(@(x) x(1:end-1), ihuman.mets, 'UniformOutput', false));

[~,~,p]=intersect(Total_uptake_mets(:,1),ihuman_met_ID(:,1));

Total_uptake_met_names=ihuman.metNames(p,1);

Total_uptake_mets_no_of_rxns=[["Total_uptake_met_IDs","Total_uptake_met_names","P1_Breast","P2_Breast","P3_Breast"];...
    Total_uptake_mets,Total_uptake_met_names,g];

save Total_uptake_mets_no_of_rxns Total_uptake_mets_no_of_rxns


