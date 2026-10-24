%% Finding the flux values of reactions corresponding to each secretory metabolites

clear
clc

% P1 Breast
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\P1_Breast_clst1_secretory_mets.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\P1_Breast_clst1_secretory_rxns.mat')

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\P1_Breast_clst2_secretory_mets.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\P1_Breast_clst2_secretory_rxns.mat')


% P2 Breast
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\P2_Breast_clst1_secretory_mets.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\P2_Breast_clst1_secretory_rxns.mat')

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\P2_Breast_clst2_secretory_mets.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\P2_Breast_clst2_secretory_rxns.mat')


% P3 Breast
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\P3_Breast_clst1_secretory_mets.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\P3_Breast_clst1_secretory_rxns.mat')

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\P3_Breast_clst2_secretory_mets.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\P3_Breast_clst2_secretory_rxns.mat')


% Load the metabolic perturbations
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Metastatic_level_changes_GSE131007\Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Metastatic_level_changes_GSE131007\Metastatic_level_pert_rxns_P3_Bone_final_GSE131007.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Metastatic_level_changes_GSE131007\Metastatic_level_pert_rxns_P3_Brain_final_GSE131007.mat')


% Total secretory metabolites
Total_secretory_mets=union(union(union(union(union(P1_Breast_clst1_secretory_mets, P1_Breast_clst2_secretory_mets), P2_Breast_clst1_secretory_mets), ...
    P2_Breast_clst2_secretory_mets), P3_Breast_clst1_secretory_mets), P3_Breast_clst2_secretory_mets);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P1_Breast_Clst1
temp1={};
for i=1:size(Total_secretory_mets,1)
    target=Total_secretory_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P1_Breast_clst1_secretory_rxns(:,2));
    if any(idx)
        ind=find(idx);
        % All reactions
        temp1{i,1}=sum(cell2mat(P1_Breast_clst1_secretory_rxns(ind,4)));
        g=P1_Breast_clst1_secretory_rxns(ind,:);
        % Perturbed reactions
        [~,~,p]=intersect(Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007(:,1),g(:,1));
        temp1{i,2}=sum(cell2mat(g(p,4)));
    end
    clear target idx ind g p
end
idx=find(cellfun(@isempty, temp1(:,1)));
temp1(idx,1)={0};
clear idx
idx=find(cellfun(@isempty, temp1(:,2)));
temp1(idx,2)={0};
clear idx
temp1=cell2mat(temp1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P1_Breast_Clst2
temp2={};
for i=1:size(Total_secretory_mets,1)
    target=Total_secretory_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P1_Breast_clst2_secretory_rxns(:,2));
    if any(idx)
        ind=find(idx);
        % All reactions
        temp2{i,1}=sum(cell2mat(P1_Breast_clst2_secretory_rxns(ind,4)));
        g=P1_Breast_clst2_secretory_rxns(ind,:);
        % Perturbed reactions
        [~,~,p]=intersect(Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007(:,1),g(:,1));
        temp2{i,2}=sum(cell2mat(g(p,4)));
    end
    clear target idx ind g p
end
idx=find(cellfun(@isempty, temp2(:,1)));
temp2(idx,1)={0};
clear idx
idx=find(cellfun(@isempty, temp2(:,2)));
temp2(idx,2)={0};
clear idx
temp2=cell2mat(temp2);


% P2_Breast_Clst1
temp3={};
for i=1:size(Total_secretory_mets,1)
    target=Total_secretory_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P2_Breast_clst1_secretory_rxns(:,2));
    if any(idx)
        ind=find(idx);
        % All reactions
        temp3{i,1}=sum(cell2mat(P2_Breast_clst1_secretory_rxns(ind,4)));
        g=P2_Breast_clst1_secretory_rxns(ind,:);
        % Perturbed reactions
        [~,~,p]=intersect(Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007(:,1),g(:,1));
        temp3{i,2}=sum(cell2mat(g(p,4)));
    end
    clear target idx ind g p
end
idx=find(cellfun(@isempty, temp3(:,1)));
temp3(idx,1)={0};
clear idx
idx=find(cellfun(@isempty, temp3(:,2)));
temp3(idx,2)={0};
clear idx
temp3=cell2mat(temp3);


% P2_Breast_Clst2
temp4={};
for i=1:size(Total_secretory_mets,1)
    target=Total_secretory_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P2_Breast_clst2_secretory_rxns(:,2));
    if any(idx)
        ind=find(idx);
        % All reactions
        temp4{i,1}=sum(cell2mat(P2_Breast_clst2_secretory_rxns(ind,4)));
        g=P2_Breast_clst2_secretory_rxns(ind,:);
        % Perturbed reactions
        [~,~,p]=intersect(Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007(:,1),g(:,1));
        temp4{i,2}=sum(cell2mat(g(p,4)));
    end
    clear target idx ind g p
end
idx=find(cellfun(@isempty, temp4(:,1)));
temp4(idx,1)={0};
clear idx
idx=find(cellfun(@isempty, temp4(:,2)));
temp4(idx,2)={0};
clear idx
temp4=cell2mat(temp4);


%%%%% Bone
% P3_Breast_Clst1
temp5={};
for i=1:size(Total_secretory_mets,1)
    target=Total_secretory_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P3_Breast_clst1_secretory_rxns(:,2));
    if any(idx)
        ind=find(idx);
        % All reactions
        temp5{i,1}=sum(cell2mat(P3_Breast_clst1_secretory_rxns(ind,4)));
        g=P3_Breast_clst1_secretory_rxns(ind,:);
        % Perturbed reactions
        [~,~,p]=intersect(Metastatic_level_pert_rxns_P3_Bone_final_GSE131007(:,1),g(:,1));
        temp5{i,2}=sum(cell2mat(g(p,4)));
    end
    clear target idx ind g p
end
idx=find(cellfun(@isempty, temp5(:,1)));
temp5(idx,1)={0};
clear idx
idx=find(cellfun(@isempty, temp5(:,2)));
temp5(idx,2)={0};
clear idx
temp5=cell2mat(temp5);


% P3_Breast_Clst2
temp6={};
for i=1:size(Total_secretory_mets,1)
    target=Total_secretory_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P3_Breast_clst2_secretory_rxns(:,2));
    if any(idx)
        ind=find(idx);
        % All reactions
        temp6{i,1}=sum(cell2mat(P3_Breast_clst2_secretory_rxns(ind,4)));
        g=P3_Breast_clst2_secretory_rxns(ind,:);
        % Perturbed reactions
        [~,~,p]=intersect(Metastatic_level_pert_rxns_P3_Bone_final_GSE131007(:,1),g(:,1));
        temp6{i,2}=sum(cell2mat(g(p,4)));
    end
    clear target idx ind g p
end
idx=find(cellfun(@isempty, temp6(:,1)));
temp6(idx,1)={0};
clear idx
idx=find(cellfun(@isempty, temp6(:,2)));
temp6(idx,2)={0};
clear idx
temp6=cell2mat(temp6);


%%%%% Brain
% P3_Breast_Clst1
temp7={};
for i=1:size(Total_secretory_mets,1)
    target=Total_secretory_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P3_Breast_clst1_secretory_rxns(:,2));
    if any(idx)
        ind=find(idx);
        % All reactions
        temp7{i,1}=sum(cell2mat(P3_Breast_clst1_secretory_rxns(ind,4)));
        g=P3_Breast_clst1_secretory_rxns(ind,:);
        % Perturbed reactiona
        [~,~,p]=intersect(Metastatic_level_pert_rxns_P3_Brain_final_GSE131007(:,1),g(:,1));
        temp7{i,2}=sum(cell2mat(g(p,4)));
    end
    clear target idx ind g p
end
idx=find(cellfun(@isempty, temp7(:,1)));
temp7(idx,1)={0};
clear idx
idx=find(cellfun(@isempty, temp7(:,2)));
temp7(idx,2)={0};
clear idx
temp7=cell2mat(temp7);


% P3_Breast_Clst2
temp8={};
for i=1:size(Total_secretory_mets,1)
    target=Total_secretory_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P3_Breast_clst2_secretory_rxns(:,2));
    if any(idx)
        ind=find(idx);
        % All reactions
        temp8{i,1}=sum(cell2mat(P3_Breast_clst2_secretory_rxns(ind,4)));
        g=P3_Breast_clst2_secretory_rxns(ind,:);
        % Perturbed reactions
        [~,~,p]=intersect(Metastatic_level_pert_rxns_P3_Brain_final_GSE131007(:,1),g(:,1));
        temp8{i,2}=sum(cell2mat(g(p,4)));
    end
    clear target idx ind g p
end
idx=find(cellfun(@isempty, temp8(:,1)));
temp8(idx,1)={0};
clear idx
idx=find(cellfun(@isempty, temp8(:,2)));
temp8(idx,2)={0};
clear idx
temp8=cell2mat(temp8);


%% Combining all of them (All reactions)
g=[sum([temp1(:,1),temp2(:,1),temp3(:,1),temp4(:,1)],2)./sum([temp1(:,1),temp2(:,1),temp3(:,1),temp4(:,1)]~=0,2),...
    sum([temp5(:,1),temp6(:,1)],2)./sum([temp5(:,1),temp6(:,1)]~=0,2),...
    sum([temp7(:,1),temp8(:,1)],2)./sum([temp7(:,1),temp8(:,1)]~=0,2)];
g(isnan(g))=0; % To avoid NaN values

h=log10(g+1);


%% Combining all of them (Perturbed reactions)
g1=[sum([temp1(:,2),temp2(:,2),temp3(:,2),temp4(:,2)],2)./sum([temp1(:,2),temp2(:,2),temp3(:,2),temp4(:,2)]~=0,2),...
    sum([temp5(:,2),temp6(:,2)],2)./sum([temp5(:,2),temp6(:,2)]~=0,2),...
    sum([temp7(:,2),temp8(:,2)],2)./sum([temp7(:,2),temp8(:,2)]~=0,2)];
g1(isnan(g1))=0; % To avoid NaN values

h1=log10(g1+1);


% %% Combining all of them
% g=[sum([temp1(:),temp2(:)],2)./sum([temp1(:),temp2(:)]~=0,2),...
%     sum([temp3(:),temp4(:)],2)./sum([temp3(:),temp4(:)]~=0,2),...
%     sum([temp5(:),temp6(:)],2)./sum([temp5(:),temp6(:)]~=0,2)];
% g(isnan(g))=0; % To avoid NaN values
% 
% h=log10(g+1);


% Finding the metabolite names of secretory metabolites
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Human-GEM-main\model\Human-GEM.mat')

ihuman_met_ID= string(cellfun(@(x) x(1:end-1), ihuman.mets, 'UniformOutput', false));

[~,~,p]=intersect(Total_secretory_mets(:,1),ihuman_met_ID(:,1));

Total_secretory_met_names=ihuman.metNames(p,1);

Total_secretory_mets_flux_values_all=[["Total_secretory_met_IDs","Total_secretory_met_names","P1_Breast","P2_Breast","P3_Breast"];...
    Total_secretory_mets,Total_secretory_met_names,h];


Total_secretory_mets_flux_values_perturbed=[["Total_secretory_met_IDs","Total_secretory_met_names","P1_Breast","P2_Breast","P3_Breast"];...
    Total_secretory_mets,Total_secretory_met_names,h1];


save Total_secretory_mets_flux_values_all Total_secretory_mets_flux_values_all
save Total_secretory_mets_flux_values_perturbed Total_secretory_mets_flux_values_perturbed







