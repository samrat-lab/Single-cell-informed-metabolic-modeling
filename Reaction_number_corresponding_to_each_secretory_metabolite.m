%% Finding the number of reactions corresponding to each secretory metabolites

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


%%%%%% P1_Breast
% P1_Breast_Clst1
temp1={};
temp11={};
for i=1:size(Total_secretory_mets,1)
    target=Total_secretory_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P1_Breast_clst1_secretory_rxns(:,2));
    if any(idx)
        ind=find(idx);
        temp1{i,1}=P1_Breast_clst1_secretory_rxns(ind,1);
        temp1{i,2}=length(P1_Breast_clst1_secretory_rxns(ind,1));
        g=P1_Breast_clst1_secretory_rxns(ind,1);
        [~,~,p]=intersect(Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007(:,1),g);
        temp11{i,1}=g(p,1);
        temp11{i,2}=length(g(p,1));
    end
    clear target idx ind g p
end
idx=find(cellfun(@isempty, temp1(:,2)));
temp1(idx,2)={0};
clear idx
idx=find(cellfun(@isempty, temp11(:,2)));
temp11(idx,2)={0};
clear idx


% P1_Breast_Clst2
temp2={};
temp22={};
for i=1:size(Total_secretory_mets,1)
    target=Total_secretory_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P1_Breast_clst2_secretory_rxns(:,2));
    if any(idx)
        ind=find(idx);
        temp2{i,1}=P1_Breast_clst2_secretory_rxns(ind,1);
        temp2{i,2}=length(P1_Breast_clst2_secretory_rxns(ind,1));
        g=P1_Breast_clst2_secretory_rxns(ind,1);
        [~,~,p]=intersect(Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007(:,1),g);
        temp22{i,1}=g(p,1);
        temp22{i,2}=length(g(p,1));
    end
    clear target idx ind g p
end
idx=find(cellfun(@isempty, temp2(:,2)));
temp2(idx,2)={0};
clear idx
idx=find(cellfun(@isempty, temp22(:,2)));
temp22(idx,2)={0};
clear idx


%%%%%% P2_Breast
% P2_Breast_Clst1
temp3={};
temp33={};
for i=1:size(Total_secretory_mets,1)
    target=Total_secretory_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P2_Breast_clst1_secretory_rxns(:,2));
    if any(idx)
        ind=find(idx);
        temp3{i,1}=P2_Breast_clst1_secretory_rxns(ind,1);
        temp3{i,2}=length(P2_Breast_clst1_secretory_rxns(ind,1));
        g=P2_Breast_clst1_secretory_rxns(ind,1);
        [~,~,p]=intersect(Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007(:,1),g);
        temp33{i,1}=g(p,1);
        temp33{i,2}=length(g(p,1));
    end
    clear target idx ind g p
end
idx=find(cellfun(@isempty, temp3(:,2)));
temp3(idx,2)={0};
clear idx
idx=find(cellfun(@isempty, temp33(:,2)));
temp33(idx,2)={0};
clear idx


% P2_Breast_Clst2
temp4={};
temp44={};
for i=1:size(Total_secretory_mets,1)
    target=Total_secretory_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P2_Breast_clst2_secretory_rxns(:,2));
    if any(idx)
        ind=find(idx);
        temp4{i,1}=P2_Breast_clst2_secretory_rxns(ind,1);
        temp4{i,2}=length(P2_Breast_clst2_secretory_rxns(ind,1));
        g=P2_Breast_clst2_secretory_rxns(ind,1);
        [~,~,p]=intersect(Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007(:,1),g);
        temp44{i,1}=g(p,1);
        temp44{i,2}=length(g(p,1));
    end
    clear target idx ind g p
end
idx=find(cellfun(@isempty, temp4(:,2)));
temp4(idx,2)={0};
clear idx
idx=find(cellfun(@isempty, temp44(:,2)));
temp44(idx,2)={0};
clear idx

%%%%% PERTURBED REACTIONS
%%%%% Combining P1_Breast_Clst1, P1_Breast_Clst2, P2_Breast_Clst1 and P2_Breast_Clst2 (Since both are liver)
U={};
for i = 1:size(temp11,1)
    U{i,1} = union(union(string(temp11{i,1}), string(temp22{i,1})), ...
        union(string(temp33{i,1}), string(temp44{i,1})));
    U{i,2} = length(U{i,1});
end


%%%%% ALL REACTIONS
%%%%% Combining P1_Breast_Clst1, P1_Breast_Clst2, P2_Breast_Clst1 and P2_Breast_Clst2 (Since both are liver)
V={};
for i = 1:size(temp1,1)
    V{i,1} = union(union(string(temp1{i,1}), string(temp2{i,1})), ...
        union(string(temp3{i,1}), string(temp4{i,1})));
    V{i,2} = length(V{i,1});
end


%%%%%% P3_Breast (With Bone)
% P3_Breast_Clst1
temp5={};
temp55={};
for i=1:size(Total_secretory_mets,1)
    target=Total_secretory_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P3_Breast_clst1_secretory_rxns(:,2));
    if any(idx)
        ind=find(idx);
        temp5{i,1}=P3_Breast_clst1_secretory_rxns(ind,1);
        temp5{i,2}=length(P3_Breast_clst1_secretory_rxns(ind,1));
        g=P3_Breast_clst1_secretory_rxns(ind,1);
        [~,~,p]=intersect(Metastatic_level_pert_rxns_P3_Bone_final_GSE131007(:,1),g);
        temp55{i,1}=g(p,1);
        temp55{i,2}=length(g(p,1));
    end
    clear target idx ind g p
end
idx=find(cellfun(@isempty, temp5(:,2)));
temp5(idx,2)={0};
clear idx
idx=find(cellfun(@isempty, temp55(:,2)));
temp55(idx,2)={0};
clear idx


% P3_Breast_Clst2
temp6={};
temp66={};
for i=1:size(Total_secretory_mets,1)
    target=Total_secretory_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P3_Breast_clst2_secretory_rxns(:,2));
    if any(idx)
        ind=find(idx);
        temp6{i,1}=P3_Breast_clst2_secretory_rxns(ind,1);
        temp6{i,2}=length(P3_Breast_clst2_secretory_rxns(ind,1));
        g=P3_Breast_clst2_secretory_rxns(ind,1);
        [~,~,p]=intersect(Metastatic_level_pert_rxns_P3_Bone_final_GSE131007(:,1),g);
        temp66{i,1}=g(p,1);
        temp66{i,2}=length(g(p,1));
    end
    clear target idx ind g p
end
idx=find(cellfun(@isempty, temp6(:,2)));
temp6(idx,2)={0};
clear idx
idx=find(cellfun(@isempty, temp66(:,2)));
temp66(idx,2)={0};
clear idx


%%%%% PERTURBED REACTIONS
%%%%% Combining P3_Breast_Clst1 with Bone, P3_Breast_Clst2 with Bone
U1={};
for i = 1:size(temp55,1)
    U1{i,1} = union(string(temp55{i,1}), string(temp66{i,1}));
    U1{i,2} = length(U1{i,1});
end


%%%%% ALL REACTIONS
%%%%% Combining P3_Breast_Clst1 with Bone, P3_Breast_Clst2 with Bone
V1={};
for i = 1:size(temp5,1)
    V1{i,1} = union(string(temp5{i,1}), string(temp6{i,1}));
    V1{i,2} = length(V1{i,1});
end

%%%%%% P3_Breast (With Brain)
% P3_Breast_Clst1
temp7={};
temp77={};
for i=1:size(Total_secretory_mets,1)
    target=Total_secretory_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P3_Breast_clst1_secretory_rxns(:,2));
    if any(idx)
        ind=find(idx);
        temp7{i,1}=P3_Breast_clst1_secretory_rxns(ind,1);        
        temp7{i,2}=length(P3_Breast_clst1_secretory_rxns(ind,1));
        g=P3_Breast_clst1_secretory_rxns(ind,1);
        [~,~,p]=intersect(Metastatic_level_pert_rxns_P3_Brain_final_GSE131007(:,1),g);
        temp77{i,1}=g(p,1);
        temp77{i,2}=length(g(p,1));
    end
    clear target idx ind g p
end
idx=find(cellfun(@isempty, temp7(:,2)));
temp7(idx,2)={0};
clear idx
idx=find(cellfun(@isempty, temp77(:,2)));
temp77(idx,2)={0};
clear idx


% P3_Breast_Clst2
temp8={};
temp88={};
for i=1:size(Total_secretory_mets,1)
    target=Total_secretory_mets(i,1);
    idx=cellfun(@(x) any(contains(x,target)),P3_Breast_clst2_secretory_rxns(:,2));
    if any(idx)
        ind=find(idx);
        temp8{i,1}=P3_Breast_clst2_secretory_rxns(ind,1);        
        temp8{i,2}=length(P3_Breast_clst2_secretory_rxns(ind,1));
        g=P3_Breast_clst2_secretory_rxns(ind,1);
        [~,~,p]=intersect(Metastatic_level_pert_rxns_P3_Brain_final_GSE131007(:,1),g);
        temp88{i,1}=g(p,1);
        temp88{i,2}=length(g(p,1));
    end
    clear target idx ind g p
end
idx=find(cellfun(@isempty, temp8(:,2)));
temp8(idx,2)={0};
clear idx
idx=find(cellfun(@isempty, temp88(:,2)));
temp88(idx,2)={0};
clear idx


%%%%% PERTURBED REACTIONS
%%%%% Combining P3_Breast_Clst1 with Brain, P3_Breast_Clst2 with Brain
U2={};
for i = 1:size(temp77,1)
    U2{i,1} = union(string(temp77{i,1}), string(temp88{i,1}));
    U2{i,2} = length(U2{i,1});
end


%%%%% ALL REACTIONS
%%%%% Combining P3_Breast_Clst1 with Brain, P3_Breast_Clst2 with Brain
V2={};
for i = 1:size(temp7,1)
    V2{i,1} = union(string(temp7{i,1}), string(temp8{i,1}));
    V2{i,2} = length(V2{i,1});
end


% %% Combining all of them
% g=[max([temp1(:),temp2(:)],[],2), max([temp3(:),temp4(:)],[],2), max([temp5(:),temp6(:)],[],2)];


%% Combining all of them
g=[U(:,2), U1(:,2), U2(:,2)];

g1=[V(:,2), V1(:,2), V2(:,2)];

% Finding the metabolite names of secretory metabolites
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Human-GEM-main\model\Human-GEM.mat')

ihuman_met_ID= string(cellfun(@(x) x(1:end-1), ihuman.mets, 'UniformOutput', false));

[~,~,p]=intersect(Total_secretory_mets(:,1),ihuman_met_ID(:,1));

Total_secretory_met_names=ihuman.metNames(p,1);

Total_secretory_mets_no_of_rxns_perturbed=[["Total_secretory_met_IDs","Total_secretory_met_names","Liver","Bone","Brain"];...
    [Total_secretory_mets,Total_secretory_met_names,g]];

Total_secretory_mets_no_of_rxns_all=[["Total_secretory_met_IDs","Total_secretory_met_names","Liver","Bone","Brain"];...
    [Total_secretory_mets,Total_secretory_met_names,g1]];

save Total_secretory_mets_no_of_rxns_perturbed Total_secretory_mets_no_of_rxns_perturbed
save Total_secretory_mets_no_of_rxns_all Total_secretory_mets_no_of_rxns_all


