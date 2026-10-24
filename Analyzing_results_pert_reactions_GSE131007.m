clear
clc

% Load the perturbed reactions
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Metastatic_level_changes_GSE131007\Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Metastatic_level_changes_GSE131007\Metastatic_level_pert_rxns_P3_Bone_final_GSE131007.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Metastatic_level_changes_GSE131007\Metastatic_level_pert_rxns_P3_Brain_final_GSE131007.mat')


%% Analysing the common perturbed reactions (common to liver, bone, brain metastasis)
[p1,~,~]=intersect(intersect(Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007(2:end,1),Metastatic_level_pert_rxns_P3_Bone_final_GSE131007(2:end,1)),...
    Metastatic_level_pert_rxns_P3_Brain_final_GSE131007(2:end,1));

[p2,~,q1]=intersect(p1,Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007(2:end,1));
[~,~,q2]=intersect(p1,Metastatic_level_pert_rxns_P3_Bone_final_GSE131007(2:end,1));
[~,~,q3]=intersect(p1,Metastatic_level_pert_rxns_P3_Brain_final_GSE131007(2:end,1));

com_pert_rxns=[p2,Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007(q1+1,3),...
    Metastatic_level_pert_rxns_P3_Bone_final_GSE131007(q2+1,3),...
    Metastatic_level_pert_rxns_P3_Brain_final_GSE131007(q3+1,3)]; % 308 reactions are common to all tropism

clear p1 p2 q1 q2 q3

n=1;
temp=zeros(size(com_pert_rxns,1),1);
for i=1:size(com_pert_rxns)
    if find(com_pert_rxns(i,2)==com_pert_rxns(i,3) & com_pert_rxns(i,3)==com_pert_rxns(i,4))
    temp(n,1)=1;
    end
    n=n+1;
end

%%%%%% Common perturbed reactions with similar reaction status
com_pert_rxns=[com_pert_rxns,temp];
find(double(com_pert_rxns(:,5))==0); 
% Out of 308 common perturbed reactions 207 reactions have different directions in other tropisms, 
% resulting in 101 reactions with similar altered status.

ans1=find(double(com_pert_rxns(:,5))==1); 
com_pert_rxns_similar_status=[["Perturbed_Reactions", "Metastatic_P1_P2_Liver_reaction_status",...
    "Metastatic_P3_Bone_reaction_status", "Metastatic_P3_Brain_reaction_status", "Similar_score"];...
    com_pert_rxns(ans1,:)]; % 101 reactions

find(com_pert_rxns_similar_status(2:end,2)=='Same_Up'); % 81 reactions
find(com_pert_rxns_similar_status(2:end,2)=='Same_Down'); % 20 reactions
find(com_pert_rxns_similar_status(2:end,2)=='Opposite_Up'); % 0 reactions
find(com_pert_rxns_similar_status(2:end,2)=='Opposite_Down'); % 0 reactions

save com_pert_rxns_similar_status com_pert_rxns_similar_status

%%%%%% Common perturbed reactions with opposite reaction status
ans2=find(double(com_pert_rxns(:,5))==0);
com_pert_rxns_opposite_status=[["Perturbed_Reactions", "Metastatic_P1_P2_Liver_reaction_status",...
    "Metastatic_P3_Bone_reaction_status", "Metastatic_P3_Brain_reaction_status", "Similar_score"];...
    com_pert_rxns(ans2,:)]; % 207 reactions

save com_pert_rxns_opposite_status com_pert_rxns_opposite_status

clearvars -except com_pert_rxns com_pert_rxns_similar_status com_pert_rxns_opposite_status...
    Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007 Metastatic_level_pert_rxns_P3_Bone_final_GSE131007...
    Metastatic_level_pert_rxns_P3_Brain_final_GSE131007


%%%%%% Common perturbed reactions with opposite reaction status (same in P1_P2_Liver, P3_Bone but opposite in P3_Brain)
temp=com_pert_rxns_opposite_status(2:end,1:4);

n=1;
temp1=zeros(size(temp,1),1);
for i=1:size(temp)
    if find(temp(i,2)==temp(i,3) && temp(i,3)~=temp(i,4))
    temp1(n,1)=1;
    end
    n=n+1;
end

ans3=find(double(temp1(:,1))==1);
com_pert_rxns_opposite_status_Liver_Bone_vs_Brain=[["Perturbed_Reactions", "Metastatic_P1_P2_Liver_reaction_status",...
    "Metastatic_P3_Bone_reaction_status", "Metastatic_P3_Brain_reaction_status", "Similar_score"];...
    temp(ans3,:),temp1(ans3,1)]; % 39 reactions

clear temp n temp1

save com_pert_rxns_opposite_status_Liver_Bone_vs_Brain com_pert_rxns_opposite_status_Liver_Bone_vs_Brain

%%%%%% Common perturbed reactions with opposite reaction status (same in P1_P2_Liver, P3_Brain but opposite in P3_Bone)
temp=com_pert_rxns_opposite_status(2:end,1:4);

n=1;
temp1=zeros(size(temp,1),1);
for i=1:size(temp)
    if find(temp(i,2)==temp(i,4) && temp(i,4)~=temp(i,3))
    temp1(n,1)=1;
    end
    n=n+1;
end

ans3=find(double(temp1(:,1))==1);
com_pert_rxns_opposite_status_Liver_Brain_vs_Bone=[["Perturbed_Reactions", "Metastatic_P1_P2_Liver_reaction_status",...
    "Metastatic_P3_Bone_reaction_status", "Metastatic_P3_Brain_reaction_status", "Similar_score"];...
    temp(ans3,:),temp1(ans3,1)]; % 33 reactions

clear temp n temp1

save com_pert_rxns_opposite_status_Liver_Brain_vs_Bone com_pert_rxns_opposite_status_Liver_Brain_vs_Bone

%%%%%% Common perturbed reactions with opposite reaction status (same in P3_Bone, P3_Brain but opposite in P1_P2_Liver)
temp=com_pert_rxns_opposite_status(2:end,1:4);

n=1;
temp1=zeros(size(temp,1),1);
for i=1:size(temp)
    if find(temp(i,3)==temp(i,4) && temp(i,4)~=temp(i,2))
    temp1(n,1)=1;
    end
    n=n+1;
end

ans3=find(double(temp1(:,1))==1);
com_pert_rxns_opposite_status_Bone_Brain_vs_Liver=[["Perturbed_Reactions", "Metastatic_P1_P2_Liver_reaction_status",...
    "Metastatic_P3_Bone_reaction_status", "Metastatic_P3_Brain_reaction_status", "Similar_score"];...
    temp(ans3,:),temp1(ans3,1)]; % 113 reactions

clear temp n temp1

save com_pert_rxns_opposite_status_Bone_Brain_vs_Liver com_pert_rxns_opposite_status_Bone_Brain_vs_Liver


%%%%% Matrix for Sankey Plot
%%% Liver_Bone_vs_Brain
react_status=string(com_pert_rxns_opposite_status_Liver_Bone_vs_Brain(2:end,3:4)); % Since column 2 and 3 have same status
reactions=string(com_pert_rxns_opposite_status_Liver_Bone_vs_Brain(2:end,1));

labels=["Same_Up", "Same_Down", "Opposite_Up", "Opposite_Down"];

results=cell(size(labels,2),size(labels,2));
for i=1:length(labels)
    for j=1:length(labels)       
        idx_i=find(react_status(:,1)==labels(i));
        idx_j=find(react_status(:,2)==labels(j));
        
        common_idx=intersect(idx_i,idx_j);
        results{i,j}=reactions(common_idx);
    end
    clear idx_i idx_j common_idx
end

results=[{'Same_Up','Same_Down','Opposite_Up','Opposite_Down'}; results];
results=[[{'Labels'},{'Same_Up','Same_Down','Opposite_Up','Opposite_Down'}]',results];

empty_idx=cellfun(@isempty, results);
results(empty_idx)={0};

clear i j react_status reactions labels results idx_i idx_j common_idx empty_idx


%%% Liver_Brain_vs_Bone
react_status=string(com_pert_rxns_opposite_status_Liver_Brain_vs_Bone(2:end,3:4)); % Since column 2 and 4 have same status
reactions=string(com_pert_rxns_opposite_status_Liver_Brain_vs_Bone(2:end,1));

labels=["Same_Up", "Same_Down", "Opposite_Up", "Opposite_Down"];

results=cell(size(labels,2),size(labels,2));
for i=1:length(labels)
    for j=1:length(labels)       
        idx_i=find(react_status(:,1)==labels(i));
        idx_j=find(react_status(:,2)==labels(j));
        
        common_idx=intersect(idx_i,idx_j);
        results{i,j}=reactions(common_idx);
    end
    clear idx_i idx_j common_idx
end

results=[{'Same_Up','Same_Down','Opposite_Up','Opposite_Down'}; results];
results=[[{'Labels'},{'Same_Up','Same_Down','Opposite_Up','Opposite_Down'}]',results];

empty_idx=cellfun(@isempty, results);
results(empty_idx)={0};

clear i j react_status reactions labels results idx_i idx_j common_idx empty_idx


%%% Bone_Brain_vs_Liver
react_status=string(com_pert_rxns_opposite_status_Bone_Brain_vs_Liver(2:end,2:4)); % Since column 3 and 4 have same status
reactions=string(com_pert_rxns_opposite_status_Bone_Brain_vs_Liver(2:end,1));

labels=["Same_Up", "Same_Down", "Opposite_Up", "Opposite_Down"];

results=cell(size(labels,2),size(labels,2));
for i=1:length(labels)
    for j=1:length(labels)       
        idx_i=find(react_status(:,1)==labels(i));
        idx_j=find(react_status(:,2)==labels(j));
        
        common_idx=intersect(idx_i,idx_j);
        results{i,j}=reactions(common_idx);
    end
    clear idx_i idx_j common_idx
end

results=[{'Same_Up','Same_Down','Opposite_Up','Opposite_Down'}; results];
results=[[{'Labels'},{'Same_Up','Same_Down','Opposite_Up','Opposite_Down'}]',results];

empty_idx=cellfun(@isempty, results);
results(empty_idx)={0};

clear i j react_status reactions labels results idx_i idx_j common_idx empty_idx


%% Analysing the specific perturbed reactions

% Specific to Metatstatic_P1_P2_Liver
temp=Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007(2:end,:);
temp1=union(Metastatic_level_pert_rxns_P3_Bone_final_GSE131007(2:end,1),...
Metastatic_level_pert_rxns_P3_Brain_final_GSE131007(2:end,1));

[a,b]=setdiff(temp(:,1),temp1(:,1));
Specific_Metastatic_P1_P2_Liver_pert_rxns=[a,temp(b,2:end)]; % 63 reactions

find(Specific_Metastatic_P1_P2_Liver_pert_rxns(:,3)=='Same_Up'); % 26 reactions
find(Specific_Metastatic_P1_P2_Liver_pert_rxns(:,3)=='Same_Down'); % 33 reactions
find(Specific_Metastatic_P1_P2_Liver_pert_rxns(:,3)=='Opposite_Up'); % 1 reaction
find(Specific_Metastatic_P1_P2_Liver_pert_rxns(:,3)=='Opposite_Down'); % 3 reactions

clear temp temp1 a b 

% Specific to Metastatic_P3_Bone
temp=Metastatic_level_pert_rxns_P3_Bone_final_GSE131007(2:end,:);
temp1=union(Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007(2:end,1),...
Metastatic_level_pert_rxns_P3_Brain_final_GSE131007(2:end,1));

[a,b]=setdiff(temp(:,1),temp1(:,1));
Specific_Metastatic_P3_Bone_pert_rxns=[a,temp(b,2:end)]; % 24 reactions

find(Specific_Metastatic_P3_Bone_pert_rxns(:,3)=='Same_Up'); % 3 reactions
find(Specific_Metastatic_P3_Bone_pert_rxns(:,3)=='Same_Down'); % 19 reactions
find(Specific_Metastatic_P3_Bone_pert_rxns(:,3)=='Opposite_Up'); % 0 reactions
find(Specific_Metastatic_P3_Bone_pert_rxns(:,3)=='Opposite_Down'); % 2 reactions

clear temp temp1 a b

% Specific to Metastatic_P3_Brain
temp=Metastatic_level_pert_rxns_P3_Brain_final_GSE131007(2:end,:);
temp1=union(Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007(2:end,1),...
Metastatic_level_pert_rxns_P3_Bone_final_GSE131007(2:end,1));

[a,b]=setdiff(temp(:,1),temp1(:,1));
Specific_Metastatic_P3_Brain_pert_rxns=[a,temp(b,2:end)]; % 19 reactions

find(Specific_Metastatic_P3_Brain_pert_rxns(:,3)=='Same_Up'); % 7 reactions
find(Specific_Metastatic_P3_Brain_pert_rxns(:,3)=='Same_Down'); % 8 reactions
find(Specific_Metastatic_P3_Brain_pert_rxns(:,3)=='Opposite_Up'); % 2 reactions
find(Specific_Metastatic_P3_Brain_pert_rxns(:,3)=='Opposite_Down'); % 2 reactions

clear temp temp1 a b 

save Specific_Metastatic_P1_P2_Liver_pert_rxns Specific_Metastatic_P1_P2_Liver_pert_rxns
save Specific_Metastatic_P3_Bone_pert_rxns Specific_Metastatic_P3_Bone_pert_rxns
save Specific_Metastatic_P3_Brain_pert_rxns Specific_Metastatic_P3_Brain_pert_rxns

