clear
clc

% Load the Human 1 model
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Human-GEM-main\model\Human-GEM.mat')

% Load the common pert reactions with opposite status
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Metastatic_level_changes_GSE131007\com_pert_rxns_opposite_status_Liver_Bone_vs_Brain.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Metastatic_level_changes_GSE131007\com_pert_rxns_opposite_status_Liver_Brain_vs_Bone.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Metastatic_level_changes_GSE131007\com_pert_rxns_opposite_status_Bone_Brain_vs_Liver.mat')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Liver Bone vs Brain %%%%%%%%%%%%%%%%%%%%%%%
%% Same Up in Liver and Bone vs Same Down in Brain
idx=find(com_pert_rxns_opposite_status_Liver_Bone_vs_Brain(2:end,2) == "Same_Up" & ...
    com_pert_rxns_opposite_status_Liver_Bone_vs_Brain(2:end,3) == "Same_Up" & ...
    com_pert_rxns_opposite_status_Liver_Bone_vs_Brain(2:end,4) == "Same_Down");

Same_up_to_Same_down_Liver_Bone_vs_Brain=com_pert_rxns_opposite_status_Liver_Bone_vs_Brain(idx+1,:);

[~,q,r]=intersect(Same_up_to_Same_down_Liver_Bone_vs_Brain(:,1),string(ihuman.rxns(:,1)));
Same_up_to_Same_down_Liver_Bone_vs_Brain=[Same_up_to_Same_down_Liver_Bone_vs_Brain(q,:),...
    string(ihuman.subSystems(r,1))];

% significant_pathways=FEA(ihuman,r,'subSystems');
clear idx q r

save Same_up_to_Same_down_Liver_Bone_vs_Brain Same_up_to_Same_down_Liver_Bone_vs_Brain

%% Same Down in Liver and Bone vs Same Up in Brain
idx=find(com_pert_rxns_opposite_status_Liver_Bone_vs_Brain(2:end,2) == "Same_Down" & ...
    com_pert_rxns_opposite_status_Liver_Bone_vs_Brain(2:end,3) == "Same_Down" & ...
    com_pert_rxns_opposite_status_Liver_Bone_vs_Brain(2:end,4) == "Same_Up");

Same_down_to_Same_up_Liver_Bone_vs_Brain=com_pert_rxns_opposite_status_Liver_Bone_vs_Brain(idx+1,:);

[~,q,r]=intersect(Same_down_to_Same_up_Liver_Bone_vs_Brain(:,1),string(ihuman.rxns(:,1)));
Same_down_to_Same_up_Liver_Bone_vs_Brain=[Same_down_to_Same_up_Liver_Bone_vs_Brain(q,:),...
    string(ihuman.subSystems(r,1))];

% significant_pathways=FEA(ihuman,r,'subSystems');
clear idx q r

save Same_down_to_Same_up_Liver_Bone_vs_Brain Same_down_to_Same_up_Liver_Bone_vs_Brain


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Liver Brain vs Bone %%%%%%%%%%%%%%%%%%%%%%%
%% Same Up in Liver and Brain vs Same Down in Bone
idx=find(com_pert_rxns_opposite_status_Liver_Brain_vs_Bone(2:end,2) == "Same_Up" & ...
    com_pert_rxns_opposite_status_Liver_Brain_vs_Bone(2:end,4) == "Same_Up" & ...
    com_pert_rxns_opposite_status_Liver_Brain_vs_Bone(2:end,3) == "Same_Down");

Same_up_to_Same_down_Liver_Brain_vs_Bone=com_pert_rxns_opposite_status_Liver_Brain_vs_Bone(idx+1,:);

[~,q,r]=intersect(Same_up_to_Same_down_Liver_Brain_vs_Bone(:,1),string(ihuman.rxns(:,1)));
Same_up_to_Same_down_Liver_Brain_vs_Bone=[Same_up_to_Same_down_Liver_Brain_vs_Bone(q,:),...
    string(ihuman.subSystems(r,1))];

% significant_pathways=FEA(ihuman,r,'subSystems');
clear idx q r

save Same_up_to_Same_down_Liver_Brain_vs_Bone Same_up_to_Same_down_Liver_Brain_vs_Bone


%% Same Down in Liver and Brain vs Same Up in Bone
idx=find(com_pert_rxns_opposite_status_Liver_Brain_vs_Bone(2:end,2) == "Same_Down" & ...
    com_pert_rxns_opposite_status_Liver_Brain_vs_Bone(2:end,4) == "Same_Down" & ...
    com_pert_rxns_opposite_status_Liver_Brain_vs_Bone(2:end,3) == "Same_Up");

Same_down_to_Same_up_Liver_Brain_vs_Bone=com_pert_rxns_opposite_status_Liver_Brain_vs_Bone(idx+1,:);

[~,q,r]=intersect(Same_down_to_Same_up_Liver_Brain_vs_Bone(:,1),string(ihuman.rxns(:,1)));
Same_down_to_Same_up_Liver_Brain_vs_Bone=[Same_down_to_Same_up_Liver_Brain_vs_Bone(q,:),...
    string(ihuman.subSystems(r,1))];

% significant_pathways=FEA(ihuman,r,'subSystems');
clear idx q r

save Same_down_to_Same_up_Liver_Brain_vs_Bone Same_down_to_Same_up_Liver_Brain_vs_Bone

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bone Brain vs Liver %%%%%%%%%%%%%%%%%%%%%%%
%% Same Up in Bone and Brain vs Same Down in Liver
idx=find(com_pert_rxns_opposite_status_Bone_Brain_vs_Liver(2:end,3) == "Same_Up" & ...
    com_pert_rxns_opposite_status_Bone_Brain_vs_Liver(2:end,4) == "Same_Up" & ...
    com_pert_rxns_opposite_status_Bone_Brain_vs_Liver(2:end,2) == "Same_Down");

Same_up_to_Same_down_Bone_Brain_vs_Liver=com_pert_rxns_opposite_status_Bone_Brain_vs_Liver(idx+1,:);

[~,q,r]=intersect(Same_up_to_Same_down_Bone_Brain_vs_Liver(:,1),string(ihuman.rxns(:,1)));
Same_up_to_Same_down_Bone_Brain_vs_Liver=[Same_up_to_Same_down_Bone_Brain_vs_Liver(q,:),...
    string(ihuman.subSystems(r,1))];

% significant_pathways=FEA(ihuman,r,'subSystems');
clear idx q r

save Same_up_to_Same_down_Bone_Brain_vs_Liver Same_up_to_Same_down_Bone_Brain_vs_Liver

%% Same Down in Bone and Brain vs Same Up in Liver
idx=find(com_pert_rxns_opposite_status_Bone_Brain_vs_Liver(2:end,3) == "Same_Down" & ...
    com_pert_rxns_opposite_status_Bone_Brain_vs_Liver(2:end,4) == "Same_Down" & ...
    com_pert_rxns_opposite_status_Bone_Brain_vs_Liver(2:end,2) == "Same_Up");

Same_down_to_Same_up_Bone_Brain_vs_Liver=com_pert_rxns_opposite_status_Bone_Brain_vs_Liver(idx+1,:);

[~,q,r]=intersect(Same_down_to_Same_up_Bone_Brain_vs_Liver(:,1),string(ihuman.rxns(:,1)));
Same_down_to_Same_up_Bone_Brain_vs_Liver=[Same_down_to_Same_up_Bone_Brain_vs_Liver(q,:),...
    string(ihuman.subSystems(r,1))];

% significant_pathways=FEA(ihuman,r,'subSystems');
clear idx q r

save Same_down_to_Same_up_Bone_Brain_vs_Liver Same_down_to_Same_up_Bone_Brain_vs_Liver

