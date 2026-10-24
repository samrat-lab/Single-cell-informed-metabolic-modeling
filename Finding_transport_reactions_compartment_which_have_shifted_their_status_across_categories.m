%%%% Using this file after "Pathway_analysis_of_reactions_shifting_their_status_across_categories"

clear
clc

% Load the pathways list
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\Same_up_to_Same_down_Liver_Bone_vs_Brain.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\Same_down_to_Same_up_Liver_Bone_vs_Brain.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\Same_up_to_Same_down_Liver_Brain_vs_Bone.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\Same_down_to_Same_up_Liver_Brain_vs_Bone.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\Same_up_to_Same_down_Bone_Brain_vs_Liver.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\Same_down_to_Same_up_Bone_Brain_vs_Liver.mat')

% Load the Human1 model
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Human-GEM-main\model\Human-GEM.mat')

% Load the secretory reactions
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\P1_Breast_clst1_secretory_rxns.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\P1_Breast_clst2_secretory_rxns.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\P2_Breast_clst1_secretory_rxns.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\P2_Breast_clst2_secretory_rxns.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\P3_Breast_clst1_secretory_rxns.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\P3_Breast_clst2_secretory_rxns.mat')

Total_secretory_rxns=union(union(union(union(union(P1_Breast_clst1_secretory_rxns(:,1),P1_Breast_clst2_secretory_rxns(:,1)),...
    P2_Breast_clst1_secretory_rxns(:,1)),P2_Breast_clst2_secretory_rxns(:,1)),...
    P3_Breast_clst1_secretory_rxns(:,1)),P3_Breast_clst2_secretory_rxns(:,1));

category={'Same_up_to_Same_down_Liver_Bone_vs_Brain', 'Same_down_to_Same_up_Liver_Bone_vs_Brain',...
    'Same_up_to_Same_down_Liver_Brain_vs_Bone', 'Same_down_to_Same_up_Liver_Brain_vs_Bone',...
    'Same_up_to_Same_down_Bone_Brain_vs_Liver', 'Same_down_to_Same_up_Bone_Brain_vs_Liver'};

S=full(ihuman.S);

for i=1:size(category,2)
    str=category{1,i};
    str1=eval(str);
    
    find(str1(:,6)=='Transport reactions');
    rxns=cellstr(str1(ans,:));
    
    % Secretory reactions
    [~,~,q]=intersect(Total_secretory_rxns(:,1),rxns(:,1));
    temp=rxns(q,:);
    
    [~,p,~]=intersect(string(ihuman.rxns),temp(:,1));
    
    % Linked secretory metabolites
    for j=1:size(p,1)
    mets{j,1}=ihuman.mets(find(S(:,p(j,1))));
    end
    
    % Linked gene catalyzing secretory reactions
    for j=1:size(p,1)
    genes{j,1}=ihuman.geneShortNames(find(ihuman.rxnGeneMat(p(j,1),:)));
    end    
    
    % Assembling all results together
    temp=[temp,mets,genes];
    
    str2=['Transport_rxns.',category{1,i},'=temp'];
    eval(str2)
    
    clear str str1 rxns str2 ans mets temp genes
    
end


%% Nesting results and saving them in excel file
for i=1:size(category,2)
    str=['Transport_rxns.', category{1,i}];
    temp=eval(str);
    for j=1:size(temp,1)
        for k=1:size(temp,2)
            if iscell(temp{j,k})
                temp{j,k}=strjoin(string(temp{j,k}), ', ');
            end
        end
    end
    str1=['Transport_rxns_string.', category{1,i}, '=temp'];
    eval(str1);
    
    clear str temp str1
end





