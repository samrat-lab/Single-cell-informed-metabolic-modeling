%%%%% Using it after the file: "Flux_values_evaluation_gp_sampler_GSE131007"
clear
clc

%%%% Evaluating parsedGPR and corrRxn
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Human-GEM-main\model\Human-GEM.mat')
% model=ihuman;
% [parsedGPR_ihuman,corrRxn_ihuman]=extractGPRs(model);

% save parsedGPR_ihuman parsedGPR_ihuman
% save corrRxn_ihuman corrRxn_ihuman


% Preparing a gene_exp_file for each model
% P1_UCD46
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\P1_Liver_clst1.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\P1_Liver_clst2.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\P1_Liver_met_genes.mat')

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\P1_Breast_clst1.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\P1_Breast_clst2.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\P1_Breast_met_genes.mat')

% P2_UCD4
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\P2_Liver_clst1.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\P2_Liver_clst2.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\P2_Liver_met_genes.mat')

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\P2_Breast_clst1.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\P2_Breast_clst2.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\P2_Breast_met_genes.mat')

% P3_UCD65
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\P3_Bone_clst1.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\P3_Bone_clst2.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\P3_Bone_met_genes.mat')

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\P3_Brain_clst1.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\P3_Brain_clst2.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\P3_Brain_met_genes.mat')

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\P3_Breast_clst1.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\P3_Breast_clst2.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\P3_Breast_met_genes.mat')

[~,p1,q1]=intersect(upper(string(ihuman.geneShortNames)), upper(P1_Liver_met_genes(:,1)));
[~,p2,q2]=intersect(upper(string(ihuman.geneShortNames)), upper(P1_Breast_met_genes(:,1)));
[~,p3,q3]=intersect(upper(string(ihuman.geneShortNames)), upper(P2_Liver_met_genes(:,1)));
[~,p4,q4]=intersect(upper(string(ihuman.geneShortNames)), upper(P2_Breast_met_genes(:,1)));
[~,p5,q5]=intersect(upper(string(ihuman.geneShortNames)), upper(P3_Bone_met_genes(:,1)));
[~,p6,q6]=intersect(upper(string(ihuman.geneShortNames)), upper(P3_Brain_met_genes(:,1)));
[~,p7,q7]=intersect(upper(string(ihuman.geneShortNames)), upper(P3_Breast_met_genes(:,1)));


str1={'P1_Liver_clst1', 'P1_Liver_clst2', 'P1_Breast_clst1', 'P1_Breast_clst2', ...
    'P2_Liver_clst1', 'P2_Liver_clst2', 'P2_Breast_clst1', 'P2_Breast_clst2', ...
    'P3_Bone_clst1', 'P3_Bone_clst2', 'P3_Brain_clst1', 'P3_Brain_clst2', 'P3_Breast_clst1', 'P3_Breast_clst2'};


met_gene_exp_union=zeros(size(ihuman.genes,1),size(str1,2));
for i=1:size(str1,2)
    met_gene_exp_union(p1,1)=P1_Liver_clst1(q1,1);
    met_gene_exp_union(p1,2)=P1_Liver_clst2(q1,1);
    
    met_gene_exp_union(p2,3)=P1_Breast_clst1(q2,1);
    met_gene_exp_union(p2,4)=P1_Breast_clst2(q2,1);
    
    met_gene_exp_union(p3,5)=P2_Liver_clst1(q3,1);
    met_gene_exp_union(p3,6)=P2_Liver_clst2(q3,1);
    
    met_gene_exp_union(p4,7)=P2_Breast_clst1(q4,1);
    met_gene_exp_union(p4,8)=P2_Breast_clst2(q4,1);
    
    met_gene_exp_union(p5,9)=P3_Bone_clst1(q5,1);
    met_gene_exp_union(p5,10)=P3_Bone_clst2(q5,1);
    
    met_gene_exp_union(p6,11)=P3_Brain_clst1(q6,1);
    met_gene_exp_union(p6,12)=P3_Brain_clst2(q6,1);
    
    met_gene_exp_union(p7,13)=P3_Breast_clst1(q7,1);
    met_gene_exp_union(p7,14)=P3_Breast_clst2(q7,1);
end

% save met_gene_exp_union met_gene_exp_union


%%%% Evaluating reaction expression array (rxn_exp_array) for each model
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P1_Breast_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P1_Liver_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P2_Breast_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P2_Liver_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P3_Bone_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P3_Brain_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P3_Breast_tINIT.mat')


% Loading parsedGPR and corrRxn of ihuman model
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\parsedGPR_ihuman.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\corrRxn_ihuman.mat')

category={'P1_Liver_tINIT.clst1', 'P1_Liver_tINIT.clst2', 'P1_Breast_tINIT.clst1', 'P1_Breast_tINIT.clst2', ...
    'P2_Liver_tINIT.clst1', 'P2_Liver_tINIT.clst2', 'P2_Breast_tINIT.clst1', 'P2_Breast_tINIT.clst2', ...
    'P3_Bone_tINIT.clst1', 'P3_Bone_tINIT.clst2', 'P3_Brain_tINIT.clst1', 'P3_Brain_tINIT.clst2', ...
    'P3_Breast_tINIT.clst1', 'P3_Breast_tINIT.clst2'};

for i=1:size(category,2)
    tmp=eval(category{i});
    model1=tmp.model;
    
    [p,q,r]=intersect(string(ihuman.genes),string(model1.genes));
    
    if isequal(string(model1.genes), string(model1.genes(r)))
        
        expValues=mapGeneToRxn(model1,model1.genes,met_gene_exp_union(q,i),parsedGPR_ihuman,corrRxn_ihuman);
        
        str=[category{i},'.rxn_exp_array', '=','expValues'];
        eval(str);
        
    end
    
    clear tmp model1 expValues str ans
    
end


save('P1_Liver_tINIT.mat', 'P1_Liver_tINIT', '-v7.3')
save('P1_Breast_tINIT.mat', 'P1_Breast_tINIT', '-v7.3')
save('P2_Liver_tINIT.mat', 'P2_Liver_tINIT', '-v7.3')
save('P2_Breast_tINIT.mat', 'P2_Breast_tINIT', '-v7.3')
save('P3_Bone_tINIT.mat', 'P3_Bone_tINIT', '-v7.3')
save('P3_Brain_tINIT.mat', 'P3_Brain_tINIT', '-v7.3')
save('P3_Breast_tINIT.mat', 'P3_Breast_tINIT', '-v7.3')

