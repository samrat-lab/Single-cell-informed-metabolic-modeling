%%%%% Using it after the file: "Flux_values_evaluation_gp_sampler_human"
clear
clc

%%%% Evaluating parsedGPR and corrRxn
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Human-GEM-main\model\Human-GEM.mat')
% model=ihuman;
% [parsedGPR_ihuman,corrRxn_ihuman]=extractGPRs(model);

% save parsedGPR_ihuman parsedGPR_ihuman
% save corrRxn_ihuman corrRxn_ihuman


%%%% Preparing a gene_exp_file for each model
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE235168_variables\Primary_breast_met_data_clst1.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE235168_variables\Primary_breast_met_data_clst2.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE235168_variables\Final_met_genes.mat')

[~,p,q]=intersect(string(ihuman.genes), Final_met_genes(:,1));

met_gene_exp_union=zeros(size(ihuman.genes,1),2);
for i=1:2
    met_gene_exp_union(p,1)=Primary_breast_met_data_clst1(1,q)';
    met_gene_exp_union(p,2)=Primary_breast_met_data_clst2(1,q)';
end


%%%% Evaluating reaction expression array (rxn_exp_array) for each model
% load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE235168_variables\Primary_breast_ftINIT.mat')

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE235168_variables\Primary_breast_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE235168_variables\parsedGPR_ihuman.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE235168_variables\corrRxn_ihuman.mat')

% category={'Primary_breast_ftINIT.clst1','Primary_breast_ftINIT.clst2'};
category={'Primary_breast_tINIT.clst1','Primary_breast_tINIT.clst2'};
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

% save('Primary_breast_ftINIT.mat', 'Primary_breast_ftINIT', '-v7.3')
save('Primary_breast_tINIT.mat', 'Primary_breast_tINIT', '-v7.3')

