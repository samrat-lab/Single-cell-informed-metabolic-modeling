clear
clc

% % Initiating the Raven Toolbox
% initCobraToolbox(false)
% changeCobraSolver('gurobi', 'LP');  
checkInstallation
setRavenSolver('gurobi')

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE235168_variables\Human-GEM.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE235168_variables\taskStruct.mat')

% %% Preparing INIT Model
% origRefModel=ihuman;
% spontRxnNames=[];
% convertGenes=false;
% customRxnsToIgnore=[];
% extComp='e'; % specifies the label for the extracellular compartment in your metabolic model, ensuring that transport and exchange reactions are identified and processed correctly.
% skipScaling=false;
% prepData = prepINITModel(origRefModel, taskStruct);%, spontRxnNames, convertGenes, customRxnsToIgnore, extComp, skipScaling);
% 
% prepData_ihuman=prepData;

% save prepData_ihuman prepData_ihuman

% Preparing a gene_exp_file for building ftINIT Models
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE235168_variables\Primary_breast_met_data_clst1.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE235168_variables\Primary_breast_met_data_clst2.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE235168_variables\Final_met_genes.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE235168_variables\prepData_ihuman.mat')

[~,p,q]=intersect(string(ihuman.genes), Final_met_genes(:,1));

met_gene_exp_union=zeros(size(ihuman.genes,1),2);
for i=1:2
    met_gene_exp_union(p,1)=Primary_breast_met_data_clst1(1,q)';
    met_gene_exp_union(p,2)=Primary_breast_met_data_clst2(1,q)';
end
% save met_gene_exp_union met_gene_exp_union

threshold=[];
for i=1:size(met_gene_exp_union,2)
    x=met_gene_exp_union(:,i);
    threshold(:,i)=min(x(x>0));
    clear x
end

%% Building models using ftINIT
data_structure.genes=ihuman.genes;
tissues={'Primary_breast_met_data_clst1', 'Primary_breast_met_data_clst2'};
data_structure.tissues=tissues;
data_structure.levels=met_gene_exp_union;

for i=1:size(met_gene_exp_union,2)
    category={'clst1','clst2'};
    prepData=prepData_ihuman;
    celltype=[];
    hpaData=[];
    transcrData=data_structure;
    th=threshold(1,i);
    transcrData.threshold=th;
    metabolomicsData=[];
    INITSteps=getINITSteps([],'1+1');
    removeGenes=true;
    useScoresForTasks=true;
    paramsFT=[];
    verbose=true;
    
    tissue=transcrData.tissues{1,i};
    
    [model, metProduction, addedRxnsForTasks, deletedRxnsInINIT, fullMipRes] = ftINIT(prepData, tissue, celltype, ...
        hpaData, transcrData, metabolomicsData, INITSteps, removeGenes, useScoresForTasks, paramsFT, verbose);
    model1 = ftINIT(prepData, data_structure.tissues{1}, [], [], data_structure, {}, getHumanGEMINITSteps('1+0'), false, true);

    str=['Primary_breast_ftINIT','.',category{i},'.model','=','model'];
    
    eval(str);
    
end

save Primary_breast_ftINIT Primary_breast_ftINIT


%% Checking Metabolic Tasks
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\Primary_breast_ftINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\taskStruct.mat')

category={'Primary_breast_ftINIT.clst1', 'Primary_breast_ftINIT.clst2'};
for i=1:2 
    model=eval(category{i});
    [taskReport, ~, ~, ~]=checkTasks(model.model, [], true, false, false, taskStruct);
    
    results=taskReport.description(find(taskReport.ok));
    
    str=[category{i},'.task_report', '=','results'];
    eval(str);
    
    clear model taskReport str results
    
end


% %% Performing Flux variability analysis (FVA)
% initCobraToolbox(false)
% 
% category={'Primary_breast_ftINIT.clst1', 'Primary_breast_ftINIT.clst2'};
% for i=1:2
%     model=eval(category{i});
%     % osenseStr='max';  % --default
%     [minFlux, maxFlux] = fluxVariability(model.model);
%     
%     str=[category{i},'.FVA', '=','[minFlux, maxFlux]'];
%     eval(str);
%     
%     clear 'model' 'minFlux' 'maxFlux' 'str'
% end


%% Metabolic functionality
category={'Primary_breast_ftINIT.clst1', 'Primary_breast_ftINIT.clst2'};
for i=1:2
    model=eval(category{i});
    % osenseStr='max';  % --default
    [totTest, numTests] = testFunctionality_human(model.model);
    
    str=[category{i},'.met_fun', '=','[totTest, numTests]'];
    eval(str);
    
    clear model totTest numTests str
end

save('Primary_breast_ftINIT.mat', 'Primary_breast_ftINIT', '-v7.3')

