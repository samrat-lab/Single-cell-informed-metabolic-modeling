clear
clc

%%%% Initiating the Raven Toolbox
% initCobraToolbox(false)
% changeCobraSolver('gurobi','LP', 1);
% changeCobraSolver('gurobi','MILP', 1);

checkInstallation
setRavenSolver('gurobi')

% load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\Human-GEM-main\model\Human-GEM.mat')

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\Human-GEM-main\model\Human-GEM.mat')

% load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\taskStruct.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\prepDataHumanGEM.mat')

taskStruct=prepDataHumanGEM.taskStruct;

% %% Preparing INIT Model
% origRefModel=ihuman;
% spontRxnNames=[];
% convertGenes=false;
% customRxnsToIgnore=[];
% extComp='e'; % specifies the label for the extracellular compartment in your metabolic model, ensuring that transport and exchange reactions are identified and processed correctly.
% skipScaling=false;
% prepData = prepINITModel(origRefModel, taskStruct);%, spontRxnNames, convertGenes, customRxnsToIgnore, extComp, skipScaling);

% prepData_ihuman=prepData;
% % save prepData_ihuman prepData_ihuman

% %% Load the prepData_ihuman
% load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\prepData_ihuman.mat')


% Preparing a gene_exp_file for building ftINIT Models
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

threshold=[];
for i=1:size(met_gene_exp_union,2)
    x=met_gene_exp_union(:,i);
    threshold(:,i)=min(x(x>0));
    clear x
end

i=1; % For this minimum value (0.0333) was not working, therefore here we used second last minimum value
% Error using dispEM (line 49)
% Failed to find good enough solution within the time frame. MIPGap: 0.088757
x=met_gene_exp_union(:,i);
y=sort(x(x>0));
threshold(:,i)=y(3,1); % from sort(x(x>0));


%% Building models using ftINIT
data_structure.genes=ihuman.genes;
tissues={'P1_Liver_clst1', 'P1_Liver_clst2', 'P1_Breast_clst1', 'P1_Breast_clst2', ...
    'P2_Liver_clst1', 'P2_Liver_clst2', 'P2_Breast_clst1', 'P2_Breast_clst2', ...
    'P3_Bone_clst1', 'P3_Bone_clst2', 'P3_Brain_clst1', 'P3_Brain_clst2', 'P3_Breast_clst1', 'P3_Breast_clst2'};
data_structure.tissues=tissues;
% met_gene_exp_union1=met_gene_exp_union;
% met_gene_exp_union1(met_gene_exp_union1==0)=min(met_gene_exp_union1(met_gene_exp_union1>0));
data_structure.levels=met_gene_exp_union;

for i=1:size(met_gene_exp_union,2)
    prepData=prepDataHumanGEM; % prepData_ihuman;
    celltype=[];
    hpaData=[];
    transcrData=data_structure;
    th=threshold(1,i);
    transcrData.threshold=th;
    metabolomicsData=[];
    INITSteps=getINITSteps([],'1+1');
    removeGenes=true;
    
    useScoresForTasks=true;
    % paramsFT=[];
    paramsFT.TimeLimit=10000;
    paramsFT.MIPGap=0.01;
    paramsFT.OutputFlag=1;
    verbose=true;
    
    tissue=transcrData.tissues{1,i};
    
    [model, metProduction, addedRxnsForTasks, deletedRxnsInINIT, fullMipRes] = ftINIT(prepData, tissue, celltype, ...
        hpaData, transcrData, metabolomicsData, INITSteps, removeGenes, useScoresForTasks, paramsFT, verbose);
    
    models = ftINIT(prepData, data_structure.tissues{1}, [], [], data_structure, {}, getHumanGEMINITSteps('1+0'), false, true);
    
    str=[tissues{1,i},'_ftINIT','.model','=','model'];
    
    eval(str);
    
end


%% Saving them in structure files
n=0;
for i=1:7
    category={'P1_Liver_ftINIT', 'P1_Breast_ftINIT', 'P2_Liver_ftINIT', 'P2_Breast_ftINIT', ...
        'P3_Bone_ftINIT', 'P3_Brain_ftINIT', 'P3_Breast_ftINIT'};
    
    for j=1:2
        category2={'clst1','clst2'};
        n=n+1;
        
        category3={'P1_Liver_clst1_ftINIT', 'P1_Liver_clst2_ftINIT', 'P1_Breast_clst1_ftINIT', 'P1_Breast_clst2_ftINIT', ...
            'P2_Liver_clst1_ftINIT', 'P2_Liver_clst2_ftINIT', 'P2_Breast_clst1_ftINIT', 'P2_Breast_clst2_ftINIT', ...
            'P3_Bone_clst1_ftINIT', 'P3_Bone_clst2_ftINIT', 'P3_Brain_clst1_ftINIT', 'P3_Brain_clst2_ftINIT', ...
            'P3_Breast_clst1_ftINIT', 'P3_Breast_clst2_ftINIT'};
        
        str=[category{i},'.',category2{j},'=',category3{n}];
        
        eval(str);
        
    end
end

%% Checking Metabolic Tasks

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\ftINIT_Models\P1_Liver_ftINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\ftINIT_Models\P1_Breast_ftINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\ftINIT_Models\P2_Liver_ftINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\ftINIT_Models\P2_Breast_ftINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\ftINIT_Models\P3_Bone_ftINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\ftINIT_Models\P3_Brain_ftINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\ftINIT_Models\P3_Breast_ftINIT.mat')

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\taskStruct.mat')

category={};
for i=1:14
    category={'P1_Liver_ftINIT.clst1', 'P1_Liver_ftINIT.clst2', 'P1_Breast_ftINIT.clst1', 'P1_Breast_ftINIT.clst2', ...
        'P2_Liver_ftINIT.clst1', 'P2_Liver_ftINIT.clst2', 'P2_Breast_ftINIT.clst1', 'P2_Breast_ftINIT.clst2', ...
        'P3_Bone_ftINIT.clst1', 'P3_Bone_ftINIT.clst2', 'P3_Brain_ftINIT.clst1', 'P3_Brain_ftINIT.clst2', ...
        'P3_Breast_ftINIT.clst1', 'P3_Breast_ftINIT.clst2'};
    
    model=eval(category{i});
    [taskReport, ~, ~, ~]=checkTasks(model.model, [], true, false, false, taskStruct);
    
    results=taskReport.description(find(taskReport.ok));
    
    str=[category{i},'.task_report', '=','results'];
    
    eval(str);
    
    clear model taskReport str results
    
end


%% Performing Flux variability analysis (FVA)

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\ftINIT_Models\P1_Liver_ftINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\ftINIT_Models\P1_Breast_ftINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\ftINIT_Models\P2_Liver_ftINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\ftINIT_Models\P2_Breast_ftINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\ftINIT_Models\P3_Bone_ftINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\ftINIT_Models\P3_Brain_ftINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\ftINIT_Models\P3_Breast_ftINIT.mat')

% Initiating the Cobra Toolbox
initCobraToolbox(false)

category={'P1_Liver_ftINIT.clst1', 'P1_Liver_ftINIT.clst2', 'P1_Breast_ftINIT.clst1', 'P1_Breast_ftINIT.clst2', ...
    'P2_Liver_ftINIT.clst1', 'P2_Liver_ftINIT.clst2', 'P2_Breast_ftINIT.clst1', 'P2_Breast_ftINIT.clst2', ...
    'P3_Bone_ftINIT.clst1', 'P3_Bone_ftINIT.clst2', 'P3_Brain_ftINIT.clst1', 'P3_Brain_ftINIT.clst2', ...
    'P3_Breast_ftINIT.clst1', 'P3_Breast_ftINIT.clst2'};

for i=2:14
    model=eval(category{i});
    % osenseStr='max';  % --default
    [minFlux, maxFlux] = fluxVariability(model.model);
    
    str=[category{i},'.FVA', '=','[minFlux, maxFlux]'];
    eval(str);
    
    clear 'model' 'minFlux' 'maxFlux' 'str'
end


%% Metabolic functionality

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\ftINIT_Models\P1_Liver_ftINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\ftINIT_Models\P1_Breast_ftINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\ftINIT_Models\P2_Liver_ftINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\ftINIT_Models\P2_Breast_ftINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\ftINIT_Models\P3_Bone_ftINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\ftINIT_Models\P3_Brain_ftINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\ftINIT_Models\P3_Breast_ftINIT.mat')

category={'P1_Liver_ftINIT.clst1', 'P1_Liver_ftINIT.clst2', 'P1_Breast_ftINIT.clst1', 'P1_Breast_ftINIT.clst2', ...
    'P2_Liver_ftINIT.clst1', 'P2_Liver_ftINIT.clst2', 'P2_Breast_ftINIT.clst1', 'P2_Breast_ftINIT.clst2', ...
    'P3_Bone_ftINIT.clst1', 'P3_Bone_ftINIT.clst2', 'P3_Brain_ftINIT.clst1', 'P3_Brain_ftINIT.clst2', ...
    'P3_Breast_ftINIT.clst1', 'P3_Breast_ftINIT.clst2'};

for i=4:14
    model=eval(category{i});
    % osenseStr='max';  % --default
    [totTest, numTests] = testFunctionality_mouse(model.model);
    
    str=[category{i},'.met_fun', '=','[totTest, numTests]'];
    eval(str);
    
    clear model totTest numTests str
end

% save P1_Liver_ftINIT P1_Liver_ftINIT
% save P1_Breast_ftINIT P1_Breast_ftINIT
% save P2_Liver_ftINIT P2_Liver_ftINIT
% save P2_Breast_ftINIT P2_Breast_ftINIT
% save P3_Bone_ftINIT P3_Bone_ftINIT
% save P3_Brain_ftINIT P3_Brain_ftINIT
% save P3_Breast_ftINIT P3_Breast_ftINIT

