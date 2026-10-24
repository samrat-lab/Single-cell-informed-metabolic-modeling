clear
clc

%%%% Initiating the Raven Toolbox
% initCobraToolbox(false)
checkInstallation
setRavenSolver('gurobi')

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Human-GEM-main\model\Human-GEM.mat')

% Preparing a gene_exp_file for building tINIT Models
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE235168_variables\Primary_breast_met_data_clst1.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE235168_variables\Primary_breast_met_data_clst2.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE235168_variables\Final_met_genes.mat')

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


%% Extract information from the table into a structure called data_structure
data_structure.genes=ihuman.genes;
tissues={'Primary_breast_met_data_clst1', 'Primary_breast_met_data_clst2'};
data_structure.tissues=tissues;
data_structure.levels=met_gene_exp_union;
data_structure.threshold=min(threshold);


%% Add boundary metabolites to the Human-GEM model
ihuman = addBoundaryMets(ihuman);

% Boundary metabolites were added to 1658 reactions.
% New (boundary) versions of 1658 metabolites were added to the model.


%% Prepare the metabolic task list
essentialTasks=parseTaskList('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Human-GEM-main\data\metabolicTasks\metabolicTasks_Essential.txt');
% checkTasks(ihuman, [], true, false, false, essentialTasks); % All 57 tasks were present


%% Building models using tINIT
% refModel = ihuman;          % the reference model from which the GEM will be extracted
% tissue = 'liver';           % must match the tissue name in data_struct.tissues
% celltype = [];              % used if tissues are subdivided into cell type, which is not the case here
% hpaData = [];               % data structure containing protein abundance information (not used here)
% arrayData = data_struct;    % data structure with gene (RNA) abundance information
% metabolomicsData = [];      % list of metabolite names if metabolomic data is available
% removeGenes = true;         % (default) remove lowly/non-expressed genes from the extracted GEM
% taskFile = [];              % we already loaded the task file, so this input is not required
% useScoresForTasks = true;   % (default) use expression data to decide which reactions to keep
% printReport = true;         % (default) print status/completion report to screen
% taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
% params = [];                % additional optimization parameters for the INIT algorithm
% paramsFT = [];              % additional optimization parameters for the fit-tasks algorithm

category={'clst1','clst2'};

for i=1:size(met_gene_exp_union,2)
    refModel=ihuman; 
    tissue=data_structure.tissues{i};
    celltype=[];
    hpaData=[];
    arrayData=data_structure;
    metabolomicsData=[];
    removeGenes=true;
    taskFile=[];
    useScoresForTasks=true;
    printReport=true;
    taskStructure=essentialTasks;
    params=[];
    paramsFT=[];
        
    model=getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);

    str=['Primary_breast_tINIT','.',category{i},'.model','=','model'];
    
    eval(str);
    
end


%% We next verified the tINIT models, and observed few structural issues in the tINIT-generated model ...
%%% that may cause downstream functions like gpSampler or solveCobraLP to fail.
%%% 1. Missing Field: rules; 2. Inconsistent Field: b (Size of b is 2 but expected is 1)

category={'Primary_breast_tINIT.clst1', 'Primary_breast_tINIT.clst2'};

for i=1:2
    model=eval([category{i},'.model']);
    verifyModel(model); % To verify models and find out missing and inconsistent fields. 
    
    model = buildRxnGeneMat(model); % This will add the missing field 'rules' to the models.
    model.b = zeros(length(model.mets), 1); % This will correct the inconsistent field.
    
    str=[category{i},'.model', '=','model'];
    
    eval(str);
    
    clear model str
     
end

save Primary_breast_tINIT Primary_breast_tINIT


%% Checking Metabolic Tasks
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE235168_variables\Primary_breast_tINIT.mat')
essentialTasks=parseTaskList('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Human-GEM-main\data\metabolicTasks\metabolicTasks_Essential.txt');

category={'Primary_breast_tINIT.clst1', 'Primary_breast_tINIT.clst2'};
for i=1:2 
    model=eval(category{i});
    [taskReport, ~, ~, ~]=checkTasks(model.model, [], true, false, false, essentialTasks);
    
    results=taskReport.description(find(taskReport.ok));
    
    str=[category{i},'.task_report', '=','results'];
    eval(str);
    
    clear model taskReport str results
    
end

% save Primary_breast_tINIT Primary_breast_tINIT


% %% Performing Flux variability analysis (FVA)
% initCobraToolbox(false)
% 
% category={'Primary_breast_tINIT.clst1', 'Primary_breast_tINIT.clst2'};
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


% %% Metabolic functionality
% category={'Primary_breast_tINIT.clst1', 'Primary_breast_tINIT.clst2'};
% for i=1:2
%     model=eval(category{i});
%     % osenseStr='max';  % --default
%     [totTest, numTests] = testFunctionality_human(model.model);
%     
%     str=[category{i},'.met_fun', '=','[totTest, numTests]'];
%     eval(str);
%     
%     clear model totTest numTests str
% end
% 
% save('Primary_breast_tINIT.mat', 'Primary_breast_tINIT', '-v7.3')

