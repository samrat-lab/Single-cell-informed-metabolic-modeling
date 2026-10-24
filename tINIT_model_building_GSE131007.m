clear
clc

%%%% Initiating the Raven Toolbox
% initCobraToolbox(false)
checkInstallation
setRavenSolver('gurobi')

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Human-GEM-main\model\Human-GEM.mat')


% Preparing a gene_exp_file for building tINIT Models
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

%% Extract information from the table into a structure called data_structure
data_structure.genes=ihuman.genes;
tissues={'P1_Liver_clst1', 'P1_Liver_clst2', 'P1_Breast_clst1', 'P1_Breast_clst2', ...
    'P2_Liver_clst1', 'P2_Liver_clst2', 'P2_Breast_clst1', 'P2_Breast_clst2', ...
    'P3_Bone_clst1', 'P3_Bone_clst2', 'P3_Brain_clst1', 'P3_Brain_clst2', 'P3_Breast_clst1', 'P3_Breast_clst2'};
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
    
    str=[tissues{1,i},'_tINIT','.model','=','model'];
    
    eval(str);
    
end


%% Saving them in structure files
% P1 Liver
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\All_tINIT_models\P1_Liver_clst1_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\All_tINIT_models\P1_Liver_clst2_tINIT.mat')

% P1 Breast
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\All_tINIT_models\P1_Breast_clst1_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\All_tINIT_models\P1_Breast_clst2_tINIT.mat')

% P2 Liver
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\All_tINIT_models\P2_Liver_clst1_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\All_tINIT_models\P2_Liver_clst2_tINIT.mat')

% P2 Breast
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\All_tINIT_models\P2_Breast_clst1_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\All_tINIT_models\P2_Breast_clst2_tINIT.mat')

% P3 Bone
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\All_tINIT_models\P3_Bone_clst1_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\All_tINIT_models\P3_Bone_clst2_tINIT.mat')

% P3 Brain
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\All_tINIT_models\P3_Brain_clst1_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\All_tINIT_models\P3_Brain_clst2_tINIT.mat')

% P3 Breast
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\All_tINIT_models\P3_Breast_clst1_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\All_tINIT_models\P3_Breast_clst2_tINIT.mat')


n=0;
for i=1:7
    category={'P1_Liver_tINIT', 'P1_Breast_tINIT', 'P2_Liver_tINIT', 'P2_Breast_tINIT', ...
        'P3_Bone_tINIT', 'P3_Brain_tINIT', 'P3_Breast_tINIT'};
    
    for j=1:2
        category2={'clst1','clst2'};
        n=n+1;
        
        category3={'P1_Liver_clst1_tINIT', 'P1_Liver_clst2_tINIT', 'P1_Breast_clst1_tINIT', 'P1_Breast_clst2_tINIT', ...
            'P2_Liver_clst1_tINIT', 'P2_Liver_clst2_tINIT', 'P2_Breast_clst1_tINIT', 'P2_Breast_clst2_tINIT', ...
            'P3_Bone_clst1_tINIT', 'P3_Bone_clst2_tINIT', 'P3_Brain_clst1_tINIT', 'P3_Brain_clst2_tINIT', ...
            'P3_Breast_clst1_tINIT', 'P3_Breast_clst2_tINIT'};
        
        str=[category{i},'.',category2{j},'=',category3{n}];
        
        eval(str);
        
    end
end



%% We next verified the tINIT models, and observed few structural issues in the tINIT-generated model ...
%%% that may cause downstream functions like gpSampler or solveCobraLP to fail.
%%% 1. Missing Field: rules; 2. Inconsistent Field: b (Size of b is 2 but expected is 1)

category={'P1_Liver_tINIT.clst1', 'P1_Liver_tINIT.clst2', 'P1_Breast_tINIT.clst1', 'P1_Breast_tINIT.clst2', ...
    'P2_Liver_tINIT.clst1', 'P2_Liver_tINIT.clst2', 'P2_Breast_tINIT.clst1', 'P2_Breast_tINIT.clst2', ...
    'P3_Bone_tINIT.clst1', 'P3_Bone_tINIT.clst2', 'P3_Brain_tINIT.clst1', 'P3_Brain_tINIT.clst2', ...
    'P3_Breast_tINIT.clst1', 'P3_Breast_tINIT.clst2'};

for i=1:size(category,2)
    model=eval([category{i},'.model']);
    verifyModel(model); % To verify models and find out missing and inconsistent fields. 
    
    model = buildRxnGeneMat(model); % This will add the missing field 'rules' to the models.
    model.b = zeros(length(model.mets), 1); % This will correct the inconsistent field.
    
    str=[category{i},'.model', '=','model'];
    
    eval(str);
    
    clear model str
     
end


save P1_Liver_tINIT P1_Liver_tINIT
save P1_Breast_tINIT P1_Breast_tINIT
save P2_Liver_tINIT P2_Liver_tINIT
save P2_Breast_tINIT P2_Breast_tINIT
save P3_Bone_tINIT P3_Bone_tINIT
save P3_Brain_tINIT P3_Brain_tINIT
save P3_Breast_tINIT P3_Breast_tINIT


%% Checking Metabolic Tasks

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P1_Liver_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P1_Breast_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P2_Liver_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P2_Breast_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P3_Bone_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P3_Brain_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P3_Breast_tINIT.mat')

% Essential Tasks
essentialTasks=parseTaskList('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Human-GEM-main\data\metabolicTasks\metabolicTasks_Essential.txt');

category={};
for i=1:size(category,2)
    category={'P1_Liver_tINIT.clst1', 'P1_Liver_tINIT.clst2', 'P1_Breast_tINIT.clst1', 'P1_Breast_tINIT.clst2', ...
        'P2_Liver_tINIT.clst1', 'P2_Liver_tINIT.clst2', 'P2_Breast_tINIT.clst1', 'P2_Breast_tINIT.clst2', ...
        'P3_Bone_tINIT.clst1', 'P3_Bone_tINIT.clst2', 'P3_Brain_tINIT.clst1', 'P3_Brain_tINIT.clst2', ...
        'P3_Breast_tINIT.clst1', 'P3_Breast_tINIT.clst2'};
    
    model=eval(category{i});
    [taskReport, ~, ~, ~]=checkTasks(model.model, [], true, false, false, essentialTasks);
    
    results=taskReport.description(find(taskReport.ok));
    
    str=[category{i},'.task_report', '=','results'];
    
    eval(str);
    
    clear model taskReport str results
    
end

% save P1_Liver_tINIT P1_Liver_tINIT
% save P1_Breast_tINIT P1_Breast_tINIT
% save P2_Liver_tINIT P2_Liver_tINIT
% save P2_Breast_tINIT P2_Breast_tINIT
% save P3_Bone_tINIT P3_Bone_tINIT
% save P3_Brain_tINIT P3_Brain_tINIT
% save P3_Breast_tINIT P3_Breast_tINIT


%% Metabolic functionality
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P1_Liver_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P1_Breast_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P2_Liver_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P2_Breast_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P3_Bone_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P3_Brain_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P3_Breast_tINIT.mat')

category={'P1_Liver_tINIT.clst1', 'P1_Liver_tINIT.clst2', 'P1_Breast_tINIT.clst1', 'P1_Breast_tINIT.clst2', ...
    'P2_Liver_tINIT.clst1', 'P2_Liver_tINIT.clst2', 'P2_Breast_tINIT.clst1', 'P2_Breast_tINIT.clst2', ...
    'P3_Bone_tINIT.clst1', 'P3_Bone_tINIT.clst2', 'P3_Brain_tINIT.clst1', 'P3_Brain_tINIT.clst2', ...
    'P3_Breast_tINIT.clst1', 'P3_Breast_tINIT.clst2'};

for i=1:size(category,2)
    model=eval(category{i});
    % osenseStr='max';  % --default
    [totTest, numTests] = testFunctionality_human(model.model);
    
    str=[category{i},'.met_fun', '=','[totTest, numTests]'];
    eval(str);
    
    clear model totTest numTests str
end

% save P1_Liver_tINIT P1_Liver_tINIT
% save P1_Breast_tINIT P1_Breast_tINIT
% save P2_Liver_tINIT P2_Liver_tINIT
% save P2_Breast_tINIT P2_Breast_tINIT
% save P3_Bone_tINIT P3_Bone_tINIT
% save P3_Brain_tINIT P3_Brain_tINIT
% save P3_Breast_tINIT P3_Breast_tINIT

