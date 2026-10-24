clear
clc

initCobraToolbox(false)
changeCobraSolver('gurobi', 'LP')
gurobi_setup

% Load the models
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P1_Liver_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P1_Breast_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P2_Liver_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P2_Breast_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P3_Bone_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P3_Brain_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P3_Breast_tINIT.mat')


%% Using GP sampler
category={'P1_Liver_tINIT.clst1', 'P1_Liver_tINIT.clst2', 'P1_Breast_tINIT.clst1', 'P1_Breast_tINIT.clst2', ...
    'P2_Liver_tINIT.clst1', 'P2_Liver_tINIT.clst2', 'P2_Breast_tINIT.clst1', 'P2_Breast_tINIT.clst2', ...
    'P3_Bone_tINIT.clst1', 'P3_Bone_tINIT.clst2', 'P3_Brain_tINIT.clst1', 'P3_Brain_tINIT.clst2', ...
    'P3_Breast_tINIT.clst1', 'P3_Breast_tINIT.clst2'};

for i=1:size(category,2)
    
    model=eval([category{i},'.model']);
    %     nPoints=2*size(model.rxns,1);
    %     warmuppts=createHRWarmup(model, nPoints);
    %     model.warmupPts=warmuppts;
    nPoints=2*size(model.rxns,1);
    %     bias=[];
    %     maxTime=600;
    %     maxSteps=1e10;
    %     threads=10;
    storedata=gpSampler(model,nPoints);%,bias,maxTime,maxSteps,threads);
    
    str=[category{i},'.storedata_gpSampler', '=','storedata'];
    eval(str);
    
    i
    
    clear model nPoints storedata str
end

% filePath = 'your/full/path/to/sampleStructTemp.mat';
% if exist(filePath, 'file')
%     delete(filePath);
% end


save('P1_Liver_tINIT.mat', 'P1_Liver_tINIT', '-v7.3')
save('P1_Breast_tINIT.mat', 'P1_Breast_tINIT', '-v7.3')
save('P2_Liver_tINIT.mat', 'P2_Liver_tINIT', '-v7.3')
save('P2_Breast_tINIT.mat', 'P2_Breast_tINIT', '-v7.3')
save('P3_Bone_tINIT.mat', 'P3_Bone_tINIT', '-v7.3')
save('P3_Brain_tINIT.mat', 'P3_Brain_tINIT', '-v7.3')
save('P3_Breast_tINIT.mat', 'P3_Breast_tINIT', '-v7.3')


