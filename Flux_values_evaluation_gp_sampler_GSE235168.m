clear
clc

initCobraToolbox(false)
changeCobraSolver('gurobi', 'LP')
gurobi_setup

% Load the model
% load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE235168_variables\Primary_breast_ftINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE235168_variables\Primary_breast_tINIT.mat')


%% Using GP sampler
% category={'Primary_breast_ftINIT.clst1', 'Primary_breast_ftINIT.clst2'};
category={'Primary_breast_tINIT.clst1', 'Primary_breast_tINIT.clst2'};

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
    
    clear i model nPoints storedata str
end

% save('Primary_breast_ftINIT.mat', 'Primary_breast_ftINIT', '-v7.3')
save('Primary_breast_tINIT.mat', 'Primary_breast_tINIT', '-v7.3')


