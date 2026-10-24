clear
clc

initCobraToolbox()
changeCobraSolver('gurobi', 'LP')

%% Load the models
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\Primary_breast_ftINIT.mat')

category={'Primary_breast_ftINIT.clst1', 'Primary_breast_ftINIT.clst2'};
for i=1:2
model=P1_Liver_ftINIT.clst1.model;
warmupPts=2*size(P1_Liver_ftINIT.clst1.model.rxns,1);
nSamples=2*size(P1_Liver_ftINIT.clst1.model.rxns,1);
nSteps=100;
nThreads=1;
solverName='gurobi';
verbose=1;
sModel = optGpSampler(model, warmupPts, nSamples, nSteps, nThreads, solverName, verbose);
str=['Primary_breast_ftINIT']
end

