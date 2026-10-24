clear
clc

% rmpath(genpath('D:\CANCER\Breast_Organotropism\RAVEN-2.10.3'))
% savepath
% 
% addpath(genpath('D:\CANCER\Breast_Organotropism\RAVEN-2.7.9'))
% savepath

addpath(genpath('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\Human-GEM-main'))
savepath

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Human-GEM-main\model\Human-GEM.mat')

model=ihuman;
convertGenes=false;
essentialTasksFilePath='D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\Human-GEM-main\data\metabolicTasks\metabolicTasks_Essential.txt';
rxnsFilePath='D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\Human-GEM-main\model\reactions.tsv';
prepDataHumanGEM = prepHumanModelForftINIT(model, convertGenes, essentialTasksFilePath, rxnsFilePath);

save prepDataHumanGEM prepDataHumanGEM