clear
clc

% Load the secretory metabolites data patient wise
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\P1_Breast_secretory_mets.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\P2_Breast_secretory_mets.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\P3_Breast_secretory_mets.mat')

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Secretory_metabolites_and_reactions\Total_secretory_mets_flux_values_all.mat')

% P1 Breast
idx=find(double(Total_secretory_mets_flux_values_all(2:end,3))==0);
secretory_mets=Total_secretory_mets_flux_values_all(2:end,1);






















