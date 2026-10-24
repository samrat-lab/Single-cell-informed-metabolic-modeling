clear
clc

%% Import the data of mouse metabolites ID across different platforms- MouseGEM
[~, ~, Mouse_mets] = xlsread('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Mouse-GEM-main\metabolites.xlsx','metabolites');

Mouse_mets = string(Mouse_mets);
Mouse_mets(ismissing(Mouse_mets)) = '';


%% Import the data of mouse metabolites ID across different platforms- Human1
[~, ~, Human1_mets] = xlsread('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Human-GEM-main\model\Human_GEM_mets.xlsx','metabolites');

Human1_mets = string(Human1_mets);
Human1_mets(ismissing(Human1_mets)) = '';


%% Separating compartment information from mouseGEM.mets
comp=strings();
y=Mouse_mets(2:end,1);
for i=1:size(y,1)
    tmp=convertStringsToChars(y(i,1));
    comp(i,1)=tmp(end);
    clear tmp
end

%% Finding corresponding metabolite IDs in Human1 from Recon3D for use in testFunctionality_human
% For testFunctionality_human file
category={'atp','adp','pi','h2o','h'};
Metabolites=strings(0,1);
for i=1:size(category,2)
    idx=find(Human1_mets(:,11)==category{1,i});
    tmp=Human1_mets(idx,1:2);
    idx1=endsWith(tmp(:,1), 'c');
    Metabolites=[Metabolites; tmp(idx1,1)];
    clear idx tmp idx1
end
Metabolites=[Metabolites,string(category)'];


%% Building biomass_rxn excel file for Human1 model to test functionality
% Import the data
[~, ~, biomass_rxn] = xlsread('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\biomass_rxn.xls','Sheet1','A1:A77');

biomass_rxn = string(biomass_rxn);
biomass_rxn(ismissing(biomass_rxn)) = '';
biomass_rxn=biomass_rxn(1:60,1);

biomass_rxn_human=strings(0,1);
for i=1:size(biomass_rxn,1)
    x=extractBefore(biomass_rxn{i,1},'[');
    idx=find(Human1_mets(:,11)==x);
    tmp=Human1_mets(idx,1:2);
    idx1=endsWith(tmp(:,1), 'c');
    biomass_rxn_human=[biomass_rxn_human; tmp(idx1,1)];
    clear idx tmp idx1
end
biomass_rxn_human=[biomass_rxn_human,biomass_rxn];

xlswrite('biomass_rxn_human.xlsx', biomass_rxn_human(:,1));


