clear
clc

%% Comparing Human1 and Mouse1 genes
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\Mouse-GEM.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Human-GEM-main\model\Human-GEM.mat')

%% Import the data
% UCD4 genes
[~, ~, ucd4_genes] = xlsread('D:\CANCER\Breast_Organotropism\Data_files\organotropism_data_18december\ucd4\GSE131007_ucd4_count_matrix.xlsx','GSE131007_ucd4_count_matrix','A1:A19087');

ucd4_genes = string(ucd4_genes);
ucd4_genes(ismissing(ucd4_genes)) = '';
ucd4_genes(1) = [];

% UCD46 genes
[~, ~, ucd46_genes] = xlsread('D:\CANCER\Breast_Organotropism\Data_files\organotropism_data_18december\ucd46\GSE131007_ucd46_count_matrix.xlsx','GSE131007_ucd46_count_matrix','A1:A16768');

ucd46_genes = string(ucd46_genes);
ucd46_genes(ismissing(ucd46_genes)) = '';
ucd46_genes(1) = [];

% UCD65 genes
[~, ~, ucd65_genes] = xlsread('D:\CANCER\Breast_Organotropism\Data_files\organotropism_data_18december\udc65\GSE131007_ucd65_count_matrix.xlsx','GSE131007_ucd65_count_matrix','A1:A16817');

ucd65_genes = string(ucd65_genes);
ucd65_genes(ismissing(ucd65_genes)) = '';
ucd65_genes(1) = [];


%% Common genes between Human1 and Mouse1 model
com=intersect(upper(string(ihuman.geneShortNames(:,1))),upper(string(mouseGEM.genes(:,1))));
% ihuman- 2889 genes, mouseGEM- 2846 genes, com- 2685 genes


%% Common genes between Human1 model and data
intersect(upper(string(ihuman.geneShortNames(:,1))),upper(ucd46_genes(:,1))); % 2221 genes
intersect(upper(string(ihuman.geneShortNames(:,1))),upper(ucd4_genes(:,1))); % 2386 genes
intersect(upper(string(ihuman.geneShortNames(:,1))),upper(ucd65_genes(:,1))); % 2200 genes


%% Common genes between Mouse1 model and data
intersect(upper(string(mouseGEM.genes(:,1))),upper(ucd46_genes(:,1))); % 2129 genes
intersect(upper(string(mouseGEM.genes(:,1))),upper(ucd4_genes(:,1))); % 2280 genes
intersect(upper(string(mouseGEM.genes(:,1))),upper(ucd65_genes(:,1))); % 2111 genes


%% Intersecting common genes of Human1 and Mouse1 model with data
intersect(upper(com),upper(ucd46_genes(:,1))); % 2127 genes
intersect(upper(com),upper(ucd4_genes(:,1))); % 2277 genes
intersect(upper(com),upper(ucd65_genes(:,1))); % 2109 genes


%% Unique genes
uniqueHuman=setdiff(upper(string(ihuman.geneShortNames(:,1))),upper(string(mouseGEM.genes(:,1)))); % 204 genes
uniqueMouse=setdiff(upper(string(mouseGEM.genes(:,1))),upper(string(ihuman.geneShortNames(:,1)))); % 161 genes
