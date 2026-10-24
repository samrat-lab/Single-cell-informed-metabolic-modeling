clear
clc

%% Load the mat data saved in R
load('D:\CANCER\Breast_Organotropism\Primary_breast_cancer_data\PDX_data\GSE235168_RAW\epithelial_chunk1.mat')
load('D:\CANCER\Breast_Organotropism\Primary_breast_cancer_data\PDX_data\GSE235168_RAW\epithelial_chunk2.mat')

%% Import Ensemble_gene_ID, gene_entrez_IDs, gene_symbols and epithelial_cell_labels data
% Ensemble_gene_ID
filename = 'D:\CANCER\Breast_Organotropism\Primary_breast_cancer_data\PDX_data\GSE235168_RAW\Ensemble_gene_ID.csv';
delimiter = ',';
formatSpec = '%q%q%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);
ensemble_gene_ID = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;


% Gene_entrez_IDs
filename = 'D:\CANCER\Breast_Organotropism\Primary_breast_cancer_data\PDX_data\GSE235168_RAW\gene_entrez_IDs.csv';
delimiter = ',';
formatSpec = '%q%q%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);
gene_entrez_IDs = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;

gene_entrez_IDs(find(gene_entrez_IDs(:,2)=='NA'),2)='';

% Gene_symbols
filename = 'D:\CANCER\Breast_Organotropism\Primary_breast_cancer_data\PDX_data\GSE235168_RAW\gene_symbols.csv';
delimiter = ',';
formatSpec = '%q%q%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);
gene_symbols = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;


% Epithelial_cell_labels
filename = 'D:\CANCER\Breast_Organotropism\Primary_breast_cancer_data\PDX_data\GSE235168_RAW\epithelial_cells.csv';
delimiter = ',';
formatSpec = '%q%q%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);
epithelial_cells = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;


%% Combining the complete data
Complete_data=[["Ensemble_ID","Gene_symbol","Gene_entrez_ID",epithelial_cells(2:end,1)'];...
    [ensemble_gene_ID(2:end,2),gene_symbols(2:end,2),gene_entrez_IDs(2:end,2),epithelial_chunk1,epithelial_chunk2]];

% clearvars -except Complete_data
% save('Complete_data.mat', 'Complete_data', '-v7.3')


%% Load the Human1 model and filter out the metabolic genes data
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\Complete_data.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Human-GEM-main\model\Human-GEM.mat')
[met_genes,p,q]=intersect(string(ihuman.genes(:,1)),Complete_data(:,1));

Filtered_data=[Complete_data(1,:);Complete_data(q,:)];
% save('Filtered_data.mat', 'Filtered_data', '-v7.3')


%% Data pre-processing
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\Filtered_data.mat')
x=double(Filtered_data(2:end,4:end));

%%%% Cells filteration (Epithelial cells expressing 20% of the metabolic genes)
d1=strings();
for i=1:size(x,2)
    find(double(x(:,i)));
    d1(1,i)=length(ans);
end

z=Filtered_data(1,4:end);
p1=strings();
d2=[];
p1=z(1,find(double(d1)>0.2*(size(x,1)))); % Filtered cells ID
d2=x(:,find(double(d1)>0.2*(size(x,1)))); % Filtered data


%%%% Genes filteration (Met genes present in atleast 20% of the cells)
g1=strings();
for i=1:size(d2,1)
    find(d2(i,:));
    g1(i,1)=length(ans);
end

y=Filtered_data(2:end,1);
p2=strings();
g2=[];
p2=y(find(double(g1)>0.2*(size(d2,2))),1); % Met genes ensemble ID
g2=d2(find(double(g1)>0.2*(size(d2,2))),:); % Filtered met genes data

Final_met_genes=p2;
% save Final_met_genes Final_met_genes

clear d1 g1 x y z i ans

%%%% Scaling
d3=[];
for i=1:size(g2,1)
    for j=1:size(g2,2)
        d3(i,j)=log(g2(i,j) + 1); 
    end
end

find(isnan(d3)) % empty array (No NAN values)

%%%% Finding outliers and fixing it
d4=[];
for i=1:size(d3,1)
    tmp=d3(i,:);
    quant = quantile(tmp,[0.025 0.25 0.50 0.75 0.975]);
    IQR = quant(1,4) - quant(1,2);
    
    lower_fence = quant(1,2)-(1.5)*IQR;
    upper_fence = quant(1,4)+(1.5)*IQR;
    
    tmp(tmp<lower_fence)=quant(1,1);
    tmp(tmp>upper_fence)=quant(1,5);
    d4(i,:)=tmp;
    
    clear 'quant' 'IQR' 'tmp' 'lower_fence' 'upper_fence'
end

find(isnan(d4)) % empty array (No NAN values)

%%%% Gene-wise normalisation
d5=[];
for i=1:size(d4,1)
    tmp=d4(i,:);
    m=max(tmp);
    d5(i,:)=tmp/m;
    clear tmp
end

find(isnan(d5)) % empty array (No NAN values)

%%%% Evaluating no. of clusters
rng('default') % For reproducibility
eva = evalclusters(d5,'kmeans','CalinskiHarabasz','KList',1:10);

%%%% Making clusters of the data
[idx,C]=kmeans(d5',eva.OptimalK);

%%%% Cluster1
Primary_breast_clst1_cell_ID=p1(find(idx==1));
Primary_breast_met_data_clst1=C(1,:);

%%%% Cluster2
Primary_breast_clst2_cell_ID=p1(find(idx==2));
Primary_breast_met_data_clst2=C(2,:);

% save Primary_breast_clst1_cell_ID Primary_breast_clst1_cell_ID
% save Primary_breast_clst2_cell_ID Primary_breast_clst2_cell_ID
% save Primary_breast_met_data_clst1 Primary_breast_met_data_clst1
% save Primary_breast_met_data_clst2 Primary_breast_met_data_clst2

plot(Primary_breast_met_data_clst1,'.r')
hold on
plot(Primary_breast_met_data_clst2,'.b')

