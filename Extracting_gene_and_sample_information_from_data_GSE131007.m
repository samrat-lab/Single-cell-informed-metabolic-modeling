%%%% Here we are importing data normalised using CLR normalization technique, because in R we already integrated 
%%%% the three patients data

clear
clc

%% Import the data (Cell ID and Genes)
% Cell ID/Patient ID and organ
[~, ~, tempRow] = xlsread('D:\CANCER\Breast_Organotropism\Data_files\Data_normalised_using_diff_techniques\counts_matrix_normclr.xlsx', ...
    'counts_matrix_normclr', '1:1');

firstRow = string(tempRow);
firstRow(ismissing(firstRow)) = '';
firstRow(1) = [];

nonEmptyIdx = find(strlength(firstRow) > 0);
sample_GSE131007 = firstRow(1,nonEmptyIdx)';

temp=strings();
for i=1:size(sample_GSE131007,1)
    regexp(sample_GSE131007(i,1),'_','split','once');
    temp(i,1)=ans(1,1);
end

uni=unique(temp(:,1));

patient_ID_GSE131007 = strings();
for i=1:length(temp)
    if temp(i,1) == uni(1,1)
        patient_ID_GSE131007(i,1)='P1_Liver'; % UCD46L
    elseif temp(i,1) == uni(2,1)
        patient_ID_GSE131007(i,1)='P1_Breast'; % UCD46T
    elseif temp(i,1) == uni(3,1)
        patient_ID_GSE131007(i,1)='P2_Liver'; % UCD4L
    elseif temp(i,1) == uni(4,1)
        patient_ID_GSE131007(i,1)='P2_Breast'; % UCD4T
    elseif temp(i,1) == uni(5,1)
        patient_ID_GSE131007(i,1)='P3_Bone'; % UCD65Bo
    elseif temp(i,1) == uni(6,1)
        patient_ID_GSE131007(i,1)='P3_Brain'; % UCD65Br
    else
        patient_ID_GSE131007(i,1)='P3_Breast'; % UCD65T
    end
end

% Gene information
[~, ~, tempCol] = xlsread('D:\CANCER\Breast_Organotropism\Data_files\Data_normalised_using_diff_techniques\counts_matrix_normclr.xlsx', ...
    'counts_matrix_normclr', 'A:A');

firstCol = string(tempCol);
firstCol(ismissing(firstCol)) = '';
firstCol(1) = [];

nonEmptyIdx = find(strlength(firstCol) > 0);
genes_GSE131007 = firstCol(nonEmptyIdx,1);
genes_GSE131007=unique(genes_GSE131007);
genes_GSE131007(1:23)=[]; % Removed first 23 rows since they contain genes in the form of dates

clearvars -except patient_ID_GSE131007 genes_GSE131007 sample_GSE131007

save sample_GSE131007 sample_GSE131007
save patient_ID_GSE131007 patient_ID_GSE131007
save genes_GSE131007 genes_GSE131007


%% Import cell type information
[~, ~, cell_type_GSE131007] = xlsread('D:\CANCER\Breast_Organotropism\Data_files\Data_normalised_using_diff_techniques\label_counts_matrix_normclr_.xlsx','label_counts_matrix_normclr_');

cell_type_GSE131007 = string(cell_type_GSE131007);
cell_type_GSE131007(ismissing(cell_type_GSE131007)) = '';
cell_type_GSE131007(1,:)=[];
cell_type_GSE131007=[patient_ID_GSE131007 cell_type_GSE131007];

save cell_type_GSE131007 cell_type_GSE131007


%% Filtering epithelial cells for all patients
idx=find(cell_type_GSE131007(:,3)=='Epithelial_cells');
Epi_cells_GSE131007=cell_type_GSE131007(idx,:);

save Epi_cells_GSE131007 Epi_cells_GSE131007


%% Filtering metabolic genes
% Removing . from the genes name since ihuman have genes without dots 
% But the genes with dots are not available present in the data GSE131007
genes1=strings();
for i=1:size(genes_GSE131007,1)
    strsplit(genes_GSE131007(i,1), '.');
    genes1(i,1)=strjoin(ans,'');
clear ans
end

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Human-GEM.mat')
[~,p,q]=intersect(upper(string(ihuman.geneShortNames(:,1))),upper(genes1(:,1))); % 2480 metabolic genes out of 2889 are present

Met_genes_GSE131007=[string(ihuman.genes(p,:)),string(ihuman.geneShortNames(p,:))];

save Met_genes_GSE131007 Met_genes_GSE131007
