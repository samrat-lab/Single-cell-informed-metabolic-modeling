%% Import the data
[~, ~, P2_UCD4] = xlsread('D:\CANCER\Breast_Organotropism\Data_files\organotropism_data_18december\ucd4\GSE131007_ucd4_count_matrix.xlsx','GSE131007_ucd4_count_matrix');

P2_UCD4 = string(P2_UCD4);
P2_UCD4(ismissing(P2_UCD4)) = '';

save('P2_UCD4.mat', 'P2_UCD4', '-v7.3')

%% Filtering data of met_genes and epithelial_cells for P2_Liver and P2_Breast
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Met_genes_GSE131007.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Epi_cells_GSE131007.mat')

% P2_Liver
idx=find(Epi_cells_GSE131007(:,1)=='P2_Liver');
tmp=Epi_cells_GSE131007(idx,2);

x=P2_UCD4(2:end,2:end);
y=P2_UCD4(1,2:end);
z=P2_UCD4(2:end,1);

[p1,~,p]=intersect(tmp,y);
tmp1=x(:,p); % 2799 Epithelial cells data filtered

[q,~,r]=intersect(upper(Met_genes_GSE131007(:,2)),upper(z)); % 2386 Metabolic genes filtered
Raw_P2_Liver=[["Met_Genes",p1'];[q,tmp1(r,:)]];

save('Raw_P2_Liver.mat', 'Raw_P2_Liver', '-v7.3')

clear idx tmp x y z tmp1 p1 p q r


% P2_Breast
idx=find(Epi_cells_GSE131007(:,1)=='P2_Breast');
tmp=Epi_cells_GSE131007(idx,2);

x=P2_UCD4(2:end,2:end);
y=P2_UCD4(1,2:end);
z=P2_UCD4(2:end,1);

[p1,~,p]=intersect(tmp,y);
tmp1=x(:,p); % 3144 Epithelial cells data filtered

[q,~,r]=intersect(upper(Met_genes_GSE131007(:,2)),upper(z)); % 2386 Metabolic genes filtered
Raw_P2_Breast=[["Met_Genes",p1'];[q,tmp1(r,:)]];

save('Raw_P2_Breast.mat', 'Raw_P2_Breast', '-v7.3')

clear idx tmp x y z tmp1 p1 p q r


