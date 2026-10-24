%% Import the data
[~, ~, P1_UCD46] = xlsread('D:\CANCER\Breast_Organotropism\Data_files\organotropism_data_18december\ucd46\GSE131007_ucd46_count_matrix.xlsx','GSE131007_ucd46_count_matrix');

P1_UCD46 = string(P1_UCD46);
P1_UCD46(ismissing(P1_UCD46)) = '';

save('P1_UCD46.mat', 'P1_UCD46', '-v7.3')

%% Filtering data of met_genes and epithelial_cells for P1_Liver and P1_Breast
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Met_genes_GSE131007.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Epi_cells_GSE131007.mat')

% P1_Liver
idx=find(Epi_cells_GSE131007(:,1)=='P1_Liver');
tmp=Epi_cells_GSE131007(idx,2);

x=P1_UCD46(2:end,2:end);
y=P1_UCD46(1,2:end);
z=P1_UCD46(2:end,1);

[p1,~,p]=intersect(tmp,y);
tmp1=x(:,p); % 223 Epithelial cells data filtered

[q,~,r]=intersect(upper(Met_genes_GSE131007(:,2)),upper(z)); % 2221 Metabolic genes filtered
Raw_P1_Liver=[["Met_Genes",p1'];[q,tmp1(r,:)]];

save('Raw_P1_Liver.mat', 'Raw_P1_Liver', '-v7.3')

clear idx tmp x y z tmp1 p1 p q r


% P1_Breast
idx=find(Epi_cells_GSE131007(:,1)=='P1_Breast');
tmp=Epi_cells_GSE131007(idx,2);

x=P1_UCD46(2:end,2:end);
y=P1_UCD46(1,2:end);
z=P1_UCD46(2:end,1);

[p1,~,p]=intersect(tmp,y);
tmp1=x(:,p); % 429 Epithelial cells data filtered

[q,~,r]=intersect(upper(Met_genes_GSE131007(:,2)),upper(z)); % 2221 Metabolic genes filtered
Raw_P1_Breast=[["Met_Genes",p1'];[q,tmp1(r,:)]];

save('Raw_P1_Breast.mat', 'Raw_P1_Breast', '-v7.3')

clear idx tmp x y z tmp1 p1 p q r

