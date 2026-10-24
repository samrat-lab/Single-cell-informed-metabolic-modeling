%% Import the data
[~, ~, P3_UCD65] = xlsread('D:\CANCER\Breast_Organotropism\Data_files\organotropism_data_18december\udc65\GSE131007_ucd65_count_matrix.xlsx','GSE131007_ucd65_count_matrix');

P3_UCD65 = string(P3_UCD65);
P3_UCD65(ismissing(P3_UCD65)) = '';

save('P3_UCD65.mat', 'P3_UCD65', '-v7.3')

%% Filtering data of met_genes and epithelial_cells for P2_Liver and P2_Breast
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Met_genes_GSE131007.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Epi_cells_GSE131007.mat')

% P3_Bone
idx=find(Epi_cells_GSE131007(:,1)=='P3_Bone');
tmp=Epi_cells_GSE131007(idx,2);

x=P3_UCD65(2:end,2:end);
y=P3_UCD65(1,2:end);
z=P3_UCD65(2:end,1);

[p1,~,p]=intersect(tmp,y);
tmp1=x(:,p); % 501 Epithelial cells data filtered

[q,~,r]=intersect(upper(Met_genes_GSE131007(:,2)),upper(z)); % 2200 Metabolic genes filtered
Raw_P3_Bone=[["Met_Genes",p1'];[q,tmp1(r,:)]];

save('Raw_P3_Bone.mat', 'Raw_P3_Bone', '-v7.3')

clear idx tmp x y z tmp1 p1 p q r

% P3_Brain
idx=find(Epi_cells_GSE131007(:,1)=='P3_Brain');
tmp=Epi_cells_GSE131007(idx,2);

x=P3_UCD65(2:end,2:end);
y=P3_UCD65(1,2:end);
z=P3_UCD65(2:end,1);

[p1,~,p]=intersect(tmp,y);
tmp1=x(:,p); % 90 Epithelial cells data filtered

[q,~,r]=intersect(upper(Met_genes_GSE131007(:,2)),upper(z)); % 2200 Metabolic genes filtered
Raw_P3_Brain=[["Met_Genes",p1'];[q,tmp1(r,:)]];

save('Raw_P3_Brain.mat', 'Raw_P3_Brain', '-v7.3')

clear idx tmp x y z tmp1 p1 p q r

% P3_Breast
idx=find(Epi_cells_GSE131007(:,1)=='P3_Breast');
tmp=Epi_cells_GSE131007(idx,2);

x=P3_UCD65(2:end,2:end);
y=P3_UCD65(1,2:end);
z=P3_UCD65(2:end,1);

[p1,~,p]=intersect(tmp,y);
tmp1=x(:,p); % 1128 Epithelial cells data filtered

[q,~,r]=intersect(upper(Met_genes_GSE131007(:,2)),upper(z)); % 2200 Metabolic genes filtered
Raw_P3_Breast=[["Met_Genes",p1'];[q,tmp1(r,:)]];

save('Raw_P3_Breast.mat', 'Raw_P3_Breast', '-v7.3')

clear idx tmp x y z tmp1 p1 p q r

