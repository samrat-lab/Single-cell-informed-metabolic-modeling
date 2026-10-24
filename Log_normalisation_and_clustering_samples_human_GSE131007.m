%%%%% Normalising the data log(count + 1)/log(max count of the gene across the samples) %%%%%
clear
clc

%% Load the data
% P1_Liver
% load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Raw_P1_Liver.mat')

% % P1_Breast
% load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Raw_P1_Breast.mat')
% 
% % P2_Liver
% load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Raw_P2_Liver.mat')
% 
% % P2_Breast
% load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Raw_P2_Breast.mat')
% 
% % P3_Bone
% load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Raw_P3_Bone.mat')
% 
% % P3_Brain
% load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Raw_P3_Brain.mat')
% 
% % P3_Breast
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Raw_P3_Breast.mat')

% Load the Human1 model
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Human-GEM.mat')


%% Data pre-processing
data=Raw_P3_Breast;
x=double(data(2:end,2:end));

%%%% Cells filteration (Epithelial cells expressing 20% of the metabolic genes)
d1=strings();
for i=1:size(x,2)
    find(double(x(:,i)));
    d1(1,i)=length(ans);
end

z=data(1,2:end);
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

y=data(2:end,1);
p2=strings();
g2=[];
p2=y(find(double(g1)>0.2*(size(d2,2))),1); % Met gene names
g2=d2(find(double(g1)>0.2*(size(d2,2))),:); % Filtered met genes data

% Final_met_genes=p2;
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
P3_Breast_clst1_cell_ID=p1(find(idx==1));
P3_Breast_clst1=C(1,:)';

%%%% Cluster2
P3_Breast_clst2_cell_ID=p1(find(idx==2));
P3_Breast_clst2=C(2,:)';

%%%% Metabolic genes
P3_Breast_met_genes=p2;

save P3_Breast_clst1_cell_ID P3_Breast_clst1_cell_ID
save P3_Breast_clst1 P3_Breast_clst1
save P3_Breast_clst2_cell_ID P3_Breast_clst2_cell_ID
save P3_Breast_clst2 P3_Breast_clst2
save P3_Breast_met_genes P3_Breast_met_genes

% plot(Primary_breast_met_data_clst1,'.r')
% hold on
% plot(Primary_breast_met_data_clst2,'.b')

