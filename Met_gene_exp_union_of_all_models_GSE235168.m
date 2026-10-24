clear
clc

% Load Human-GEM
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Human-GEM-main\model\Human-GEM.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\Final_met_genes.mat')

[~,p1,q1]=intersect(string(ihuman.genes), P1_Liver_met_genes(:,1));

str1={'P1_Liver_clst1', 'P1_Liver_clst2', 'P1_Breast_clst1', 'P1_Breast_clst2', ...
    'P2_Liver_clst1', 'P2_Liver_clst2', 'P2_Breast_clst1', 'P2_Breast_clst2', ...
    'P3_Bone_clst1', 'P3_Bone_clst2', 'P3_Brain_clst1', 'P3_Brain_clst2', 'P3_Breast_clst1', 'P3_Breast_clst2'};

met_gene_exp_union=zeros(size(mouseGEM.genes,1),size(str1,2));
for i=1:size(str1,2)
    met_gene_exp_union(p1,1)=P1_Liver_clst1(1,q1)';
    met_gene_exp_union(p1,2)=P1_Liver_clst2(1,q1)';

    met_gene_exp_union(p2,3)=P1_Breast_clst1(1,q2)';
    met_gene_exp_union(p2,4)=P1_Breast_clst2(1,q2)';

    met_gene_exp_union(p3,5)=P2_Liver_clst1(1,q3)';
    met_gene_exp_union(p3,6)=P2_Liver_clst2(1,q3)';

    met_gene_exp_union(p4,7)=P2_Breast_clst1(1,q4)';
    met_gene_exp_union(p4,8)=P2_Breast_clst2(1,q4)';

    met_gene_exp_union(p5,9)=P3_Bone_clst1(1,q5)';
    met_gene_exp_union(p5,10)=P3_Bone_clst2(1,q5)';

    met_gene_exp_union(p6,11)=P3_Brain_clst1(1,q6)';
    met_gene_exp_union(p6,12)=P3_Brain_clst2(1,q6)';

    met_gene_exp_union(p7,13)=P3_Breast_clst1(1,q7)';
    met_gene_exp_union(p7,14)=P3_Breast_clst2(1,q7)';
end

% met_gene_exp_union=[['genes';string(ihuman.genes)],[str1;string(met_gene_exp_union)]];

% save met_gene_exp_union met_gene_exp_union

