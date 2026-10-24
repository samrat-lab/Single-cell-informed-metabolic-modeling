clear
clc

% Load the models
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P1_Liver_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P1_Breast_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P2_Liver_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P2_Breast_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P3_Bone_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P3_Brain_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P3_Breast_tINIT.mat')

category={'P1_Liver_tINIT.clst1', 'P1_Liver_tINIT.clst2', 'P1_Breast_tINIT.clst1', 'P1_Breast_tINIT.clst2', ...
    'P2_Liver_tINIT.clst1', 'P2_Liver_tINIT.clst2', 'P2_Breast_tINIT.clst1', 'P2_Breast_tINIT.clst2', ...
    'P3_Bone_tINIT.clst1', 'P3_Bone_tINIT.clst2', 'P3_Brain_tINIT.clst1', 'P3_Brain_tINIT.clst2', ...
    'P3_Breast_tINIT.clst1', 'P3_Breast_tINIT.clst2'};


%% Final flux values for further analysis
for i=1:size(category,2)
    
    str=[category{i},'.mic_max(1,2)'];
    tmp=eval(str);
    
    str1=[category{i},'.storedata_gpSampler.points(:,tmp)'];
    tmp1=eval(str1);
    
    str2=[category{i},'.final_flux_values', '=','tmp1'];
    eval(str2);
    
    clear 'str' 'tmp' 'str1' 'tmp1' 'str2' 'ans' 'i'
    
end

save('P1_Liver_tINIT.mat', 'P1_Liver_tINIT', '-v7.3')
save('P1_Breast_tINIT.mat', 'P1_Breast_tINIT', '-v7.3')
save('P2_Liver_tINIT.mat', 'P2_Liver_tINIT', '-v7.3')
save('P2_Breast_tINIT.mat', 'P2_Breast_tINIT', '-v7.3')
save('P3_Bone_tINIT.mat', 'P3_Bone_tINIT', '-v7.3')
save('P3_Brain_tINIT.mat', 'P3_Brain_tINIT', '-v7.3')
save('P3_Breast_tINIT.mat', 'P3_Breast_tINIT', '-v7.3')

