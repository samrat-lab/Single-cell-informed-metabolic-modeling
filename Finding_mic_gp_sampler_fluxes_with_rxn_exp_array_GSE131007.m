%%%%% Using it after the file: "Flux_values_evaluation_gp_sampler_human_GSE131007" and 
%%%%% "Evaluating_rxn_exp_array_for_each_model_GSE131007"
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


%% Finding mic of each set of flux values with the rxn_exp_array
mic=[];
for i=1:size(category,2)
    tmp=eval(category{i});
    exp=tmp.rxn_exp_array;
    
    c=tmp.storedata_gpSampler.points;
    
    k=5;
    for j=1:size(c,2)
        mic(j,1)=mi_cont_cont(exp, c(:,j), k);
        j
    end
    
    str1=[category{i},'.mic_gp_samp_fluxes_with_exp', '=','mic'];
    eval(str1);
     
    clear 'tmp' 'exp' 'str' 'c' 'mic' 'str1' 'ans'
    i
end


% %% Finding distance correlation of mean, median with the rxn_exp_array
% for i=1:size(category,2)
%     tmp=eval(category{i});
%     exp=tmp.rxn_exp_array;
%     tmp1=eval(category1{i});
%     
%     mean=tmp1.mean_flux_gp_samp;
%     k=5;
%     mic=mi_cont_cont(exp, mean, k);
%     str1=[category1{i},'.mic_mean', '=','mic'];
%     eval(str1);
%     
%     clear mean k mic str1 ans
%     
%     median=tmp1.median_flux_gp_samp;
%     k=5;
%     mic=mi_cont_cont(exp, median, k);
%     str1=[category1{i},'.mic_median', '=','mic'];
%     eval(str1);
%     
%     clear tmp tmp1 k median mic str1 ans
%     
% end


%% Finding the set of flux values with maximum mic with the rxn_exp_array
for i=1:size(category,2)
    
    str=[category{i},'.mic_gp_samp_fluxes_with_exp(:,1)'];
    tmp=eval(str);
    [a,b]=sort(tmp(:,1),'descend');
    clear str tmp
    
        
    str1=[category{i},'.mic_max', '=','[a(1,1) b(1,1)]'];
    eval(str1);
    
    clear 'a' 'b' 'a1' 'b1' 'str1'
    
end

save('P1_Liver_tINIT.mat', 'P1_Liver_tINIT', '-v7.3')
save('P1_Breast_tINIT.mat', 'P1_Breast_tINIT', '-v7.3')
save('P2_Liver_tINIT.mat', 'P2_Liver_tINIT', '-v7.3')
save('P2_Breast_tINIT.mat', 'P2_Breast_tINIT', '-v7.3')
save('P3_Bone_tINIT.mat', 'P3_Bone_tINIT', '-v7.3')
save('P3_Brain_tINIT.mat', 'P3_Brain_tINIT', '-v7.3')
save('P3_Breast_tINIT.mat', 'P3_Breast_tINIT', '-v7.3')

