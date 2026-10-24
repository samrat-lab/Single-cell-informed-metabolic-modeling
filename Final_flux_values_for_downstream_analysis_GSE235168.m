clear
clc

% Load the models
% load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE235168_variables\Primary_breast_ftINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE235168_variables\Primary_breast_tINIT.mat')


%% Final flux values for further analysis
% category={'Primary_breast_ftINIT.clst1', 'Primary_breast_ftINIT.clst2'};
category={'Primary_breast_tINIT.clst1', 'Primary_breast_tINIT.clst2'};

for i=1:size(category,2)
    
    str=[category{i},'.mic_max(1,2)'];
    tmp=eval(str);
    
    str1=[category{i},'.storedata_gpSampler.points(:,tmp)'];
    tmp1=eval(str1);
    
    str2=[category{i},'.final_flux_values', '=','tmp1'];
    eval(str2);
    
    clear 'str' 'tmp' 'str1' 'tmp1' 'str2' 'ans' 'i'
    
end

% save('Primary_breast_ftINIT.mat', 'Primary_breast_ftINIT', '-v7.3')
save('Primary_breast_tINIT.mat', 'Primary_breast_tINIT', '-v7.3')



