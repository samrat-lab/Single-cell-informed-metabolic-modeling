clear
clc

% Load the models
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P1_Breast_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P3_Breast_tINIT.mat')


%% Changes between patient 1 (UCD46) and patient 3 (UCD65)
category1={'P1_Breast_tINIT.clst1', 'P1_Breast_tINIT.clst2'}; % liver organotropism
category2={'P3_Breast_tINIT.clst1', 'P3_Breast_tINIT.clst2'}; % bone and brain organotropism

Rxns_P1_P3={};
Flux_foldchange={};
pert_rxns2={};
for i=1:size(category1,2)
    
    for j=1:size(category2,2)
        
        str1=[category1{i},'.model.rxns'];
        rxns1=eval(str1);
        
        str2=[category2{j},'.model.rxns'];
        rxns2=eval(str2);
        
        [c1,c2,c3]=intersect(string(rxns1),string(rxns2));
        
        str3=[category1{i},'.final_flux_values'];
        flux1=eval(str3);
        
        str4=[category2{j},'.final_flux_values'];
        flux2=eval(str4);
        
        %%%% Filtering very small fluxes that arise from solver tolerance or floating-point precision %%%%
        fxc2=flux1(c2);
        fxc2(fxc2<1e-6 & fxc2>0)=1e-6;
        fxc2(fxc2>-1e-6 & fxc2<0)=-1e-6;
        
        fxc3=flux2(c3);
        fxc3(fxc3<1e-6 & fxc3>0)=1e-6;
        fxc3(fxc3>-1e-6 & fxc3<0)=-1e-6;
        
        %%%% Flux fold change %%%%
        fold_change=fxc3./fxc2;
        Flux_foldchange{i,j}=fold_change;
        
        Rxns_P1_P3{i,j}(:,1)=string(c1);
        Rxns_P1_P3{i,j}(:,2)=num2cell(fold_change);
        
        %%%% Reaction directionality %%%%
        same_sign_idx=((fxc2>=0 & fxc3>=0) | (fxc2<=0 & fxc3<=0)) &...
            ~(fxc2==0 & fxc3==0); % Same sign idx
        
        forward_idx=fxc2>0 & fxc3>0; % Forward direction
        Rxns_P1_P3{i,j}(find(forward_idx),3)={'Forward_reaction'};
        Rxns_P1_P3{i,j}(find(forward_idx),4)={'Same_sign_reaction'};
        
        backward_idx=fxc2<0 & fxc3<0; % Backward direction
        Rxns_P1_P3{i,j}(find(backward_idx),3)={'Backward_reaction'};
        Rxns_P1_P3{i,j}(find(backward_idx),4)={'Same_sign_reaction'};
        
        activation_idx=fxc2==0 & fxc3~=0; % Activation
        Rxns_P1_P3{i,j}(find(activation_idx),3)={'Activation'};
        Rxns_P1_P3{i,j}(find(activation_idx),4)={'Same_sign_reaction'};
        
        inactivation_idx=fxc2~=0 & fxc3==0; % Inactivation
        Rxns_P1_P3{i,j}(find(inactivation_idx),3)={'Inactivation'};
        Rxns_P1_P3{i,j}(find(inactivation_idx),4)={'Same_sign_reaction'};
        
        opp_sign_idx=(fxc2<0 & fxc3>0) | (fxc2>0 & fxc3<0); % Direction switch
        Rxns_P1_P3{i,j}(find(opp_sign_idx),3)={'Direction_switch'};
        Rxns_P1_P3{i,j}(find(opp_sign_idx),4)={'Opposite_sign_reaction'};
        
        
        pert_rxns2{i,j}=c1(or(abs(fxc3./fxc2)>=2,abs(fxc3./fxc2)<=0.5)); % Bone and Brain vs Liver

        
        clear str1 str2 str3 str4 c1 c2 c3 rxns1 rxns2 flux1 flux2 fxc2 fxc3 fold_change same_sign_idx...
            forward_idx backward_idx activation_idx inactivation_idx opp_sign_idx
        
    end
end

P1_P3_pert_rxns_Breast_same_sign=union(union(union(pert_rxns2{1,1},pert_rxns2{1,2}),...
    pert_rxns2{2,1}),pert_rxns2{2,2});

P1_P3_pert_rxns_Breast_=intersect(intersect(intersect(Rxns_P1_P3{1,1},Rxns_P1_P3{1,2}),...
    Rxns_P1_P3{2,1}),Rxns_P1_P3{2,2});


P1_P3_pert_rxns_Breast_same_sign_int=intersect(intersect(intersect(pert_rxns2{1,1},pert_rxns2{1,2}),...
    pert_rxns2{2,1}),pert_rxns2{2,2});