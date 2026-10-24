%% Capturing changes occurring at primary level of each patient
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
        same_sign_idx=((flux1(c2)>=0 & flux2(c3)>=0) | (flux1(c2)<=0 & flux2(c3)<=0)) &...
            ~(flux1(c2)==0 & flux2(c3)==0); % Same sign idx
        
        forward_idx=flux1(c2)>0 & flux2(c3)>0; % Forward direction
        Rxns_P1_P3{i,j}(find(forward_idx),3)={'Forward_reaction'};
        Rxns_P1_P3{i,j}(find(forward_idx),4)={'Same_sign_reaction'};
        
        backward_idx=flux1(c2)<0 & flux2(c3)<0; % Backward direction
        Rxns_P1_P3{i,j}(find(backward_idx),3)={'Backward_reaction'};
        Rxns_P1_P3{i,j}(find(backward_idx),4)={'Same_sign_reaction'};
        
        activation_idx=flux1(c2)==0 & flux2(c3)~=0; % Activation
        Rxns_P1_P3{i,j}(find(activation_idx),3)={'Activation'};
        Rxns_P1_P3{i,j}(find(activation_idx),4)={'Same_sign_reaction'};
        
        inactivation_idx=flux1(c2)~=0 & flux2(c3)==0; % Inactivation
        Rxns_P1_P3{i,j}(find(inactivation_idx),3)={'Inactivation'};
        Rxns_P1_P3{i,j}(find(inactivation_idx),4)={'Same_sign_reaction'};
        
        opp_sign_idx=(flux1(c2)<0 & flux2(c3)>0) | (flux1(c2)>0 & flux2(c3)<0); % Direction switch
        Rxns_P1_P3{i,j}(find(opp_sign_idx),3)={'Direction_switch'};
        Rxns_P1_P3{i,j}(find(opp_sign_idx),4)={'Opposite_sign_reaction'};
        
        clear str1 str2 str3 str4 c1 c2 c3 rxns1 rxns2 flux1 flux2 fxc2 fxc3 fold_change same_sign_idx...
            forward_idx backward_idx activation_idx inactivation_idx opp_sign_idx
        
    end
end


%% Filtering perturbed reactions
Pert_Rxns_P1_P3_same_sign={};
Pert_Rxns_P1_P3_opp_sign={};
Rxns_P1_P3_opp_sign={}; % Not perturbed but direction changed
for i=1:size(Rxns_P1_P3,2)
    for j=1:size(Rxns_P1_P3,2)
        
        temp=Rxns_P1_P3{i,j};
        
        %%%% Perturbed and same sign
        
        % Upregulated & same: fold change > 1.2 and same sign
        up_idx=double(temp(:,2))>=(1.2) & temp(:,4)=='Same_sign_reaction';
        temp(up_idx,5)='Up-regulated';
        
        % Downregulated & same: fold change < (1/(2)) and same sign
        down_idx=double(temp(:,2))<=(1/(1.2)) & temp(:,4)=='Same_sign_reaction';
        temp(down_idx,5)='Down-regulated';
        
        Pert_Rxns_P1_P3_same_sign{i,j}=temp(up_idx | down_idx,:);
        
        %%%% Perturbed and opposite sign
        
        % Upregulated & opp: fold change > 1.2 and opposite sign
        up_idx1=abs(double(temp(:,2)))>=(1.2) & temp(:,4)=='Opposite_sign_reaction';
        temp(up_idx1,5)='Up-regulated';
        
        % Downregulated & opp: fold change < (1/(2)) and opposite sign
        down_idx1=abs(double(temp(:,2)))<=(1/(1.2)) & temp(:,4)=='Opposite_sign_reaction';
        temp(down_idx1,5)='Down-regulated';
        
        Pert_Rxns_P1_P3_opp_sign{i,j}=temp(up_idx1 | down_idx1,:);
        
        idx=double(temp(:,2))==-1 & temp(:,4)=='Opposite_sign_reaction';
        Rxns_P1_P3_opp_sign{i,j}=temp(idx);
        
    end
    clear temp up_idx up_idx1 down_idx down_idx1 idx
end


%% Concatenating all the perturbed reactions obtained in each case ---- OPPOSITE SIGN
all_rxns={};
for i=1:size(Pert_Rxns_P1_P3_opp_sign,1)
    for j=1:size(Pert_Rxns_P1_P3_opp_sign,2)
        rxns=Pert_Rxns_P1_P3_opp_sign{i,j}(:,1);
        all_rxns=[all_rxns;rxns];
    end
end
unique_rxns=unique(all_rxns);

% Finding in how many cluster cases, does each perturbed reaction is present
l=[];
for i=1:size(unique_rxns,1)
    idx=find(unique_rxns(i,1)==all_rxns);
    l(i,1)=length(idx);
    clear idx
end


% Filtering reactions present in more than 50% of the cluster cases
filt_pert_rxns_P1_P3=unique_rxns(l>1); % or unique_rxns(l>=2)


% Fold change in all four cases
k=0;
Tot_Pert_Rxns_P1_P3_opp_sign=zeros(size(filt_pert_rxns_P1_P3,1),4);
for i=1:2
    for j=1:2
        k=k+1;
        [p1,p2,p3]=intersect(string(Pert_Rxns_P1_P3_opp_sign{i,j}(:,1)),filt_pert_rxns_P1_P3);
        Tot_Pert_Rxns_P1_P3_opp_sign(p3,k)=Pert_Rxns_P1_P3_opp_sign{i,j}(p2,2);
    end
    clear p1 p2 p3
end
Tot_Pert_Rxns_P1_P3_opp_sign=[filt_pert_rxns_P1_P3,Tot_Pert_Rxns_P1_P3_opp_sign];

% save Tot_Pert_Rxns_P1_P3_opp_sign Tot_Pert_Rxns_P1_P3_opp_sign % saved in Primary_breast_changes_GSE131007


% Direction information in all four cases
k=0;
Direction_Tot_Pert_Rxns_P1_P3_opp_sign=strings();
for i=1:2
    for j=1:2
        k=k+1;
        [p1,p2,p3]=intersect(string(Pert_Rxns_P1_P3_opp_sign{i,j}(:,1)),filt_pert_rxns_P1_P3);
        Direction_Tot_Pert_Rxns_P1_P3_opp_sign(p3,k)=Pert_Rxns_P1_P3_opp_sign{i,j}(p2,5);
    end
    clear p1 p2 p3
end
Direction_Tot_Pert_Rxns_P1_P3_opp_sign(Direction_Tot_Pert_Rxns_P1_P3_opp_sign=="")=missing; % Find empty strings and replace with missing
Direction_Tot_Pert_Rxns_P1_P3_opp_sign(find(ismissing(Direction_Tot_Pert_Rxns_P1_P3_opp_sign)))='Not_present';
Direction_Tot_Pert_Rxns_P1_P3_opp_sign=[filt_pert_rxns_P1_P3,Direction_Tot_Pert_Rxns_P1_P3_opp_sign];

% save Direction_Tot_Pert_Rxns_P1_P3_opp_sign Direction_Tot_Pert_Rxns_P1_P3_opp_sign % saved in Primary_breast_changes_GSE131007


%% Common between P1_P3_pert_rxns_Breast and P2_P3_pert_rxns_Breast
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Primary_breast_changes_GSE131007\Tot_Pert_Rxns_P1_P3_opp_sign.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Primary_breast_changes_GSE131007\Direction_Tot_Pert_Rxns_P1_P3_opp_sign.mat')

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Primary_breast_changes_GSE131007\Tot_Pert_Rxns_P2_P3_opp_sign.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Primary_breast_changes_GSE131007\Direction_Tot_Pert_Rxns_P2_P3_opp_sign.mat')

[p1,p2,p3]=union(Tot_Pert_Rxns_P1_P3_opp_sign(:,1), Tot_Pert_Rxns_P2_P3_opp_sign(:,1));
Primary_level_pert_rxns(:,1)=p1;
Primary_level_pert_rxns(:,2:5)=Tot_Pert_Rxns_P1_P3_opp_sign(p2,2:end);
Primary_level_pert_rxns(:,6:9)=Tot_Pert_Rxns_P2_P3_opp_sign(p3,2:end);

Primary_level_pert_rxns_direction(:,1)=p1;
Primary_level_pert_rxns_direction(:,2:5)=Direction_Tot_Pert_Rxns_P1_P3_opp_sign(p2,2:end);
Primary_level_pert_rxns_direction(:,6:9)=Direction_Tot_Pert_Rxns_P2_P3_opp_sign(p3,2:end);

mean_fold_change=[];
final_rxns=strings();
nat=strings();
n = 1;
for i=1:size(Primary_level_pert_rxns_direction, 1)
    idx=find(Primary_level_pert_rxns_direction(i, 2:end) == 'Up-regulated');
    idx1=find(Primary_level_pert_rxns_direction(i, 2:end) == 'Down-regulated');
    
    l=length(idx);
    l1=length(idx1);
    
    % Skip if both are exactly 4
    if ~(l==4 && l1==4)
        if l>=4
            final_rxns(n,1)=Primary_level_pert_rxns(i,1);
            mean_fold_change(n,1)=mean(double(Primary_level_pert_rxns(i,idx+1)));
            nat(n,1)='Up-regulated';
            n=n+1;
        elseif l1>=4
            final_rxns(n,1)=Primary_level_pert_rxns(i,1);
            mean_fold_change(n,1)=mean(double(Primary_level_pert_rxns(i,idx1+1)));
            nat(n,1)='Down-regulated';
            n=n+1;
        end
    end
end

Primary_level_pert_rxns_final_opp_sign=[final_rxns, mean_fold_change, nat];

save Primary_level_pert_rxns_final_opp_sign Primary_level_pert_rxns_final_opp_sign


%% Finding how many are up-regulated and down-regulated
length(find(Primary_level_pert_rxns_final_opp_sign(:,3)=='Up-regulated')); % 75 reactions are Up

length(find(Primary_level_pert_rxns_final_opp_sign(:,3)=='Down-regulated')); % 161 reactions are Down


%% Pathway analysis
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Human-GEM-main\model\Human-GEM.mat');
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Primary_breast_changes_GSE131007\Primary_level_pert_rxns_final_opp_sign.mat')

model=ihuman;
[~,~,rxnSet]=intersect(Primary_level_pert_rxns_final_opp_sign(:,1),string(model.rxns(:,1)));
group='subSystems';
pathways= FEA(model, rxnSet, group);

find(cell2mat(pathways(2:end,1))<0.05);
sig_path=[pathways(1,:); pathways(ans+1,:)];


%% Finding distribution (up, down) of perturbed reactions in each perturbed pathways
[~,~,rxnSet]=intersect(Primary_level_pert_rxns_final_opp_sign(:,1),string(model.rxns(:,1)));
rxns=[string(ihuman.rxns(rxnSet)),string(ihuman.subSystems(rxnSet))];
path=sig_path(2:end,3);

for i=1:size(path,1)
    find(rxns(:,2)==string(path{i,1}));
    path{i,2}=rxns(ans,1);
    [~,~,idx]=intersect(string(path{i,2}),Primary_level_pert_rxns_final_opp_sign(:,1));
    path{i,3}=Primary_level_pert_rxns_final_opp_sign(idx,3);
    l1=length(find(string(path{i,3})=='Up-regulated'));
    l2=length(find(string(path{i,3})=='Down-regulated'));
    path{i,4}=l1;
    path{i,5}=l2;
end




