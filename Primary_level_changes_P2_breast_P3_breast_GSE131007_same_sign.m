%% Capturing changes occurring at primary level of each patient
clear
clc

% Load the models
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P2_Breast_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P3_Breast_tINIT.mat')


%% Changes between patient 2 (UCD4) and patient 3 (UCD65)
category1={'P2_Breast_tINIT.clst1', 'P2_Breast_tINIT.clst2'}; % liver organotropism
category2={'P3_Breast_tINIT.clst1', 'P3_Breast_tINIT.clst2'}; % bone and brain organotropism

Rxns_P2_P3={};
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
        
        Rxns_P2_P3{i,j}(:,1)=string(c1);
        Rxns_P2_P3{i,j}(:,2)=num2cell(fold_change);
        
        
        %%%% Reaction directionality %%%%
        same_sign_idx=((fxc2>=0 & fxc3>=0) | (fxc2<=0 & fxc3<=0)) &...
            ~(fxc2==0 & fxc3==0); % Same sign idx
        
        forward_idx=fxc2>0 & fxc3>0; % Forward direction
        Rxns_P2_P3{i,j}(find(forward_idx),3)={'Forward_reaction'};
        Rxns_P2_P3{i,j}(find(forward_idx),4)={'Same_sign_reaction'};
        
        backward_idx=fxc2<0 & fxc3<0; % Backward direction
        Rxns_P2_P3{i,j}(find(backward_idx),3)={'Backward_reaction'};
        Rxns_P2_P3{i,j}(find(backward_idx),4)={'Same_sign_reaction'};
        
        activation_idx=fxc2==0 & fxc3~=0; % Activation
        Rxns_P2_P3{i,j}(find(activation_idx),3)={'Activation'};
        Rxns_P2_P3{i,j}(find(activation_idx),4)={'Same_sign_reaction'};
        
        inactivation_idx=fxc2~=0 & fxc3==0; % Inactivation
        Rxns_P2_P3{i,j}(find(inactivation_idx),3)={'Inactivation'};
        Rxns_P2_P3{i,j}(find(inactivation_idx),4)={'Same_sign_reaction'};
        
        opp_sign_idx=(fxc2<0 & fxc3>0) | (fxc2>0 & fxc3<0); % Direction switch
        Rxns_P2_P3{i,j}(find(opp_sign_idx),3)={'Direction_switch'};
        Rxns_P2_P3{i,j}(find(opp_sign_idx),4)={'Opposite_sign_reaction'};
        
        
        pert_rxns{i,j}=c1(or(abs(fxc3./fxc2)>=1.2,abs(fxc3./fxc2)<=(1/(1.2)))); % Complete Perturbed (Bone and Brain vs Liver)
      
        clear str1 str2 str3 str4 c1 c2 c3 rxns1 rxns2 flux1 flux2 fxc2 fxc3 fold_change same_sign_idx...
            forward_idx backward_idx activation_idx inactivation_idx opp_sign_idx
        
    end
end


%% Filtering perturbed reactions
Rxns_P2_P3_opp_sign={}; % Not perturbed but direction changed
Pert_Rxns_P2_P3={}; % Total (same_sign + opp_sign)
for i=1:size(Rxns_P2_P3,2)
    for j=1:size(Rxns_P2_P3,2)
        
        temp=Rxns_P2_P3{i,j};
        
        %%%% Perturbed and same sign
        
        % Upregulated & same: fold change > 2 and same sign
        up_idx=double(temp(:,2))>=(1.2) & temp(:,4)=='Same_sign_reaction';
        temp(up_idx,5)='Up-regulated';
        
        % Downregulated & same: fold change < 0.5 and same sign
        down_idx=double(temp(:,2))<=(1/(1.2)) & temp(:,4)=='Same_sign_reaction';
        temp(down_idx,5)='Down-regulated';
        
        % Pert_Rxns_P2_P3_same_sign{i,j}=temp(up_idx | down_idx,:);
        
        %%%% Perturbed and opposite sign
        
        % Upregulated & opp: fold change > 1.2 and opposite sign
        up_idx1=abs(double(temp(:,2)))>=(1.2) & temp(:,4)=='Opposite_sign_reaction';
        temp(up_idx1,5)='Up-regulated';
        
        % Downregulated & opp: fold change < (1/(2)) and opposite sign
        down_idx1=abs(double(temp(:,2)))<=(1/(1.2)) & temp(:,4)=='Opposite_sign_reaction';
        temp(down_idx1,5)='Down-regulated';
        
        % Pert_Rxns_P1_P3_opp_sign{i,j}=temp(up_idx1 | down_idx1,:);
        
        Pert_Rxns_P2_P3{i,j}=temp(up_idx | down_idx | up_idx1 | down_idx1, :);
        
        idx=double(temp(:,2))==-1 & temp(:,4)=='Opposite_sign_reaction';
        Rxns_P2_P3_opp_sign{i,j}=temp(idx);
        
    end
    clear temp up_idx up_idx1 down_idx down_idx1
end


%% Taking union of all reactions in each case (Substitute to above)
filt_pert_rxns_P2_P3=union(union(union(Pert_Rxns_P2_P3{1,1}(:,1),Pert_Rxns_P2_P3{1,2}(:,1)),...
    Pert_Rxns_P2_P3{2,1}(:,1)),Pert_Rxns_P2_P3{2,2}(:,1));

% Fold change in all four cases
k=0;
Tot_Pert_Rxns_P2_P3=zeros(size(filt_pert_rxns_P2_P3,1),4);
for i=1:2
    for j=1:2
            k=k+1;
            [p1,p2,p3]=intersect(string(Pert_Rxns_P2_P3{i,j}(:,1)),filt_pert_rxns_P2_P3);
            Tot_Pert_Rxns_P2_P3(p3,k)=Pert_Rxns_P2_P3{i,j}(p2,2);
    end
    clear p1 p2 p3
end
Tot_Pert_Rxns_P2_P3=[filt_pert_rxns_P2_P3,Tot_Pert_Rxns_P2_P3];

% save Tot_Pert_Rxns_P2_P3 Tot_Pert_Rxns_P2_P3 % saved in Primary_breast_changes


% Reaction direction information in all four cases
k=0;
Direction_Tot_Pert_Rxns_P2_P3=strings();
for i=1:2
    for j=1:2
            k=k+1;
            [p1,p2,p3]=intersect(string(Pert_Rxns_P2_P3{i,j}(:,1)),filt_pert_rxns_P2_P3);
            Direction_Tot_Pert_Rxns_P2_P3(p3,k)=Pert_Rxns_P2_P3{i,j}(p2,4);
    end
    clear p1 p2 p3
end
Direction_Tot_Pert_Rxns_P2_P3(find(ismissing(Direction_Tot_Pert_Rxns_P2_P3)))='Not_present';
Direction_Tot_Pert_Rxns_P2_P3=[filt_pert_rxns_P2_P3,Direction_Tot_Pert_Rxns_P2_P3];

% save Direction_Tot_Pert_Rxns_P2_P3 Direction_Tot_Pert_Rxns_P2_P3 % saved in Primary_breast_changes


% Reaction nature information in all four cases
k=0;
Nature_Tot_Pert_Rxns_P2_P3=strings();
for i=1:2
    for j=1:2
        k=k+1;
        [p1,p2,p3]=intersect(string(Pert_Rxns_P2_P3{i,j}(:,1)),filt_pert_rxns_P2_P3);
        Nature_Tot_Pert_Rxns_P2_P3(p3,k)=Pert_Rxns_P2_P3{i,j}(p2,5);
    end
    clear p1 p2 p3
end
Nature_Tot_Pert_Rxns_P2_P3(find(ismissing(Nature_Tot_Pert_Rxns_P2_P3)))='Not_present';
Nature_Tot_Pert_Rxns_P2_P3=[filt_pert_rxns_P2_P3,Nature_Tot_Pert_Rxns_P2_P3];

% save Nature_Tot_Pert_Rxns_P2_P3 Nature_Tot_Pert_Rxns_P2_P3 % saved in Primary_breast_changes


%% Common between P1_P3_pert_rxns_Breast and P2_P3_pert_rxns_Breast
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Primary_breast_changes_GSE131007\Tot_Pert_Rxns_P1_P3.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Primary_breast_changes_GSE131007\Direction_Tot_Pert_Rxns_P1_P3.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Primary_breast_changes_GSE131007\Nature_Tot_Pert_Rxns_P1_P3.mat')

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Primary_breast_changes_GSE131007\Tot_Pert_Rxns_P2_P3.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Primary_breast_changes_GSE131007\Direction_Tot_Pert_Rxns_P2_P3.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Primary_breast_changes_GSE131007\Nature_Tot_Pert_Rxns_P2_P3.mat')

[p1,p2,p3]=intersect(Tot_Pert_Rxns_P1_P3(:,1), Tot_Pert_Rxns_P2_P3(:,1));
Primary_level_pert_rxns(:,1)=p1;
Primary_level_pert_rxns(:,2:5)=Tot_Pert_Rxns_P1_P3(p2,2:end);
Primary_level_pert_rxns(:,6:9)=Tot_Pert_Rxns_P2_P3(p3,2:end);

Primary_level_pert_rxns_direction(:,1)=p1;
Primary_level_pert_rxns_direction(:,2:5)=Direction_Tot_Pert_Rxns_P1_P3(p2,2:end);
Primary_level_pert_rxns_direction(:,6:9)=Direction_Tot_Pert_Rxns_P2_P3(p3,2:end);

Primary_level_pert_rxns_nature(:,1)=p1;
Primary_level_pert_rxns_nature(:,2:5)=Nature_Tot_Pert_Rxns_P1_P3(p2,2:end);
Primary_level_pert_rxns_nature(:,6:9)=Nature_Tot_Pert_Rxns_P2_P3(p3,2:end);


nat = strings(size(p1,1),2);
for i=1:size(p1,1)
    values=[l(i,1), l1(i,1), l2(i,1)];
    [max_val, idx]=max(values);
    
    if idx == 1
        nat(i,1)="Up-regulated=" + max_val; % l represent Up-regulated
        nat(i,2)="Up-regulated";
        
    elseif idx == 2
        nat(i,1)="Down-regulated=" + max_val; % l1 represent Down-regulated
        nat(i,2)="Down-regulated";
        
    else % idx == 3 --> Not_present is max
        % Check if Up or Down exists (non-zero), if so pick the larger of Up/Down
        if l(i,1) > 0 || l1(i,1) > 0
            if l(i,1) >= l1(i,1)
                nat(i,1)="Up-regulated=" + l(i,1);
                nat(i,2)="Up-regulated";
            else
                nat(i,1)="Down-regulated=" + l1(i,1);
                nat(i,2)="Down-regulated";
            end
        else
            % Only if Up and Down are both zero, keep Not_present
            nat(i,1)="Not_present=" + l2(i,1);
            nat(i,2)="Not_present";
        end
    end
    
    clear values max_val idx
end

find(nat(:,2)=="Up-regulated") % 191 reactions
find(nat(:,2)=="Down-regulated") % 213 reactions
nat(find(nat(:,2)=="Not_present"),:) % 0 reactions

%%%% However, we observed 43 reactions, the number of upregulated == 4 and number of downregulated ==4
idx=find(l==4 & l1==4); % Difficult to give them a particular reaction nature.
nat(idx,:)=[];




mean_fold_change=[];
final_rxns=strings();
nat=strings();
n=1;
for i=1:size(Primary_level_pert_rxns_direction,1)
    idx=find(Primary_level_pert_rxns_direction(i,2:end) == 'Up-regulated'); 
    idx1=find(Primary_level_pert_rxns_direction(i,2:end) == 'Down-regulated');
    
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

Primary_level_pert_rxns_final=[final_rxns, mean_fold_change, nat];

save Primary_level_pert_rxns_final Primary_level_pert_rxns_final


%% Finding how many are up-regulated and down-regulated
length(find(Primary_level_pert_rxns_final(:,3)=='Up-regulated')); % 75 reactions are Up

length(find(Primary_level_pert_rxns_final(:,3)=='Down-regulated')); % 161 reactions are Down


%% Pathway analysis
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Human-GEM-main\model\Human-GEM.mat');
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Primary_breast_changes_GSE131007\Primary_level_pert_rxns_final.mat')

model=ihuman;
[~,~,rxnSet]=intersect(Primary_level_pert_rxns_final(:,1),string(model.rxns(:,1)));
group='subSystems';
pathways= FEA(model, rxnSet, group);

find(cell2mat(pathways(2:end,1))<0.05);
sig_path=[pathways(1,:); pathways(ans+1,:)];


%% Finding distribution (up, down) of perturbed reactions in each perturbed pathways
[~,~,rxnSet]=intersect(Primary_level_pert_rxns_final(:,1),string(model.rxns(:,1)));
rxns=[string(ihuman.rxns(rxnSet)),string(ihuman.subSystems(rxnSet))];
path=sig_path(2:end,3);

for i=1:size(path,1)
    find(rxns(:,2)==string(path{i,1}));
    path{i,2}=rxns(ans,1);
    [~,~,idx]=intersect(string(path{i,2}),Primary_level_pert_rxns_final(:,1));
    path{i,3}=Primary_level_pert_rxns_final(idx,3);
    l1=length(find(string(path{i,3})=='Up-regulated'));
    l2=length(find(string(path{i,3})=='Down-regulated'));
    path{i,4}=l1;
    path{i,5}=l2;
end






