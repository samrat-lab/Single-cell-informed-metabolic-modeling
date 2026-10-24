%% Capturing changes occurring at metastatic level of each organotropism compared to its primary level
clear
clc

% Load the models
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P1_Breast_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P1_Liver_tINIT.mat')


%% Changes between primary and metastatic level of patient 1 (UCD46)
category1={'P1_Breast_tINIT.clst1', 'P1_Breast_tINIT.clst2'}; % Primary breast
category2={'P1_Liver_tINIT.clst1', 'P1_Liver_tINIT.clst2'}; % Liver organotropism

Rxns_P1B_P1L={};
Flux_foldchange={};
pert_rxns={};
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
        fxc2=flux1(c2); % Primary Breast
        fxc2(fxc2<1e-6 & fxc2>0)=1e-6;
        fxc2(fxc2>-1e-6 & fxc2<0)=-1e-6;
        
        fxc3=flux2(c3); % Liver Organotropism
        fxc3(fxc3<1e-6 & fxc3>0)=1e-6;
        fxc3(fxc3>-1e-6 & fxc3<0)=-1e-6;
        
        %%%% Flux fold change %%%%
        fold_change=fxc3./fxc2;
        Flux_foldchange{i,j}=fold_change;
        
        Rxns_P1B_P1L{i,j}(:,1)=string(c1);
        Rxns_P1B_P1L{i,j}(:,2)=num2cell(fold_change);
        
        %%%% Reaction directionality %%%%
        same_sign_idx=((fxc2>=0 & fxc3>=0) | (fxc2<=0 & fxc3<=0)) &...
            ~(fxc2==0 & fxc3==0); % Same sign idx
        
        forward_idx=fxc2>0 & fxc3>0; % Forward direction
        Rxns_P1B_P1L{i,j}(find(forward_idx),3)={'Forward_reaction'};
        Rxns_P1B_P1L{i,j}(find(forward_idx),4)={'Same_sign_reaction'};
        
        backward_idx=fxc2<0 & fxc3<0; % Backward direction
        Rxns_P1B_P1L{i,j}(find(backward_idx),3)={'Backward_reaction'};
        Rxns_P1B_P1L{i,j}(find(backward_idx),4)={'Same_sign_reaction'};
        
        activation_idx=fxc2==0 & fxc3~=0; % Activation
        Rxns_P1B_P1L{i,j}(find(activation_idx),3)={'Activation'};
        Rxns_P1B_P1L{i,j}(find(activation_idx),4)={'Same_sign_reaction'};
        
        inactivation_idx=fxc2~=0 & fxc3==0; % Inactivation
        Rxns_P1B_P1L{i,j}(find(inactivation_idx),3)={'Inactivation'};
        Rxns_P1B_P1L{i,j}(find(inactivation_idx),4)={'Same_sign_reaction'};
        
        opp_sign_idx=(fxc2<0 & fxc3>0) | (fxc2>0 & fxc3<0); % Direction switch
        Rxns_P1B_P1L{i,j}(find(opp_sign_idx),3)={'Direction_switch'};
        Rxns_P1B_P1L{i,j}(find(opp_sign_idx),4)={'Opposite_sign_reaction'};
        
        
        pert_rxns{i,j}=c1(or(abs(fxc3./fxc2)>=1.2,abs(fxc3./fxc2)<=(1/(1.2)))); % Complete Perturbed
        
        
        clear str1 str2 str3 str4 c1 c2 c3 rxns1 rxns2 flux1 flux2 fxc2 fxc3 fold_change same_sign_idx...
            forward_idx backward_idx activation_idx inactivation_idx opp_sign_idx
        
    end
end


%% Filtering perturbed reactions
Rxns_P1B_P1L_opp_sign={}; % Not perturbed but direction changed
Pert_Rxns_P1B_P1L={}; % Total (same_sign + opp_sign)
for i=1:size(Rxns_P1B_P1L,2)
    for j=1:size(Rxns_P1B_P1L,2)
        
        temp=Rxns_P1B_P1L{i,j};
        
        %%%% Perturbed and same sign
        
        % Upregulated & same: fold change > 2 and same sign
        up_idx=double(temp(:,2))>=(1.2) & temp(:,4)=='Same_sign_reaction';
        temp(up_idx,5)='Up-regulated';
        
        % Downregulated & same: fold change < 0.5 and same sign
        down_idx=double(temp(:,2))<=(1/(1.2)) & temp(:,4)=='Same_sign_reaction';
        temp(down_idx,5)='Down-regulated';
        
        % Pert_Rxns_P1_P3_same_sign{i,j}=temp(up_idx | down_idx,:);
        
        %%%% Perturbed and opposite sign
        
        % Upregulated & opp: fold change > 1.2 and opposite sign
        up_idx1=abs(double(temp(:,2)))>=(1.2) & temp(:,4)=='Opposite_sign_reaction';
        temp(up_idx1,5)='Up-regulated';
        
        % Downregulated & opp: fold change < (1/(2)) and opposite sign
        down_idx1=abs(double(temp(:,2)))<=(1/(1.2)) & temp(:,4)=='Opposite_sign_reaction';
        temp(down_idx1,5)='Down-regulated';
        
        % Pert_Rxns_P1_P3_opp_sign{i,j}=temp(up_idx1 | down_idx1,:);
        
        Pert_Rxns_P1B_P1L{i,j}=temp(up_idx | down_idx | up_idx1 | down_idx1, :);
        
        idx=double(temp(:,2))==-1 & temp(:,4)=='Opposite_sign_reaction';
        Rxns_P1B_P1L_opp_sign{i,j}=temp(idx);
        
    end
    clear temp up_idx up_idx1 down_idx down_idx1 idx
end


%% Taking union of all reactions in each case (Substitute to above)
filt_pert_rxns_P1B_P1L=union(union(union(Pert_Rxns_P1B_P1L{1,1}(:,1),Pert_Rxns_P1B_P1L{1,2}(:,1)),...
    Pert_Rxns_P1B_P1L{2,1}(:,1)),Pert_Rxns_P1B_P1L{2,2}(:,1));

% Fold change in all four cases
k=0;
Tot_Pert_Rxns_P1B_P1L=zeros(size(filt_pert_rxns_P1B_P1L,1),4);
for i=1:2
    for j=1:2
        k=k+1;
        [p1,p2,p3]=intersect(string(Pert_Rxns_P1B_P1L{i,j}(:,1)),filt_pert_rxns_P1B_P1L);
        Tot_Pert_Rxns_P1B_P1L(p3,k)=Pert_Rxns_P1B_P1L{i,j}(p2,2);
    end
    clear p1 p2 p3
end
Tot_Pert_Rxns_P1B_P1L=[filt_pert_rxns_P1B_P1L,Tot_Pert_Rxns_P1B_P1L];

% save Tot_Pert_Rxns_P1B_P1L Tot_Pert_Rxns_P1B_P1L % saved in Metastatic_level_changes_GSE131007


% Reaction direction information in all four cases
k=0;
Direction_Tot_Pert_Rxns_P1B_P1L=strings();
for i=1:2
    for j=1:2
        k=k+1;
        [p1,p2,p3]=intersect(string(Pert_Rxns_P1B_P1L{i,j}(:,1)),filt_pert_rxns_P1B_P1L);
        Direction_Tot_Pert_Rxns_P1B_P1L(p3,k)=Pert_Rxns_P1B_P1L{i,j}(p2,4);
    end
    clear p1 p2 p3
end
Direction_Tot_Pert_Rxns_P1B_P1L(find(ismissing(Direction_Tot_Pert_Rxns_P1B_P1L)))='Not_present';
Direction_Tot_Pert_Rxns_P1B_P1L=[filt_pert_rxns_P1B_P1L,Direction_Tot_Pert_Rxns_P1B_P1L];

% save Direction_Tot_Pert_Rxns_P1B_P1L Direction_Tot_Pert_Rxns_P1B_P1L % saved in Metastatic_level_changes_GSE131007


% Reaction nature information in all four cases
k=0;
Nature_Tot_Pert_Rxns_P1B_P1L=strings();
for i=1:2
    for j=1:2
        k=k+1;
        [p1,p2,p3]=intersect(string(Pert_Rxns_P1B_P1L{i,j}(:,1)),filt_pert_rxns_P1B_P1L);
        Nature_Tot_Pert_Rxns_P1B_P1L(p3,k)=Pert_Rxns_P1B_P1L{i,j}(p2,5);
    end
    clear p1 p2 p3
end
Nature_Tot_Pert_Rxns_P1B_P1L(find(ismissing(Nature_Tot_Pert_Rxns_P1B_P1L)))='Not_present';
Nature_Tot_Pert_Rxns_P1B_P1L=[filt_pert_rxns_P1B_P1L,Nature_Tot_Pert_Rxns_P1B_P1L];

% save Nature_Tot_Pert_Rxns_P1B_P1L Nature_Tot_Pert_Rxns_P1B_P1L % saved in Metastatic_level_changes_GSE131007


%% Common between P1B_P1L_pert_rxns_Breast and P2B_P2L_pert_rxns_Breast (Taking intersection)
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Metastatic_level_changes_GSE131007\Tot_Pert_Rxns_P1B_P1L.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Metastatic_level_changes_GSE131007\Direction_Tot_Pert_Rxns_P1B_P1L.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Metastatic_level_changes_GSE131007\Nature_Tot_Pert_Rxns_P1B_P1L.mat')

load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Metastatic_level_changes_GSE131007\Tot_Pert_Rxns_P2B_P2L.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Metastatic_level_changes_GSE131007\Direction_Tot_Pert_Rxns_P2B_P2L.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Metastatic_level_changes_GSE131007\Nature_Tot_Pert_Rxns_P2B_P2L.mat')

[p1,p2,p3]=intersect(Tot_Pert_Rxns_P1B_P1L(:,1), Tot_Pert_Rxns_P2B_P2L(:,1));
Metastatic_level_pert_rxns_P1_P2(:,1)=p1;
Metastatic_level_pert_rxns_P1_P2(:,2:5)=Tot_Pert_Rxns_P1B_P1L(p2,2:end);
Metastatic_level_pert_rxns_P1_P2(:,6:9)=Tot_Pert_Rxns_P2B_P2L(p3,2:end);

Metastatic_level_pert_rxns_direction_P1_P2(:,1)=p1;
Metastatic_level_pert_rxns_direction_P1_P2(:,2:5)=Direction_Tot_Pert_Rxns_P1B_P1L(p2,2:end);
Metastatic_level_pert_rxns_direction_P1_P2(:,6:9)=Direction_Tot_Pert_Rxns_P2B_P2L(p3,2:end);

Metastatic_level_pert_rxns_nature_P1_P2(:,1)=p1;
Metastatic_level_pert_rxns_nature_P1_P2(:,2:5)=Nature_Tot_Pert_Rxns_P1B_P1L(p2,2:end);
Metastatic_level_pert_rxns_nature_P1_P2(:,6:9)=Nature_Tot_Pert_Rxns_P2B_P2L(p3,2:end);

% save Metastatic_level_pert_rxns_P1_P2 Metastatic_level_pert_rxns_P1_P2
% save Metastatic_level_pert_rxns_direction_P1_P2 Metastatic_level_pert_rxns_direction_P1_P2
% save Metastatic_level_pert_rxns_nature_P1_P2 Metastatic_level_pert_rxns_nature_P1_P2

nat={};
for i=1:size(Metastatic_level_pert_rxns_P1_P2,1)
    str1=Metastatic_level_pert_rxns_direction_P1_P2(i,2:end);
    str2=Metastatic_level_pert_rxns_nature_P1_P2(i,2:end);
    for j=1:size(str1,2)
        g1=string(extractBefore(str1(1,j), "_"));
        g2=string(extractBefore(str2(1,j), "-"));
        A=[g1, g2];
        nat{i,j}=strjoin(A, "_");
    end
    clear g1 g2 A
end
ans1=find(ismissing(string(nat)));
nat(ans1)={"Not_present"};
nat=string(nat);

% % Checking
% ans1=find(nat=='Not_present'); % 589 positions
% ans2=find(Metastatic_level_pert_rxns_direction_P1_P2(:,2:end)=='Not_present'); % 589 positions
% ans3=find(Metastatic_level_pert_rxns_nature_P1_P2(:,2:end)=='Not_present'); % 589 positions
% intersect(ans1,ans2) % same 589
% intersect(ans1,ans3) % same 589
% intersect(ans2,ans3) % same 589


%% Finding the number of reactions in each case above
idx1={};
idx2={};
idx3={};
idx4={};
idx5={};
l1=[];
l2=[];
l3=[];
l4=[];
l5=[];

for i=1:size(nat,1)
    
    idx1{i,1}=find(nat(i,:) == 'Same_Up');
    idx2{i,1}=find(nat(i,:) == 'Same_Down');
    idx3{i,1}=find(nat(i,:) == 'Opposite_Up');
    idx4{i,1}=find(nat(i,:) == 'Opposite_Down');
    idx5{i,1}=find(nat(i,:) == 'Not_present');
    
    l1(i,1)=length(idx1{i,1});
    l2(i,1)=length(idx2{i,1});
    l3(i,1)=length(idx3{i,1});
    l4(i,1)=length(idx4{i,1});
    l5(i,1)=length(idx5{i,1});
    
end

% % Checking
% find((l1+l2+l3+l4+l5)~=8) % No reaction

%% Deciding direction/nature of the reactions across 8 cluster cases
Reaction_class=strings(size(nat,1),2);
for i=1:size(nat,1)
    
    values=[l1(i,1), l2(i,1), l3(i,1), l4(i,1), l5(i,1)];
    [max_val, idx]=max(values);
    
    if idx == 1
        Reaction_class(i,1)="Same_Up=" + max_val; % l1 represent Same_Up
        Reaction_class(i,2)="Same_Up";
        
    elseif idx == 2
        Reaction_class(i,1)="Same_Down=" + max_val; % l2 represent Same_Down
        Reaction_class(i,2)="Same_Down";
        
    elseif idx == 3
        Reaction_class(i,1)="Opposite_Up=" + max_val; % l3 represent Opposite_Up
        Reaction_class(i,2)="Opposite_Up";
        
    elseif idx == 4
        Reaction_class(i,1)="Opposite_Down=" + max_val; % l4 represent Opposite_Down
        Reaction_class(i,2)="Opposite_Down";
              
    else % idx == 5 --> Not_present is max % l5 represent Not_present
        % Check if Same_Up/Same_Down/Opposite_Up/Opposite_Down exists (non-zero), we picked the larger of them
        if l1(i,1) > 0 || l2(i,1) > 0 || l3(i,1) > 0 || l4(i,1) > 0
            if l1(i,1) == max([l1(i,1), l2(i,1), l3(i,1), l4(i,1)])
                Reaction_class(i,1)="Same_Up=" + l1(i,1);
                Reaction_class(i,2)="Same_Up";
                
            elseif l2(i,1) == max([l1(i,1), l2(i,1), l3(i,1), l4(i,1)])
                Reaction_class(i,1)="Same_Down=" + l2(i,1);
                Reaction_class(i,2)="Same_Down";
                
            elseif l3(i,1) == max([l1(i,1), l2(i,1), l3(i,1), l4(i,1)])
                Reaction_class(i,1)="Opposite_Up=" + l3(i,1);
                Reaction_class(i,2)="Opposite_Up";
                
            elseif l4(i,1) == max([l1(i,1), l2(i,1), l3(i,1), l4(i,1)])
                Reaction_class(i,1)="Opposite_Down=" + l4(i,1);
                Reaction_class(i,2)="Opposite_Down";
            end
        else
            % Only if Same_Up, Same_Down, Opposite_Up and Opposite_Down are all zero, keep Not_present
            Reaction_class(i,1)="Not_present=" + l5(i,1);
            Reaction_class(i,2)="Not_present";
        end
    end
    
    clear values max_val idx
end

find(Reaction_class(:,2)=="Same_Up") % 185 reactions
find(Reaction_class(:,2)=="Same_Down") % 159 reactions
find(Reaction_class(:,2)=="Opposite_Up"); % 29 reactions
find(Reaction_class(:,2)=="Opposite_Down"); % 23 reactions
find(Reaction_class(:,2)=="Not_present"); % 0 reactions

% save Reaction_class Reaction_class


%% Finalizing the flux foldchange of these reactions
mean_fold_change=[];
g=double(Metastatic_level_pert_rxns_P1_P2(:,2:end));
for i=1:size(Reaction_class,1)
    ind=find(nat(i,:)==Reaction_class(i,2));
    mean_fold_change(i,1)=mean(g(i,ind));
    clear ind
end

Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007=[["Reaction_ID","Reaction_class","Reaction_class","Mean_fold_change"];...
    [Metastatic_level_pert_rxns_P1_P2(:,1),Reaction_class,mean_fold_change]];

% save Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007 Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007

%% Pathway analysis
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Human-GEM-main\model\Human-GEM.mat');
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\Metastatic_level_changes_GSE131007\Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007.mat')

model=ihuman;
[~,~,rxnSet]=intersect(Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007(:,1),string(model.rxns(:,1)));
group='subSystems';
pathways= FEA(model, rxnSet, group);

find(cell2mat(pathways(2:end,1))<0.05);
sig_path=[pathways(1,:); pathways(ans+1,:)];


%% Finding distribution (up, down) of perturbed reactions in each perturbed pathways
[~,~,rxnSet]=intersect(Metastatic_level_pert_rxns_P1_P2_Liver_final_GSE131007(:,1),string(model.rxns(:,1)));
rxns=[string(ihuman.rxns(rxnSet)),string(ihuman.subSystems(rxnSet))];
path=sig_path(2:end,3);

for i=1:size(path,1)
    find(rxns(:,2)==string(path{i,1}));
    path{i,2}=rxns(ans,1);
    [~,~,idx]=intersect(string(path{i,2}),Metastatic_level_pert_rxns_P1_P2(:,1));
    path{i,3}=Metastatic_level_pert_rxns_P1_P2(idx,3);
    l1=length(find(string(path{i,3})=='Up-regulated'));
    l2=length(find(string(path{i,3})=='Down-regulated'));
    path{i,4}=l1;
    path{i,5}=l2;
end


