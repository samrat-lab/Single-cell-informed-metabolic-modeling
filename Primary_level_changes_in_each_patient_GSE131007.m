%% Capturing changes occurring at primary level of each patient
clear
clc

% Load the models
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P1_Breast_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P2_Breast_tINIT.mat')
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism_human\GSE131007_variables\tINIT_models\P3_Breast_tINIT.mat')


%% Changes between patient 1 (UCD46) and patient 3 (UCD65)
category1={'P1_Breast_tINIT.clst1', 'P1_Breast_tINIT.clst2'}; % liver organotropism
category2={'P3_Breast_tINIT.clst1', 'P3_Breast_tINIT.clst2'}; % bone and brain organotropism

rxns={};
pert_rxns={};
pert_rxns1={};
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
        
        % same_sign_idx=flux2(c3).*flux1(c2)>=0; % True if both fluxes have the same sign
        % opp_sign_idx=flux2(c3).*flux1(c2)<=0;
        
        same_sign_idx=(flux1(c2)>=0 & flux2(c3)>=0) | (flux2(c3)<=0 & flux1(c2)<=0);
        opp_sign_idx=(flux1(c2)<0 & flux2(c3)>0) | (flux2(c3)>0 & flux1(c2)<0);
        
        rxns{i,j}=c1;
        
        fxc3=flux2(c3);
        fxc3(fxc3<1e-6 & fxc3>0)=1e-6;
        fxc3(fxc3>-1e-6 & fxc3<0)=-1e-6;

        fxc2=flux1(c2);
        fxc2(fxc2<1e-6 & fxc2>0)=1e-6;
        fxc2(fxc2>-1e-6 & fxc2<0)=-1e-6;
        
        % Fold change
%         fold_change=flux2(c3)./flux1(c2);
        fold_change=fxc3./fxc2;

        Flux_foldchange{i,j}=fold_change;
        
        
        %% Perturbed and same sign
        
        % Upregulated & same: fold change > 2 and same sign
        up_idx=fold_change>=2 & same_sign_idx;
        
        % Downregulated & same: fold change < 0.5 and same sign
        down_idx=fold_change<=0.5 & same_sign_idx;
        
        up_rxns{i,j}=c1(up_idx);
        down_rxns{i,j}=c1(down_idx);
        
        % Perturbed (up or down, same sign only)
        pert_rxns{i,j}=c1(up_idx | down_idx);
        
        
        %% Perturbed and opposite sign
        
        % Upregulated & opp: fold change > 2 and opposite sign
        up_idx1=fold_change>=2 & opp_sign_idx;
        
        % Downregulated & opp: fold change < 0.5 and opposite sign
        down_idx1=fold_change<=0.5 & opp_sign_idx;
        
        up_rxns1{i,j}=c1(up_idx1);
        down_rxns1{i,j}=c1(down_idx1);
        
        % Perturbed (up or down, opp sign only)
        pert_rxns1{i,j}=c1(up_idx1 | down_idx1);
        
        
        %% All Perturbed
        % All up reactions
        %         up_rxns2{i,j}=c1(abs(flux2(c3)./flux1(c2))>=2);
        
        % All down reactions
        %         down_rxns2{i,j}=c1(abs(flux2(c3)./flux1(c2))<=0.5);
        
        % All perturbed (up or down, all sign)
        %         pert_rxns2{i,j}=c1(or(abs(flux2(c3)./flux1(c2))>=2,abs(flux2(c3)./flux1(c2))<=0.5)); % Bone and Brain vs Liver
        
    end
    
    clear str1 str2 str3 str4 c1 c2 c3 rxns1 rxns2 up_idx up_idx1 down_idx down_idx1 flux1 flux2 i j
    
end

clear category1 category2 up_idx up_idx1 down_idx down_idx1 same_sign_idx opp_sign_idx fold_change

%%%% Same sign
P1_P3_pert_rxns_Breast_same_sign=intersect(intersect(intersect(pert_rxns{1,1},pert_rxns{1,2}),pert_rxns{2,1}),pert_rxns{2,2});

k=1;
for i=1:2
    for j=1:2
        k=k+1;
        [p1,~,p2]=intersect(P1_P3_pert_rxns_Breast_same_sign,string(rxns{i,j}));
        P1_P3_pert_rxns_Breast_same_sign(:,1)=p1;
        P1_P3_pert_rxns_Breast_same_sign(:,k)=Flux_foldchange{i,j}(p2,1);
    end
end

clear i j k p1 p2

P1_P3_up_rxns_Breast_same_sign=intersect(intersect(intersect(up_rxns{1,1},up_rxns{1,2}),up_rxns{2,1}),up_rxns{2,2});

k=1;
for i=1:2
    for j=1:2
        k=k+1;
        [p1,~,p2]=intersect(P1_P3_up_rxns_Breast_same_sign,string(rxns{i,j}));
        P1_P3_up_rxns_Breast_same_sign(:,1)=p1;
        P1_P3_up_rxns_Breast_same_sign(:,k)=Flux_foldchange{i,j}(p2,1);
    end
end

clear i j k p1 p2

P1_P3_down_rxns_Breast_same_sign=intersect(intersect(intersect(down_rxns{1,1},down_rxns{1,2}),down_rxns{2,1}),down_rxns{2,2});

k=1;
for i=1:2
    for j=1:2
        k=k+1;
        [p1,~,p2]=intersect(P1_P3_down_rxns_Breast_same_sign,string(rxns{i,j}));
        P1_P3_down_rxns_Breast_same_sign(:,1)=p1;
        P1_P3_down_rxns_Breast_same_sign(:,k)=Flux_foldchange{i,j}(p2,1);
    end
end

clear i j k p1 p2


%%%% Opposite sign
P1_P3_pert_rxns_Breast_opp_sign=intersect(intersect(intersect(pert_rxns1{1,1},pert_rxns1{1,2}),pert_rxns1{2,1}),pert_rxns1{2,2});

k=1;
for i=1:2
    for j=1:2
        k=k+1;
        [p1,~,p2]=intersect(P1_P3_pert_rxns_Breast_opp_sign,string(rxns{i,j}));
        P1_P3_pert_rxns_Breast_opp_sign(:,1)=p1;
        P1_P3_pert_rxns_Breast_opp_sign(:,k)=Flux_foldchange{i,j}(p2,1);
    end
end

clear i j k p1 p2


P1_P3_up_rxns_Breast_opp_sign=intersect(intersect(intersect(up_rxns1{1,1},up_rxns1{1,2}),up_rxns1{2,1}),up_rxns1{2,2});

k=1;
for i=1:2
    for j=1:2
        k=k+1;
        [p1,~,p2]=intersect(P1_P3_up_rxns_Breast_opp_sign,string(rxns{i,j}));
        P1_P3_up_rxns_Breast_opp_sign(:,1)=p1;
        P1_P3_up_rxns_Breast_opp_sign(:,k)=Flux_foldchange{i,j}(p2,1);
    end
end

clear i j k p1 p2

P1_P3_down_rxns_Breast_opp_sign=intersect(intersect(intersect(down_rxns1{1,1},down_rxns1{1,2}),down_rxns1{2,1}),down_rxns1{2,2});


k=1;
for i=1:2
    for j=1:2
        k=k+1;
        [p1,~,p2]=intersect(P1_P3_down_rxns_Breast_opp_sign,string(rxns{i,j}));
        P1_P3_down_rxns_Breast_opp_sign(:,1)=p1;
        P1_P3_down_rxns_Breast_opp_sign(:,k)=Flux_foldchange{i,j}(p2,1);
    end
end

clear i j k p1 p2

%%%% Total perturbation
P1_P3_pert_rxns_Breast_total=[P1_P3_pert_rxns_Breast_same_sign; P1_P3_pert_rxns_Breast_opp_sign];

P1_P3_up_rxns_Breast_total=[P1_P3_up_rxns_Breast_same_sign; P1_P3_up_rxns_Breast_opp_sign];

P1_P3_down_rxns_Breast_total=[P1_P3_down_rxns_Breast_same_sign; P1_P3_down_rxns_Breast_opp_sign];


% P1_P3_pert_rxns2_Breast=intersect(intersect(intersect(pert_rxns2{1,1},pert_rxns2{1,2}),pert_rxns2{2,1}),pert_rxns2{2,2});
% 
% P1_P3_up_rxns2_Breast=intersect(intersect(intersect(up_rxns2{1,1},up_rxns2{1,2}),up_rxns2{2,1}),up_rxns2{2,2});
% 
% P1_P3_down_rxns2_Breast=intersect(intersect(intersect(down_rxns2{1,1},down_rxns2{1,2}),down_rxns2{2,1}),down_rxns2{2,2});


% save P1_P3_pert_rxns_Breast_same_sign P1_P3_pert_rxns_Breast_same_sign
% save P1_P3_up_rxns_Breast_same_sign P1_P3_up_rxns_Breast_same_sign
% save P1_P3_down_rxns_Breast_same_sign P1_P3_down_rxns_Breast_same_sign
% 
% save P1_P3_pert_rxns_Breast_opp_sign P1_P3_pert_rxns_Breast_opp_sign
% save P1_P3_up_rxns_Breast_opp_sign P1_P3_up_rxns_Breast_opp_sign
% save P1_P3_down_rxns_Breast_opp_sign P1_P3_down_rxns_Breast_opp_sign
% 
% save P1_P3_pert_rxns_Breast_total P1_P3_pert_rxns_Breast_total
% save P1_P3_up_rxns_Breast_total P1_P3_up_rxns_Breast_total
% save P1_P3_down_rxns_Breast_total P1_P3_down_rxns_Breast_total

clear rxns up_rxns up_rxns1 down_rxns down_rxns1 pert_rxns pert_rxns1 Flux_foldchange


%% Changes between patient 2 (UCD4) and patient 3 (UCD65)
category1={'P2_Breast_tINIT.clst1', 'P2_Breast_tINIT.clst2'}; % liver organotropism
category2={'P3_Breast_tINIT.clst1', 'P3_Breast_tINIT.clst2'}; % bone and brain organotropism

rxns={};
pert_rxns={};
pert_rxns1={};
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
        
        % same_sign_idx=flux2(c3).*flux1(c2)>=0; % True if both fluxes have the same sign
        % opp_sign_idx=flux2(c3).*flux1(c2)<=0;
        
        same_sign_idx=(flux1(c2)>0 & flux2(c3)>=0) | (flux2(c3)>0 & flux1(c2)>=0);
        opp_sign_idx=(flux1(c2)<0 & flux2(c3)<=0) | (flux2(c3)<0 & flux1(c2)<=0);
        
        rxns{i,j}=c1;
        
        % Fold change
        fold_change=flux2(c3)./flux1(c2);
        
        Flux_foldchange{i,j}=fold_change;
        
        
        %% Perturbed and same sign
        
        % Upregulated & same: fold change > 2 and same sign
        up_idx=fold_change>=2 & same_sign_idx;
        
        % Downregulated & same: fold change < 0.5 and same sign
        down_idx=fold_change<=0.5 & same_sign_idx;
        
        up_rxns{i,j}=c1(up_idx);
        down_rxns{i,j}=c1(down_idx);
        
        % Perturbed (up or down, same sign only)
        pert_rxns{i,j}=c1(up_idx | down_idx);
        
        
        %% Perturbed and opposite sign
        
        % Upregulated & opp: fold change > 2 and opposite sign
        up_idx1=fold_change>=2 & opp_sign_idx;
        
        % Downregulated & opp: fold change < 0.5 and opposite sign
        down_idx1=fold_change<=0.5 & opp_sign_idx;
        
        up_rxns1{i,j}=c1(up_idx1);
        down_rxns1{i,j}=c1(down_idx1);
        
        % Perturbed (up or down, opp sign only)
        pert_rxns1{i,j}=c1(up_idx1 | down_idx1);
        
        
        %% All Perturbed
        % All up reactions
        %         up_rxns2{i,j}=c1(abs(flux2(c3)./flux1(c2))>=2);
        
        % All down reactions
        %         down_rxns2{i,j}=c1(abs(flux2(c3)./flux1(c2))<=0.5);
        
        % All perturbed (up or down, all sign)
        %         pert_rxns2{i,j}=c1(or(abs(flux2(c3)./flux1(c2))>=2,abs(flux2(c3)./flux1(c2))<=0.5)); % Bone and Brain vs Liver
        
    end
    
    clear str1 str2 str3 str4 c1 c2 c3 rxns1 rxns2 up_idx up_idx1 down_idx down_idx1 flux1 flux2 i j
    
end

clear category1 category2 up_idx up_idx1 down_idx down_idx1 same_sign_idx opp_sign_idx fold_change


%%%% Same sign
P2_P3_pert_rxns_Breast_same_sign=intersect(intersect(intersect(pert_rxns{1,1},pert_rxns{1,2}),pert_rxns{2,1}),pert_rxns{2,2});

k=1;
for i=1:2
    for j=1:2
        k=k+1;
        [p1,~,p2]=intersect(P2_P3_pert_rxns_Breast_same_sign,string(rxns{i,j}));
        P2_P3_pert_rxns_Breast_same_sign(:,1)=p1;
        P2_P3_pert_rxns_Breast_same_sign(:,k)=Flux_foldchange{i,j}(p2,1);
    end
end

clear i j k p1 p2

P2_P3_up_rxns_Breast_same_sign=intersect(intersect(intersect(up_rxns{1,1},up_rxns{1,2}),up_rxns{2,1}),up_rxns{2,2});

k=1;
for i=1:2
    for j=1:2
        k=k+1;
        [p1,~,p2]=intersect(P2_P3_up_rxns_Breast_same_sign,string(rxns{i,j}));
        P2_P3_pert_up_Breast_same_sign(:,1)=p1;
        P2_P3_pert_up_Breast_same_sign(:,k)=Flux_foldchange{i,j}(p2,1);
    end
end

clear i j k p1 p2

P2_P3_down_rxns_Breast_same_sign=intersect(intersect(intersect(down_rxns{1,1},down_rxns{1,2}),down_rxns{2,1}),down_rxns{2,2});


k=1;
for i=1:2
    for j=1:2
        k=k+1;
        [p1,~,p2]=intersect(P2_P3_down_rxns_Breast_same_sign,string(rxns{i,j}));
        P2_P3_down_rxns_Breast_same_sign(:,1)=p1;
        P2_P3_down_rxns_Breast_same_sign(:,k)=Flux_foldchange{i,j}(p2,1);
    end
end

clear i j k p1 p2

%%%% Opposite sign
P2_P3_pert_rxns_Breast_opp_sign=intersect(intersect(intersect(pert_rxns1{1,1},pert_rxns1{1,2}),pert_rxns1{2,1}),pert_rxns1{2,2});

k=1;
for i=1:2
    for j=1:2
        k=k+1;
        [p1,~,p2]=intersect(P2_P3_pert_rxns_Breast_opp_sign,string(rxns{i,j}));
        P2_P3_pert_rxns_Breast_opp_sign(:,1)=p1;
        P2_P3_pert_rxns_Breast_opp_sign(:,k)=Flux_foldchange{i,j}(p2,1);
    end
end

clear i j k p1 p2


P2_P3_up_rxns_Breast_opp_sign=intersect(intersect(intersect(up_rxns1{1,1},up_rxns1{1,2}),up_rxns1{2,1}),up_rxns1{2,2});

k=1;
for i=1:2
    for j=1:2
        k=k+1;
        [p1,~,p2]=intersect(P2_P3_up_rxns_Breast_opp_sign,string(rxns{i,j}));
        P2_P3_up_rxns_Breast_opp_sign(:,1)=p1;
        P2_P3_up_rxns_Breast_opp_sign(:,k)=Flux_foldchange{i,j}(p2,1);
    end
end

clear i j k p1 p2


P2_P3_down_rxns_Breast_opp_sign=intersect(intersect(intersect(down_rxns1{1,1},down_rxns1{1,2}),down_rxns1{2,1}),down_rxns1{2,2});


k=1;
for i=1:2
    for j=1:2
        k=k+1;
        [p1,~,p2]=intersect(P2_P3_down_rxns_Breast_opp_sign,string(rxns{i,j}));
        P2_P3_down_rxns_Breast_opp_sign(:,1)=p1;
        P2_P3_down_rxns_Breast_opp_sign(:,k)=Flux_foldchange{i,j}(p2,1);
    end
end

clear i j k p1 p2

%%%% Total perturbation
P2_P3_pert_rxns_Breast_total=union(P2_P3_pert_rxns_Breast_same_sign, P2_P3_pert_rxns_Breast_opp_sign);

P2_P3_up_rxns_Breast_total=union(P2_P3_up_rxns_Breast_same_sign, P2_P3_up_rxns_Breast_opp_sign);

P2_P3_down_rxns_Breast_total=union(P2_P3_down_rxns_Breast_same_sign, P2_P3_down_rxns_Breast_opp_sign);


% P2_P3_pert_rxns2_Breast=intersect(intersect(intersect(pert_rxns2{1,1},pert_rxns2{1,2}),pert_rxns2{2,1}),pert_rxns2{2,2});
% 
% P2_P3_up_rxns2_Breast=intersect(intersect(intersect(up_rxns2{1,1},up_rxns2{1,2}),up_rxns2{2,1}),up_rxns2{2,2});
% 
% P2_P3_down_rxns2_Breast=intersect(intersect(intersect(down_rxns2{1,1},down_rxns2{1,2}),down_rxns2{2,1}),down_rxns2{2,2});


% save P2_P3_pert_rxns_Breast_same_sign P2_P3_pert_rxns_Breast_same_sign
% save P2_P3_up_rxns_Breast_same_sign P2_P3_up_rxns_Breast_same_sign
% save P2_P3_down_rxns_Breast_same_sign P2_P3_down_rxns_Breast_same_sign
% 
% save P2_P3_pert_rxns_Breast_opp_sign P2_P3_pert_rxns_Breast_opp_sign
% save P2_P3_up_rxns_Breast_opp_sign P2_P3_up_rxns_Breast_opp_sign
% save P2_P3_down_rxns_Breast_opp_sign P2_P3_down_rxns_Breast_opp_sign
% 
% save P2_P3_pert_rxns_Breast_total P2_P3_pert_rxns_Breast_total
% save P2_P3_up_rxns_Breast_total P2_P3_up_rxns_Breast_total
% save P2_P3_down_rxns_Breast_total P2_P3_down_rxns_Breast_total

clear rxns up_rxns up_rxns1 down_rxns down_rxns1 pert_rxns pert_rxns1 Flux_foldchange


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Common between P1_P3_pert_rxns_Breast and P2_P3_pert_rxns_Breast
[p1,p2,p3]=intersect(P1_P3_pert_rxns_Breast_same_sign(:,1), P2_P3_pert_rxns_Breast_same_sign(:,1));
Primary_level_pert_rxns(:,1)=p1;
Primary_level_pert_rxns(:,2:5)=P1_P3_pert_rxns_Breast_same_sign(p2,2:end);
Primary_level_pert_rxns(:,6:9)=P2_P3_pert_rxns_Breast_same_sign(p3,2:end);

l1=[];
for i=1:size(Primary_level_pert_rxns,1)
    l1(i,1)=length(find(double(Primary_level_pert_rxns(i,2:end))>2));
end

mean_fold_change=[];
final_rxns=strings();
nat=strings();
n=1;
for i=1:size(Primary_level_pert_rxns,1)
    idx=find(double(Primary_level_pert_rxns(i,2:end))>2);
    idx1=find(double(Primary_level_pert_rxns(i,2:end))<(1/2));
    l=length(idx);
    if l>=5
        final_rxns(n,1)=Primary_level_pert_rxns(i,1);
        mean_fold_change(n,1)=mean(double(Primary_level_pert_rxns(i,idx+1)));
        nat(n,1)='Up-regulated';
        n=n+1;
    elseif l<4
        final_rxns(n,1)=Primary_level_pert_rxns(i,1);
        mean_fold_change(n,1)=mean(double(Primary_level_pert_rxns(i,idx1+1)));
        nat(n,1)='Down-regulated';
        n=n+1;
    end
end
Primary_level_pert_rxns_final=[final_rxns, mean_fold_change, nat];

save Primary_level_pert_rxns_final Primary_level_pert_rxns_final


%% Finding how many are up-regulated and down-regulated
length(find(Primary_level_pert_rxns_final(:,3)=='Up-regulated')); % 694 reactions are Up

length(find(Primary_level_pert_rxns_final(:,3)=='Down-regulated')); % 248 reactions are Down


%% Pathway analysis
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Metabolic_models\Mouse-GEM.mat');

model=mouseGEM;
[~,~,rxnSet]=intersect(Primary_level_pert_rxns_final(:,1),string(model.rxns(:,1)));
group='subSystems';
pathways= FEA(model, rxnSet, group);

find(cell2mat(pathways(2:end,1))<0.05);
sig_path=[pathways(1,:); pathways(ans+1,:)];
Pathway_analysis.Pathways_primary_level=sig_path;


%% Finding distribution (up, down) of perturbed reactions in each perturbed pathways
[~,~,rxnSet]=intersect(Primary_level_pert_rxns_final(:,1),string(model.rxns(:,1)));
rxns=[string(mouseGEM.rxns(rxnSet)),string(mouseGEM.subSystems(rxnSet))];
path=Pathway_analysis.Pathways_primary_level(2:end,3);

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

