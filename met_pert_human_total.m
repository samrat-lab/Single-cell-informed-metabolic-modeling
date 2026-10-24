function [final_pathlist,pval,overlap]=met_pert_human_total(input)
load('D:\CANCER\Breast_Organotropism\Matlab_variables_Organotropism\Human-GEM-main\model\Human-GEM.mat')

funn_sys=unique(string(ihuman.subSystems));
aa=string(ihuman.subSystems);
for i=1:length(funn_sys)
    nu(i)=length(strmatch(funn_sys(i),aa,'exact'));
end

nu_11={};
nu_12={};
aa1=string(ihuman.subSystems(input));
aa1_un=unique(aa1);
for i=1:length(aa1_un)
    
    nu_1(i)=length(strmatch(aa1_un(i),aa1,'exact'));
    olap{i}=[num2str(length(strmatch(aa1_un(i),aa1,'exact'))),'','\',num2str(length(strmatch(aa1_un(i),aa,'exact')))];
end

for i=1:length(aa1_un)
    [~,~,ix]=intersect(aa1_un(i),funn_sys); % How many times the unique pathways repeating in ihuman.subSystems
    hyp=hygepdf(1:nu_1(i),size(ihuman.rxns,1),nu(ix),length(input)); % Hypergeometric proability density function
    % Y= hygepdf(X,M,K,N) - This above function computes the hypergeometric probability density function
    % at each of the values in X using the corresponding size of the population M, number of items with the desired
    % characteristics in the population K, and the number of samples drawn, N.
    
    hyp1(i)=hyp(end);
end
% pval=hyp1(find(hyp1<.05))';
pval=hyp1';
[pval,b1]=sort(pval);
% final_pathlist=aa1_un(find(hyp1<.05));
final_pathlist=aa1_un;
final_pathlist=final_pathlist(b1);
% overlap=olap(find(hyp1<.05))';
overlap=olap';
overlap=overlap(b1);

return