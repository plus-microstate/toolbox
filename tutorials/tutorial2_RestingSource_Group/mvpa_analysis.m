function [p,accuracy] = mvpa_analysis(nets) ; 

% Get degree distributions
degree_dist = [] ; 
tbl = table([],[],[],'VariableNames',{'Dataset','Microstate','Degree_Distribution'}) ; 
for i=1:size(nets,1)
    for j=1:size(nets,2)
        degree_dist_ij = squeeze(nansum(nets{i,j},1))' ;
        if size(degree_dist_ij,2)==1 && size(degree_dist_ij,1)==230
            degree_dist_ij = degree_dist_ij' ; 
        end
        N = size(tbl,1) ; 
        M = size(degree_dist_ij,1) ;
        tbl_ij = table(i*ones(M,1),j*ones(M,1),degree_dist_ij,'VariableNames',{'Dataset','Microstate','Degree_Distribution'}) ; 
        tbl(N+(1:M),:) = tbl_ij ; 
    end
end

%% Do machine learning

Npermutations = 200 ; 
acc = nan(Npermutations,1) ; 

mdl = fitcdiscr(tbl.Degree_Distribution,tbl.Microstate,'CrossVal','on','KFold',10) ; 
acc(1) = 1-kfoldLoss(mdl,'LossFun','ClassifError') ; 

for perm = 2:Npermutations
    
    permuted_Microstate = tbl.Microstate(randperm(size(tbl,1))) ;
    
    mdl_perm = fitcdiscr(tbl.Degree_Distribution,permuted_Microstate,'CrossVal','on','KFold',10) ; 
    acc(perm) = 1-kfoldLoss(mdl_perm,'LossFun','ClassifError') ; 

end
    
p = 1-sum(acc(1)>acc)/Npermutations ; 
accuracy = acc(1) ; 
end