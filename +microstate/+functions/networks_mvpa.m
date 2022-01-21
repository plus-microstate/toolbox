function [p,acc,confusionmat] = networks_mvpa(nets) ; 

% If multiple participtants, append
if size(nets,1)>1
    for i = 1:size(nets,2)
        
        if isempty(nets{1,i})
            j = 1 ; 
            while isempty(nets{j,i})
                j=j+1 ; 
            end
            nets{1,i} = nets{j,i} ; nets{j,i} = [] ; 
        end
        for j = 2:size(nets,1)
            nets{1,i} = cat(3,nets{1,i},nets{j,i}) ; 
            nets{j,i} = [] ; 
        end
    end
    nets = nets(1,:) ; 
end
            

% Get degree distributions
degree_dist = [] ; 
tbl = table([],[],'VariableNames',{'Microstate','Degree_Distribution'}) ; 
for i=1:length(nets)
    degree_dist_i = squeeze(nansum(nets{i},1))' ;
    if size(degree_dist_i,2)==1
        degree_dist_i = degree_dist_i' ; 
    end
    N = size(tbl,1) ; 
    M = size(degree_dist_i,1) ;
    tbl_i = table(i*ones(M,1),degree_dist_i,'VariableNames',{'Microstate','Degree_Distribution'}) ; 
    tbl(N+(1:M),:) = tbl_i ; 
end

%% Do machine learning

Npermutations = 200 ; 
acc = nan(Npermutations,1) ; 

mdl = fitcdiscr(tbl.Degree_Distribution,tbl.Microstate,'CrossVal','on','KFold',10) ; 
acc(1) = 1-kfoldLoss(mdl,'LossFun','ClassifError') ; 

% Make confusion matrix
pred = kfoldPredict(mdl) ; 
confusionmat = zeros(length(nets)) ; 
for i = 1:length(nets)
    for j = 1:length(nets)
        confusionmat(i,j) = sum( tbl.Microstate==i & pred==j ) ; 
    end
end

for perm = 2:Npermutations
    
    permuted_Microstate = tbl.Microstate(randperm(size(tbl,1))) ;
    
    mdl_perm = fitcdiscr(tbl.Degree_Distribution,permuted_Microstate,'CrossVal','on','KFold',10) ; 
    acc(perm) = 1-kfoldLoss(mdl_perm,'LossFun','ClassifError') ; 

end
    
p = 1-sum(acc(1)>acc)/Npermutations ; 
accuracy = acc(1) ; 
end