function coh_avg = erp_grandaverage(obj)
% Calculate grand average ERPs

% Number of conditions
Ncond = length(obj.conditionlabels) ; 

% initialize output
coh_avg = microstate.cohort ; 

% Calculate group average ERPs
for i = 1:Ncond % Loop over conditions
    ind = find(obj.condition == i) ;
    ERP = [obj.individual(ind).data] ;
    ERP = reshape(ERP,size(ERP,1),size(ERP,2)/length(ind),length(ind)) ;
    ERP = mean(ERP,3) ; 
    
    ms = microstate.individual(ERP,obj.individual(1).modality,obj.individual(1).time) ; 
    coh_avg = coh_avg.add_individuals(ms,obj.conditionlabels{i},1) ; 
end    
 

    
   