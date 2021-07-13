function [p,info] = erp_clusterperm_TANOVA(obj,design,varargin)
% Perform cluster permutation TANOVA analysis

%% Check inputs
% Check multiple conditions appear in the data
if length(obj.conditionlabels) == 1
    error('The cohort structure must contain multiple conditions'); 
end

% Get pairs of conditions - check if specfied
conditions = [] ; 
idx = find(strcmp(varargin,'conditions')) ; 
if ~isempty(idx)
    conditions = varargin{idx+1} ; 
elseif length(obj.conditionlabels) == 2
    % If conditions are not specified, check if only two conditions
    conditions = obj.conditionlabels ; 
else
    error('>2 conditions in data, so the pair of conditions to be tested must be specified') ; 
end

% Convert conditions to numbers if not already
if iscell(conditions)
    condlbl = conditions ; 
    conditions = zeros(1,2) ; 
    conditions(1) = find(strcmp(obj.conditionlabels,condlbl{1})) ; 
    conditions(2) = find(strcmp(obj.conditionlabels,condlbl{2})) ; 
end

% Check if number of permutations is input
idx = find(strcmp(varargin,'numpermuation')) ; 
if ~isempty(idx)
    Nperms = varargin{idx+1} ; 
else
    Nperms = 2000 ; 
end

% Check if tstart or tend are input
idx = find(strcmp(varargin,'tstart')) ; 
if ~isempty(idx)
    tstart = varargin{idx+1} ; 
else
    tstart = 0 ; 
end

idx = find(strcmp(varargin,'tend')) ; 
if ~isempty(idx)
    tend = varargin{idx+1} ; 
else
    tend = inf ; 
end

% Make a within or between design flag
switch design
    case 'between'
        withinflag = false ; 
    case 'within'
        withinflag = true ; 
    otherwise
        error('design must be either between or within')  
end


% Select the correct time samples
obj = microstate.functions.select_time(obj,[tstart,tend]) ; 
    
    
%% Run the analysis
        
% Get the maximal cluster statistic
[t,thresh] = DISSstat(obj,conditions) ; 
[T,samples] = max_cluster_stat(t,thresh) ; 

% Run permutations
Ts = nan(1,Nperms) ; Ts(1) = T ;
ts = nan(Nperms,length(t)) ; ts(1,:) = t ; 
msg = [] ; 
for perm = 2:Nperms

    if ~mod(perm,100)
        fprintf(repmat('\b',1,length(msg)))
        msg = sprintf('Cluster permutation testing (TANOVA), permutation %d of %d',perm,Nperms) ; 
        fprintf(msg) ; 
    end

    % Permute condition labels
    if ~withinflag
        objp = obj ; 
        ind = find( (objp.condition == conditions(1)) | (objp.condition == conditions(2)) ) ; 
        cond = objp.condition(ind) ;
        cond = cond(randperm(length(cond))) ; 
        objp.condition(ind) = cond ; 
    else
        objp = obj ; 
        ind1 = find((objp.condition == conditions(1))) ; 
        ind2 = find((objp.condition == conditions(2))) ; 
        if length(ind1)~=length(ind2)
            error('For within design, the same number of individuals are required in each condition (and should be in the same order)')
        end
        cond = [objp.condition(ind1) , objp.condition(ind2)] ; 
        ind_perm = rand(size(cond,1),1)<0.5 ;
        cond_perm = cond ; 
        cond_perm(ind_perm,1) = cond(ind_perm,2) ; 
        cond_perm(ind_perm,2) = cond(ind_perm,1) ; 
        objp.condition(ind1) = cond_perm(:,1) ; 
        objp.condition(ind2) = cond_perm(:,2) ; 
    end

    % Get the maximal cluster statistic
    [ts(perm,:),~] = DISSstat(objp,conditions) ; 
    Ts(perm) = max_cluster_stat(ts(perm,:),thresh) ; 

end
fprintf(repmat('\b',1,length(msg)))
p = 1-(sum(T>Ts)/Nperms) ; 
pms = 1-(sum(repmat(t,Nperms,1)>ts))/Nperms ; 

info = struct ; 
info.time = obj.individual(1).time ; 
info.DISS = t ; 
info.threshold = thresh ; 
info.maxcluster_time = obj.individual(1).time([samples(1),samples(end)]) ; 
info.p_sample = pms ; 



end


%% Helper functions
function [t,thresh] = DISSstat(obj,conditions)

    % Calculate group average ERPs
    for i = 1:2 % Loop over conditions
        ind = find(obj.condition == conditions(i)) ;
        ERP{i} = [obj.individual(ind).data] ;
        ERP{i} = reshape(ERP{i},size(ERP{i},1),size(ERP{i},2)/length(ind),length(ind)) ;
        ERP{i} = mean(ERP{i},3) ; 
    end
    
    % Calculate DISS
    t = diag(microstate.functions.DISS(ERP{1},ERP{2},obj.individual(1).modality))' ;

    % Get threshold
    thresh = min(t) + 0.5*(max(t)-min(t)) ;     
    
end
    
function [T,samples] = max_cluster_stat(t,thresh)

    ind = find(t>thresh) ;
    if isempty(ind)
        T = 0 ; samples = [] ; 
        return
    end
    dx = [0,find(diff(ind)>1),length(ind)] ; 
    
    for i = 1:length(dx)-1
        
        indi = ind((dx(i)+1):dx(i+1)) ; 
        T(i) = sum(t(indi)) ; 
        
        index_onoff(:,i) = [indi(1);indi(end)] ; 
        
    end
    
    [T,index_max] = max(T)  ; 
    index_onoff = index_onoff(:,index_max) ; 
    samples = index_onoff(1):index_onoff(2); 
    
end

    
   