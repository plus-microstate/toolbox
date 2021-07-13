function info = erp_chi2stateprobs(obj,expectedConditions,observedConditions,varargin)
% Perform chi2 test for responses

%% Check inputs
% Check multiple conditions appear in the data
if length(obj.conditionlabels) == 1
    error('The cohort structure must contain multiple conditions'); 
end

% Check if number of permutations is input
idx = find(strcmp(varargin,'numpermutation')) ; 
if ~isempty(idx)
    Nperms = varargin{idx+1} ; 
else
    Nperms = 2000 ; 
end

% Convert expected and baseline conditios to cell if only one specified
if ~iscell(expectedConditions); expectedConditions = {expectedConditions} ; end
if ~iscell(observedConditions); observedConditions = {observedConditions} ; end

% Check there are trials with expected conditions and convert to numeric
if length(expectedConditions) == 1 & strcmpi(expectedConditions{1},'prestim')
    expectedConditions = -1 ; 
elseif length(expectedConditions) > 1 & any(strcmpi(expectedConditions{1},'prestim'))
    error('Cannot specify both prestimulus and a post-stimulus condition as expected')
else
    tmp = nan(1,length(expectedConditions)) ;  ; 
    for i = 1:length(expectedConditions)
        ind = find(strcmpi(expectedConditions{i},obj.conditionlabels)) ; 
        if isempty(ind)
            error('Cannot find condition %s',expectedConditions{i}) 
        end
        tmp(i) = ind ; 
    end
    expectedConditions = tmp ;
end
 

% Check there are trials with observed conditions and convert to numeric
tmp = nan(1,length(observedConditions)) ;  ; 
for i = 1:length(observedConditions)
    ind = find(strcmpi(observedConditions{i},obj.conditionlabels)) ; 
    if isempty(ind)
        error('Cannot find condition %s',observedConditions{i}) 
    end
    tmp(i) = ind ;
end
observedConditions = tmp ; 

% Get latencies
idx = find(strcmpi(varargin,'latency'));
if isempty(idx) ;
    latency = [0,inf] ; 
else 
    latency = varargin{idx+1} ; 
end

idx = find(strcmpi(varargin,'prestim_latency')) ; 
if isempty(idx)
    prestim_latency = [-inf,0] ; 
else
    prestim_latency = varargin{idx+1} ; 
end

%% Check all trials have the same time axis and concatenate sequences

time = {obj.individual(:).time} ;
try
    time = cell2mat(time(:)) ; 
    dtime = diff(time) ; 
    if any(dtime(:)~=0)
        error('tmperror')
    end
    time = time(1,:) ; 
catch
    error('All trials must be sampled along the same time axis') ; 
end


msseq = {obj.individual(:).label} ; 
msseq = cell2mat(msseq(:)) ; 



    
%% Observed Counts

% Extract the correct trials
observedTrials = any(repmat(obj.condition,1,length(observedConditions)) == observedConditions,2) ; 
numObservedTrials = sum(observedTrials) ; 

% Loop over microstates
numMaps = max(msseq(:)) ; 
indTime = (time >= latency(1))&(time<=latency(2)) ; 
msseqObs = msseq(observedTrials,:) ;
for i = 1:numMaps
    observedCount(i,:) = sum(msseqObs == i) ; 
end
observedCount(:,~indTime) = nan ;   

%% Expected Counts

% Deal with pre-stim case
if expectedConditions == -1
    
    % Loop over microstates
    indTime = (time >= prestim_latency(1))&(time<=prestim_latency(2)) ; 
    for i = 1:numMaps
        expectedCount(i,:) = sum(msseq == i) ; 
    end
    
    % Just get prestim period
    expectedCount = expectedCount(:,indTime) ; 
    expectedCount = sum(expectedCount,2) ; 
    
    % Normalize to numObservedTrials
    expectedCount = numObservedTrials*(expectedCount/sum(expectedCount,1)) ;
    
    % Make correct size
    expectedCount = repmat(expectedCount,1,size(observedCount,2)) ; 
    indTime = (time >= latency(1))&(time<=latency(2)) ; 
    expectedCount(:,~indTime) = nan ;
    expectedCount(expectedCount == 0) = nan ; 
    
    
else
    % post-stim case
    
    % Extract the correct trials
    expectedTrials = any(repmat(obj.condition,1,length(expectedConditions)) == expectedConditions,2) ; 
    numExpectedTrials = sum(expectedTrials) ; 

    % Loop over microstates
    indTime = (time >= latency(1))&(time<=latency(2)) ; 
    msseqExp = msseq(expectedTrials,:) ;
    for i = 1:numMaps
        expectedCount(i,:) = sum(msseqExp == i) ; 
    end
    
    % Normalize to numObservedTrials
    expectedCount = (expectedCount/numExpectedTrials)*numObservedTrials ; 
    expectedCount(:,~indTime) = nan ;
    expectedCount(expectedCount == 0) = nan ; 
    
end

%% Calculate Pearson residuals and chi2 distance

pearsonResiduals = (observedCount-expectedCount)./sqrt(expectedCount) ;
chi2 = nansum(pearsonResiduals.^2) ; 
chi2(~indTime) = nan ; 
pval = chi2cdf(chi2,numMaps-1,'upper') ; 

%% Do cluster permutation analysis

% Get threshold - uncorrected p=0.05
chi2Threshold = fminsearch(@(chi2) abs(chi2cdf(chi2,numMaps-1,'upper')-0.05),20) ;

info = struct ; 
info.time = time ; 
info.pearsonResiduals = pearsonResiduals ; 
info.chi2 = chi2 ; 
info.threshold = chi2Threshold ; 
info.uncorrected_pvals = pval ; 

if Nperms >= 0 && expectedConditions ~= -1 


    % Find maximal cluster statistic
    [CHI2,samples] = max_cluster_stat(chi2,chi2Threshold) ; 
    CHI2Perm = nan(1,Nperms) ; CHI2Perm(1) = CHI2 ; 

    % Do permutations
    msg = [] ; 
    for perm = 2:Nperms
        
        if ~mod(perm,100)
            fprintf(repmat('\b',1,length(msg)))
            msg = sprintf('Cluster permutation testing (chi2), permutation %d of %d',perm,Nperms) ; 
            fprintf(msg) ; 
        end

        % Permute trials
        isaTrial = find(expectedTrials | observedTrials) ; 
        perminds = randperm(length(isaTrial)) ; 
        perminds = isaTrial(perminds) ; 

        expectedPerm = false(size(expectedTrials)) ; 
        expectedPerm(perminds(1:numExpectedTrials)) = true ; 

        observedPerm = false(size(observedTrials)) ; 
        observedPerm(perminds(numExpectedTrials+(1:numObservedTrials))) = true ; 

        % Observed count 
        msseqObs = msseq(observedPerm,:) ;
        for i = 1:numMaps
            observedCountPerm(i,:) = sum(msseqObs == i) ; 
        end
        observedCountPerm(:,~indTime) = nan ;   

        % Expected count
        msseqExp = msseq(expectedPerm,:) ;
        for i = 1:numMaps
            expectedCountPerm(i,:) = sum(msseqExp == i) ; 
        end

        % Normalize to numObservedTrials
        expectedCountPerm = (expectedCountPerm/numExpectedTrials)*numObservedTrials ; 
        expectedCountPerm(:,~indTime) = nan ;
        expectedCountPerm(expectedCountPerm == 0) = nan ; 

        % Chi2
        pearsonResidualsPerm = (observedCountPerm-expectedCountPerm)./sqrt(expectedCountPerm) ;
        chi2Perm = nansum(pearsonResidualsPerm.^2) ; 
        chi2Perm(~indTime) = nan ; 


        % Find maximal cluster statistic
        [CHI2Perm(perm),~] = max_cluster_stat(chi2Perm,chi2Threshold) ; 

    end

    p = 1-(sum(CHI2>CHI2Perm)/Nperms) ;
    info.maxcluster_time = time([samples(1),samples(end)]) ; 
    info.maxcluster_pval = p ; 
end





function [T,samples] = max_cluster_stat(t,thresh)

    ind = find(t>thresh) ;
    if isempty(ind)
        T = 0 ; samples = [nan,nan] ; 
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
end