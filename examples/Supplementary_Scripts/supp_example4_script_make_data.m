%% Step 1: 5000 GFP peaks from first scan and perform clustering

cohort_cluster = microstate.cohort ; % make microstate cohort
for i = 1:length(subject_IDs)
    % load the source data
    source = load(subject_IDs{i}) ;
    
    % get the bad samples from fieldtrip's automated artifact detection and
    % convert from a fieldtrip artfctdef structure to an array of samples
    load([subject_IDs{i},'_artfcts'],'artfctdef') % load artfctdef structure
    bad_samples = [] ; % initialize empty array of bad samples
    for mth = {'clip','jump','zscore'} % loop over types of artfct
        for j = 1:size(artfctdef.(mth{1}).artifact,1) % append artifacts
            bad_samples = [bad_samples , (artfctdef.(mth{1}).artifact(j,1)-5):(artfctdef.(mth{1}).artifact(j,2)+5)] ; 
        end
    end
    bad_samples = unique(bad_samples) ; % in case of overlap
      
    % for clustering, use only scan 1
    % make a microstate object
    ms = microstate.individual(source.trial{1}','source',source.time{1}) ; % make microstate individual object
    bad_samples_scan1 = bad_samples(...
                    (bad_samples>=source.sampleinfo(1,1))...
                    &(bad_samples<=source.sampleinfo(1,2))...
                ) - source.sampleinfo(1,1)+1 ; 
    ms = ms.add_bad_samples(bad_samples_scan1) ; 
            
    % add the microstate object to the cohort with 5000 GFP peaks
    cohort_cluster = cohort_cluster.add_individuals(ms,sesh,5000) ; 
end
  
% Perform clustering
cohort_cluster = cohort_cluster.cluster_globalkoptimum('kvec',2:40) ; 
maps = cohort_cluster.globalmaps ; 
clear cohort_cluster

%% Step 2: Backfit to all scans

cohort_task = microstate.cohort ; 
for i = 1:length(subject_IDs) 
    
    % load the source data
    source = load(subject_IDs{i}) ;

    for scan = 1:4
        % make a microstate object
        ms = microstate.individual(source.trial{scan}','source',source.time{scan}) ; % make microstate individual object

        % backfit maps
        ms.maps = maps ;
        ms = ms.cluster_alignmaps ; 
                 
        % add to cohort
        cohort_task = cohort_task.add_individuals(ms,sprintf('Continuous data, participant %d',i),0) ;
    end
end