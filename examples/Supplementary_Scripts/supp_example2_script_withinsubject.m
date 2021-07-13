%% Format the data into a microstate cohort
% For a within-individual design, each participant and each condition is an
% individual. Importantly for within-individual design, in the cohort you
% MUST make sure to upload one of each condition for each participant *in
% the same order for all conditions*, e.g. if there are two conditions, the
% first individual object from condition 1 and the first individual object
% from condition 2 must be from the same participant. 
% 
% The (preprocessed) data for this example can be downloaded from: 
% ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/tutorial/cluster_permutation_timelock/ERF_orig.mat

clear, clc, close all
rng('default')

% Load in the data
load('ERF_orig')

% Initialize a cohort
coh = microstate.cohort ; 

% Loop over participants, and add each of their conditions
for i = 1:length(allsubjFC)
    % Add incongruent trials
    ms = microstate.individual(allsubjFIC{i}.avg','meg',allsubjFIC{i}.time) ; 
    coh = coh.add_individuals(ms,'Incongruent',1) ; 
    
    % Add congruent trials
    ms = microstate.individual(allsubjFC{i}.avg','meg',allsubjFC{i}.time) ; 
    coh = coh.add_individuals(ms,'Congruent',1) ; 
end

% Following Fieldtrip, we will analyse 0-1 seconds
coh = microstate.functions.select_time(coh,[0,1]) ; 

%% perform topographic ERP millisecond-by-millisecond GFP analysis

coh = coh.ind_calculate_gfp ; % calculate GFP
[p_GFP,info_GFP] = coh.erp_clusterperm_gfp('within') ; % cluster permutation test on GFP

%% perform topographic ERP millisecond-by-millisecond TANOVA analysis

[p_TANOVA,info_TANOVA] = coh.erp_clusterperm_TANOVA('within') ; % cluster permutation test on GFP