%% Fieldtrip processing
% NOTE: This script does not contain any tutorials for +microstate toolbox,
% and is just used for pre-processing the data. For the tutorial on
% +microstate, see the script example_script_EEGcanonicalmaps.m
% 
% This script follows the Fieldtrip tutorial for cluster permutation
% analysis of ERPs [1] to pre-process and perform ERP analysis on  the
% example MEG ERP data set. You must have a version of Fieldtrip
% installed and the raw data downloaded [2] to run this script.
% 
% Full descriptions of the steps and each line of code are given in [1]. 
% This script is just included to get the data in a format ready for
% microstate analysis, and to plot ERP results. 
% 
% [1] https://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/
% [2] ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/tutorial/Subject01.zip

clear, clc, close all
rng('default') % for reproducibility

% Do the trial definition for all conditions together
cfg                         = [];
cfg.dataset                 = 'Subject01.ds';
cfg.trialfun                = 'ft_trialfun_general'; % this is the default
cfg.trialdef.eventtype      = 'backpanel trigger';
cfg.trialdef.eventvalue     = [3 5 9]; % the values of the stimulus trigger for the three conditions
% 3 = fully incongruent (FIC), 5 = initially congruent (IC), 9 = fully congruent (FC)
cfg.trialdef.prestim        = 1; % in seconds
cfg.trialdef.poststim       = 2; % in seconds

cfg = ft_definetrial(cfg);

% Cleaning
% remove the trials that have artifacts from the trl
cfg.trl([2, 5, 6, 8, 9, 10, 12, 39, 43, 46, 49, 52, 58, 84, 102, 107, 114, 115, 116, 119, 121, 123, 126, 127, 128, 133, 137, 143, 144, 147, 149, 158, 181, 229, 230, 233, 241, 243, 245, 250, 254, 260],:) = [];

% preprocess the data
cfg.channel    = {'MEG', '-MLP31', '-MLO12'};        % read all MEG channels except MLP31 and MLO12
cfg.demean     = 'yes';
cfg.baselinewindow  = [-0.2 0];
cfg.lpfilter   = 'yes';                              % apply lowpass filter
cfg.lpfreq     = 35;                                 % lowpass at 35 Hz.

data_all = ft_preprocessing(cfg);

% Split into Fully Congruent and Fully Incongruent
cfg = [];
cfg.trials = data_all.trialinfo == 3;
dataFIC_LP = ft_redefinetrial(cfg, data_all);
save('fieldtrip_output','dataFIC_LP')

cfg = [];
cfg.trials = data_all.trialinfo == 9;
dataFC_LP = ft_redefinetrial(cfg, data_all);
save('fieldtrip_output','dataFC_LP','-append')

%% Run statistics

cfg = [];
cfg.keeptrials = 'yes';
timelockFIC    = ft_timelockanalysis(cfg, dataFIC_LP);
timelockFC     = ft_timelockanalysis(cfg, dataFC_LP);


cfg                  = [];
cfg.method           = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability
cfg.statistic        = 'indepsamplesT'; % use the independent samples T-statistic as a measure to
                                   % evaluate the effect at the sample level
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;       % alpha level of the sample-specific test statistic that
                                   % will be used for thresholding
cfg.clusterstatistic = 'maxsum';   % test statistic that will be evaluated under the
                                   % permutation distribution.
cfg.minnbchan        = 2;          % minimum number of neighborhood channels that is
                                   % required for a selected sample to be included
                                   % in the clustering algorithm (default=0).
% cfg.neighbours     = neighbours; % see below
cfg.tail             = 0;          % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail      = 0;
cfg.alpha            = 0.025;      % alpha level of the permutation test
cfg.numrandomization = 100;        % number of draws from the permutation distribution

n_fc  = size(timelockFC.trial, 1);
n_fic = size(timelockFIC.trial, 1);

cfg.design           = [ones(1,n_fic), ones(1,n_fc)*2]; % design matrix
cfg.ivar             = 1; % number or list with indices indicating the independent variable(s)

cfg_neighb        = [];
cfg_neighb.method = 'distance';
neighbours        = ft_prepare_neighbours(cfg_neighb, dataFC_LP);

cfg.neighbours    = neighbours;  % the neighbours specify for each sensor with
                                 % which other sensors it can form clusters
cfg.channel       = {'MEG'};     % cell-array with selected channel labels
cfg.latency       = [0 1];       % time interval over which the experimental
                                 % conditions must be compared (in seconds)


[stat] = ft_timelockstatistics(cfg, timelockFIC, timelockFC);
save('fieldtrip_output','stat','-append')

%% Plot the data

cfg    = [];
avgFIC = ft_timelockanalysis(cfg, dataFIC_LP);
avgFC  = ft_timelockanalysis(cfg, dataFC_LP);

% Then take the difference of the averages using ft_math
cfg           = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
raweffectFICvsFC = ft_math(cfg, avgFIC, avgFC);

pos_cluster_pvals = [stat.posclusters(:).prob];

% Then, find which clusters are deemed interesting to visualize, here we use a cutoff criterion based on the
% cluster-associated p-value, and take a 5% two-sided cutoff (i.e. 0.025 for the positive and negative clusters,
% respectively
pos_clust = find(pos_cluster_pvals < 0.02);
pos       = ismember(stat.posclusterslabelmat, pos_clust);

% and now for the negative clusters...
neg_cluster_pvals = [stat.negclusters(:).prob];
neg_clust         = find(neg_cluster_pvals < 0.02);
neg               = ismember(stat.negclusterslabelmat, neg_clust);


% 
% pos = stat.posclusterslabelmat == 1; % or == 2, or 3, etc.
% neg = stat.negclusterslabelmat == 1;


timestep      = 0.05; % timestep between time windows for each subplot (in seconds)
sampling_rate = dataFC_LP.fsample; % Data has a temporal resolution of 300 Hz
sample_count  = length(stat.time);
% number of temporal samples in the statistics object
j = [0.55,0.75] ; %  [0:timestep:1]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [466,526]-300 ; % [1:timestep*sampling_rate:sample_count]; % temporal endpoints in M/EEG samples

% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(raweffectFICvsFC.label, stat.label);

cfg = [];
cfg.xlim = j;   % time interval of the subplot
cfg.zlim = [-2.5e-13 2.5e-13];
% If a channel is in a to-be-plotted cluster, then
% the element of pos_int with an index equal to that channel
% number will be set to 1 (otherwise 0).

% Next, check which channels are in the clusters over the
% entire time interval of interest.
pos_int = zeros(numel(raweffectFICvsFC.label),1);
neg_int = zeros(numel(raweffectFICvsFC.label),1);
pos_int(i1) = all(pos(i2, m), 2);
neg_int(i1) = all(neg(i2, m), 2);

cfg.highlight   = 'on';
% Get the index of the to-be-highlighted channel
cfg.highlightchannel = find(pos_int | neg_int);
cfg.comment     = 'no';
%    cfg.commentpos  = 'title';
cfg.layout      = 'CTF151_helmet.mat';
cfg.interactive = 'no';
ft_topoplotER(cfg, raweffectFICvsFC);

savefig('ERPanalysis_Fieldtrip')

%% Make layout for plotting

ft_path = fileparts(which('ft_defaults')) ; % find 

load(fullfile(ft_path,'template','layout','CTF151_helmet.mat'))
clear pos
for i = 1:length(dataFC_LP.label)
    ind = find(strcmp(dataFC_LP.label{i},lay.label)) ; 
    pos(i,:) = lay.pos(ind,:) ; 
end
save('fieldtrip_output','pos','-append')
