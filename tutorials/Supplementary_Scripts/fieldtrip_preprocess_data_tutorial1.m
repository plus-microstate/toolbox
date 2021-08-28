%% Fieldtrip preprocessing of data
% NOTE: This script does not contain any tutorials for +microstate toolbox,
% and is just used for pre-processing the data. For the tutorial on
% +microstate, see the script example_script_EEGcanonicalmaps.m
% 
% This script follows the Fieldtrip tutorial for resting-state data
% cleaning [1] to pre-process the example resting-state EEG data set.
% You must have a version of Fieldtrip installed and the raw data
% downloaded [2] to run this script. 
% 
% Full descriptions of the steps and each line of code are given in [1]. 
% This script is just included to get the data in a format ready for
% microstate analysis. 
% 
% [1] https://www.fieldtriptoolbox.org/workshop/madrid2019/tutorial_cleaning/
% [2] ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/workshop/madrid2019/tutorial_cleaning/

clear, clc, close all
rng('default') ; % for reproducibility 

subj = 'sub-22';

% minimal preprocessing
cfg = [] ; 
cfg.dataset =  ['./madrid2019/tutorial_cleaning/single_subject_resting/' subj '_task-rest_run-3_eeg.vhdr'];
cfg.channel    = 'all';
cfg.demean     = 'yes';
cfg.detrend    = 'no';
cfg.reref      = 'yes';
cfg.refchannel = 'all';
cfg.refmethod  = 'avg';
data = ft_preprocessing(cfg);

% prepare electrodes
cd('madrid2019/tutorial_cleaning')
data.elec = prepare_elec_chennu2016(data.label) ; 
cd .. ; cd ..

% load artifact definitions
load('./madrid2019/tutorial_cleaning/sub-22_run-03_eeg_artif')

% load neighbours
load('./madrid2019/tutorial_cleaning/cfg_neighbours', 'neighbours');

% Interpolate bad channels
cfg = [];
cfg.badchannel     = {artif.badchannel};
cfg.method         = 'weighted';
cfg.neighbours     = neighbours;
data_fixed = ft_channelrepair(cfg,data);

% Interpolate bad channels at certain segments
artpadding  = 0.1;
begart      = artif.artfctdef.badchannel.artifact(:,1)-round(artpadding.*data.fsample);
endart      = artif.artfctdef.badchannel.artifact(:,2)+round(artpadding.*data.fsample);
offset      = zeros(size(endart));
% do not go before the start of the recording or the end
begart(begart<1) = 1;
endart(endart>max(data.sampleinfo(:,2))) = max(data.sampleinfo(:,2));

cfg      = [];
cfg.trl  = [begart endart offset];
data_bad = ft_redefinetrial(cfg, data);

% The parameters to detect artifacts are:
proportion  = 0.4; % criterion proportion of bad samples
thresh1     = 3;   % threshold in units of median absolute value over all data
data_fixed  = {};
for k=1:size(data_bad.trial,2)
    w            = ones(size(data_bad.trial{1,k}));
    md           = median(abs(data_bad.trial{1,k}(:)));
    w(find(abs(data_bad.trial{1,k}) > thresh1*md)) = 0;
    iBad         = find(mean(1-w,2)>proportion);
    [val iBad_a] = max(max(abs(data_bad.trial{1,k}.*(1-w)),[],2));
    if isempty(iBad)
        iBad     = find(mean(1-w,2)==max(mean(1-w,2)));
        warning(['decreasing threshold to: ' num2str(max(mean(1-w,2)))]);
    end

    % we use ft_channelrepair to interpolate these short selected artifacts
    cfg = [];
    cfg.badchannel = data_bad.label([iBad;iBad_a]);
    cfg.method     = 'weighted';
    cfg.neighbours = neighbours;
    cfg.trials     = k;
    data_fixed{1,k} = ft_channelrepair(cfg, data_bad);
end

data_fixed = ft_appenddata([], data_fixed{:});

clear data_bad

cfg                               = [];
cfg.artfctdef.minaccepttim        = 0.010;
cfg.artfctdef.reject              = 'partial';
cfg.artfctdef.badchannel.artifact = [begart endart];
data_rejected = ft_rejectartifact(cfg, data);

data = ft_appenddata([], data_rejected, data_fixed);

% clear these variables from memory to avoid confusion later on
clear data_rejected data_fixed

cfg = [];
cfg.trl = [min(data.sampleinfo(:,1)) max(data.sampleinfo(:,2)) 0];
data = ft_redefinetrial(cfg,data);

% Reject visual and muscle artifacts
cfg = [];
cfg.artfctdef.minaccepttim = 0.010;
cfg.artfctdef.reject       = 'partial';
if isfield(artif.artfctdef,'visual')
    cfg.artfctdef.visual.artifact = artif.artfctdef.visual.artifact;
end
if isfield(artif.artfctdef,'muscle')
    cfg.artfctdef.muscle.artifact = artif.artfctdef.muscle.artifact;
end
data = ft_rejectartifact(cfg, data);

% Detrend
cfg             = [];
cfg.channel     = 'all';
cfg.demean      = 'yes';
cfg.polyremoval = 'yes';
cfg.polyorder   = 1; % with cfg.polyorder = 1 is equivalent to cfg.detrend = 'yes'
data = ft_preprocessing(cfg, data);

% Re-reference
cfg            = [];
cfg.channel    = 'all';
cfg.demean     = 'yes';
cfg.reref      = 'yes';
cfg.refchannel = 'all';
cfg.refmethod  = 'avg';
data = ft_preprocessing(cfg, data);

% Visual artifacts for ICA
comp = load(['./madrid2019/tutorial_cleaning/',subj,'_comp']);

cfg = [];
cfg.demean    = 'no';           % This has to be explicitly stated, as the default is to demean.
cfg.unmixing  = comp.unmixing;  % Supply the matrix necessary to 'unmix' the channel-series data into components
cfg.topolabel = comp.topolabel; % Supply the original channel label information
comp = ft_componentanalysis(cfg, data);

cfg = [];
cfg.component = [2,7];
data = ft_rejectcomponent(cfg, comp, data);

% Remove line noise
cfg = [];
cfg.channel    = 'all';
cfg.demean     = 'yes';
cfg.dftfilter  = 'yes';
cfg.dftfreq    = [50];
% cfg.bsfilter  = 'no'; % band-stop method
% cfg.bsfreq    = [48 52];
data = ft_preprocessing(cfg,data);

% Only take 1020 channels
cfg = [] ; 
cfg.channel = 'eeg1020' ; 
data = ft_preprocessing(cfg,data) ; 

% Get 2d layout of data using Fieldtrip
cd('./madrid2019/tutorial_cleaning')
data.elec = prepare_elec_chennu2016(data.label) ;
cd .. ; cd ..
layout = ft_prepare_layout([],data); 
pos = layout.pos(1:end-2,:) ;
clearvars -except data pos

% Save
save('fieldtrip_output','pos','data') ; 
