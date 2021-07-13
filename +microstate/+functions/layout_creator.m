function layout = layout_creator(varargin) ; 
%UNTITLED2 Summary of this function goes here 
%   Detailed explanation goes here

narginchk(0,2)

switch nargin 
    case 0 
        microstate.functions.layout_creator_gui ; % call the GUI
    case 1
        if ~strcmpi(varargin{1},'gui')
            error("Unexpected input. Useage microstate.functions.layout_creator('gui') to call GUI, or microstate.fucntions.layout_creator(template,labels) to call from the command line.")
        end
        microstate.functions.layout_creator_gui ; % call the GUI
    case 2
        template = varargin{1} ; 
        labels = varargin{2} ; 
        
        % Decide the modality
        switch lower(template)
            case {'eeg','eeg1005','eeg1020','eeg1010'}
                dir = 'layout';
                file = 'EEG1005.mat' ; 
            case 'biosemi'
                dir = 'layout' ; 
                file = 'biosemi256.mat' ; 
            case 'ctf151'
                dir = 'layout' ; 
                file = 'CTF151.mat' ; 
            case 'ctf275' 
                dir = 'layout' ; 
                file = 'CTF275.mat' ; 
            case 'neuromag'
                dir = 'layout' ; 
                file = 'neuromag306.mat' ; 
            case 'aal'
                dir = 'atlas' ; 
                file = 'aal.mat' ; 
            case 'afni_brick0'
                dir = 'atlas' ; 
                file = 'afni_brick0.mat' ; 
            case 'afni_brick1'
                dir = 'atlas' ; 
                file = 'afni_brick1.mat' ; 
            case 'brainnetome'
                dir = 'atlas' ; 
                file = 'brainnetome.mat' ; 
            case 'hcp230'
                dir = 'atlas' ; 
                file = 'hcp230.mat' ; 
            otherwise
                error('First input was not a valid EEG/MEG sensor template or source space atlas')
        end
            
        % make layout structure 
        mspath = microstate.functions.toolbox_path ; 
        if strcmpi(template,'hcp230')
            tmplayout = load(fullfile(mspath,'+external','HCP230',dir,file)) ;
        else
            tmplayout = load(fullfile(mspath,'+external','fieldtrip','template',dir,file)) ;
        end
        if strcmp(dir,'layout')
            tmplayout.type = 'sensor' ; 
        elseif strcmp(dir,'atlas')
            tmplayout.type = 'source' ; 
        end
        
        % make layout
        layout = struct ; 
        switch tmplayout.type
            case 'sensor'
                layout.label = labels ; 
                for i = 1:length(layout.label)
                    ind = find(cellfun(@(c) contains(c,layout.label{i}),tmplayout.label)) ; 
                    if length(ind)>1
                        ind = find(strcmpi(layout.label{i},tmplayout.label)) ; 
                    end
                    if isempty(ind)
                        error('Cannot find electrode labelled %s',layout.label{i}) ; 
                    end
                    layout.pos(i,:) = tmplayout.pos(ind,:) ; 
                end
                layout.outline = tmplayout.outline ; 
            case 'source'
                layout.tissue = labels ; 
                for i = 1:length(layout.tissue)
                    ind = find(strcmp(layout.tissue{i},tmplayout.tissuelabel)) ; 
                    if isempty(ind)
                        error('Cannot find electrode labelled %s',layout.label{i}) ; 
                    end
                    layout.Vertices{i,1} = find(tmplayout.tissue == ind) ; 
                    layout.Centroid = tmplayout.centroid(ind,:) ; 
                end 
        end
        



end

