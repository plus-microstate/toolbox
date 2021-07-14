classdef individual
	%% ------ PROPERTIES -------    
    properties
        data double % Number of samples x number of channels (or ROIs) array of MEG/EEG/source/amplitude envelope time courses
        time (1,:) double {mustBeReal} % 1 x number of samples array containing the time axis
        sample (1,:) double {mustBeReal, mustBeInteger} % 1 x number of samples array specifying which samples are included from the original dataset
        bad_samples (1,:) double {mustBeReal, mustBeInteger} % Array specifying which samples are bad samples (e.g. artifacts)
        modality char % Character array, either 'eeg', 'meg', 'source', or 'ampenv'
        gfp (1,:) double {mustBeReal} % 1 x Number of samples array containing the GFP values
        maps % Number of channels (or ROIs) x number of microstates array containing microstate maps
        label (1,:) double % 1 x number of samples array containing microstate labels of each sample
        stats = struct ; % Structure containing microstate statistics
        process = table([],[],'VariableNames',{'Process','Info'}); % Table listing every process performed on data
        
    end
    
    %% ------ METHODS ------
    
    methods
        
        % Initialization --------------------------------------------------
        function obj = individual(data,modality,time) 
            % Make an individual object with inputs: data, modality, and time
            
            % 3 argument input
            if nargin < 1; data = [] ; end
            if nargin < 2; modality = ''; end
            if nargin < 3; time = [] ; end
                
            if nargin == 1 && isstruct(data)
                obj = obj.import_struct(data) ; 
            elseif nargin == 3
                obj = obj.add_data(data,modality,time) ; 
            elseif nargin == 0
                % do nothing
            else
                error('Useage: data = mstate.Data makes an empty Data object, data = mstate.Data(data,modality,time) inputs data')
            end
        end
        
        function obj = add_data(obj,data,modality,time)
            % Add data to empty microstate.individual object.
            
            % --- add data ---
            if nargin < 2
                error('No time series specified in first argument')
            end
            if nargin < 3
                error('No modality specified in second argument')
            end
            if nargin < 4
                error('Either a time axis or sampling rate must be included as third argument')
            end
            
            % Check is time series
            validateattributes(data,{'numeric'},{'2d','nonsparse','real'})
    
%             % Check is column
%             if size(data,2)>size(data,1)
%                 warning("Each time series should be a column. Possibly the data should be transposed.")
%             end
            
            % Check more than one sampling point
            if size(data,1) == 1
                error('Only a single sample point, microstate analysis not possible')
            end
            
            obj.data = data ;
            obj.sample = 1:size(data,1) ; 
            
            % --- add modality ---
            validmodalities = {'eeg','meg','source','ampenv'}; 
            modality = lower(modality) ; 
            if ~strcmp(modality,validmodalities)
                error('Valid modalities are eeg, meg, source, or ampenv')
            end
            obj.modality = modality ;
            
            % --- add time axis ---
            if length(time) == 1
                fs = time ; 
                time = (1/fs)*(0:(size(data,1)-1)) ; 
            end
            if ~isempty(data) && (length(time) ~= size(data,1))
                error('Time axis must be of equal length to time series')
            end
            obj.time = time ; 
            
            % --- update process ---
            obj = microstate.functions.process_append(obj,'Added data') ; 
            
            % --- calculate some statistics of the data ---
            obj = obj.calculate_gfp() ;    
            
        end
        
        function obj = add_bad_samples(obj,bad_samples)
            % Specify bad samples in microstate.individual object. 
            
            % --- add bad samples ---
            obj.bad_samples = bad_samples ; 
            
            % --- update process ---
            obj = microstate.functions.process_append(obj,'Added bad samples') ; 
            
            % --- calculate some statistics of the data ---
            obj = obj.calculate_gfp() ; 
            
        end
        
        
        %% ----------- IMPORT FROM OTHER TOOLBOXES ------------------------
        function obj = import_struct(obj,data)
            % Imports data saved using individual.export_struct
            f = fieldnames(data) ;
            for i = 1:length(f)
                obj.(f{i}) = data.(f{i}) ; 
            end
            obj = microstate.functions.process_append(obj,'Imported +microstate structure') ; 
        end
        
        function s = export_struct(obj)
	    % Export microstate.individual object's properties to a MATLAB structure.
            f = fieldnames(obj) ;
            s = struct ; 
            for i = 1:length(f)
                s.(f{i}) = obj.(f{i}) ; 
            end
            obj = microstate.functions.process_append(obj,'Exported +microstate structure') ; 
        end
        
        % FIELDTRIP ---
        function obj = import_fieldtrip(obj,data,varargin)
            % Function to import data from Fieldtrip.
	    
            % Check data is a struct
            if ~isstruct(data)
                error('Data type is not a fieldtrip raw, timelock, source, or dip structure. See Fieldtrip functions ft_datatype_raw, ft_datatype_timelock, ft_datatype_source, ft_datatype dip for details.')
            end
            
            % Make options array
            options = struct ; 
            
            % Check if modality was specified
            if length(varargin)~=2*floor(length(varargin)/2) % check even
                error('Additional inputs must be in name-value pairs') ; 
            end
            idx_modality = find(strcmp(varargin,'modality')) ; 
            if ~isempty(idx_modality)
                modality = varargin{idx_modality+1} ; 
                if ~any(strcmp(modality,{'eeg','meg','source','amplitude'}))
                    error('modality must be eeg, meg, source, or amplitude')
                end
                hasmodality = true ; 
                options.modality = modality ; 
            else
                hasmodality = false ; 
            end
            
            % Find the datatype - from ft_datatype in Fieldtrip 20200607
            israw          =  isfield(data, 'label') && isfield(data, 'time') && isa(data.time, 'cell') && isfield(data, 'trial') && isa(data.trial, 'cell') && ~isfield(data,'trialtime');
            istimelock     =  isfield(data, 'label') && isfield(data, 'time') && ~isfield(data, 'freq') && ~isfield(data,'timestamp') && ~isfield(data,'trialtime') && ~(isfield(data, 'trial') && iscell(data.trial)) && ~isfield(data, 'pos') && ~isfield(data,'dip'); 
            issource       = (isfield(data, 'pos') || isfield(data, 'pnt')) && isstruct(data) && numel(data)==1; 
            
            % Check there is only one data type
            datatype = [israw,istimelock,issource] ; 
            if sum(datatype)==0
                error('Data type is not a fieldtrip raw, timelock, or source structure. See Fieldtrip documentation for functions ft_datatype_raw, ft_datatype_timelock, and ft_datatype_source for details.')
            elseif sum(datatype)>1
                error('Cannot determine datatype. See Fieldtrip functions ft_datatype_raw, ft_datatype_timelock, ft_datatype_source, ft_datatype dip for details.')
            end
            
            datatype = find(datatype) ; 
            switch datatype
                case 1 % raw
                    % Check if multiple trials or single trial, and give
                    % warning in the case of multiple trials
                    if length(data.trial)>1
                        warning('%d trials in data, concatenating to a single continuous trial.',length(data.trial)) 
                    end
                    
                    % If the modality isn't specified, try to estimate it
                    % from the data. This is done by looking for sensor
                    % definitions: For EEG data, Fieldtrip calls the sensor
                    % definition data.elec, while for MEG data Fieldtrip
                    % calls the sensors data.grad. 
                    if ~hasmodality
                        if ~isfield(data,'hdr') ; data.hdr.chantype = {'unknown'} ; end
                        iseeg = isfield(data,'elec') || all(strcmp(data.hdr.chantype,'EEG')) ;
                        ismeg = isfield(data,'grad') || all(strcmp(data.hdr.chantype,'MEG')) ; 
                        if (iseeg && ismeg) || (~iseeg && ~ismeg) % neither or both are specified
                            error('Cannot identify data modality, try calling individual.import_fieldtrip(data,modality)')
                        elseif iseeg
                            modality = 'eeg' ; 
                        elseif ismeg
                            modality = 'meg' ; 
                        end
                        options.modality = modality ; 
                    elseif strcmp('modality','source')
                        warning('Input modality was source, but we estimate data to be raw Fieldtrip data')
                    end
                    
                    % Add the data to the object
                    options.labels = data.label ; 
                    obj = microstate.functions.process_append(obj,'Imported fieldtrip raw data structure',options) ; 
                    obj = obj.add_data(cell2mat(data.trial)',modality,1/mean(diff(data.time{1}))) ;
                    
                case 2 % timelock
                    % Check if multiple trials or single trial, and give
                    % warning in the case of multiple trials
                    if ~isfield(data,'trial') ; data.trial = [] ; end
                    if ~isfield(data,'avg') ; data.avg = [] ; end
                    
                    if size(data.trial,3)>1
                        dimord = strsplit(data.dimord,'_') ; 
                        newdim(3) = find(strcmp(dimord,'rpt')) ;
                        warning('%d trials in data, concatenating to a single continuous trial. Use cohort.import_fieldtrip to import trials individually.',size(data.trial,newdim(3)))
                        newdim(1) = find(strcmp(dimord,'chan')) ; 
                        newdim(2) = find(strcmp(dimord,'time')) ; 
                        data.avg = permute(data.trial,newdim) ; 
                        data.time = repmat(data.time,1,size(data.avg,3)) ; 
                        data.avg = reshape(data.avg,size(data.avg,1),size(data.avg,2)*size(data.avg,3)) ; 
                    end
                    
                    % If the modality isn't specified, try to estimate it
                    % from the data. This is done by looking for sensor
                    % definitions: For EEG data, Fieldtrip calls the sensor
                    % definition data.elec, while for MEG data Fieldtrip
                    % calls the sensors data.grad. 
                    if ~hasmodality
                        iseeg = isfield(data,'elec') ;
                        ismeg = isfield(data,'grad') ; 
                        if (iseeg && ismeg) || (~iseeg && ~ismeg) % neither or both are specified
                            error('Cannot identify data modality, try calling individual.import_fieldtrip(data,modality)')
                        elseif iseeg
                            modality = 'eeg' ; 
                        elseif ismeg
                            modality = 'meg' ; 
                        end
                        options.modality = modality ; 
                    elseif strcmp('modality','source')
                        warning('Input modality was source, but we estimate data to be timelock Fieldtrip data')
                    end
                    
                    % Add the data to the object
                    options.labels = data.label ; 
                    obj = microstate.functions.process_append(obj,'Imported fieldtrip timelock data structure',options) ; 
                    obj = obj.add_data(data.avg',modality,data.time') ;                   
                    
                case 3 % source
                    
                    % Modality will be default be source - if amplitude
                    % data is used this needs to be specified
                    if ~hasmodality
                        modality = 'source' ; 
                        options.modality = modality ; 
                    elseif any(strcmp('modality',{'meg','eeg'}))
                        warning('Input modality was MEG or EEG, but we estimate data to be raw Fieldtrip data')
                    end
                    
                    % Include positions of sources in options
                    options.dipole_positions = data.pos(data.inside,:) ; 
                    if isfield(data,'tri') ; options.dipole_faces = data.tri ; end
                    
                    switch data.method
                        case 'rawtrial' % trials are source reconstructed
                            
                            % Get number of trials
                            Ntrial = length(data.trial) ; 
                            if Ntrial>1
                                warning('%d trials in data, concatenating to a single continuous trial. Use cohort.import_fieldtrip to import trials individually.',Ntrial)
                            end
                            
                            % Check that dipole moments are column arrays
                            % and remove dipoles outside of the grid
                            for i = 1:Ntrial
                                data.trial(i).mom = data.trial(i).mom(:) ;
                                data.trial(i).mom = data.trial(i).mom(data.inside) ; 
                            end
                            
                            % concatenate trials
                            mom = [data.trial.mom] ; 
                            data.time = repmat(data.time,1,Ntrial) ; 
                            
                            % Check if scalar or vector source
                            % reconstruction
                            len = cellfun(@(x) size(x,1),data.avg.mom) ; 
                            isscalar = all(len(:) == 1) ; 
                            
                            % Use SVD to reduce to a single dimension
                            if ~isscalar
                                warning('A vector leadfield was used for source reconstruction, meaning source dynamics are >1 dimensional. Using SVD to reduce to 1 dimensional dynamics.')
                                
                                U = cellfun(@(C) svd(C),mom,'UniformOutput',false) ; 
                                mom = cellfun(@(C,U) real(U(:,1)'*C) , mom, U,'UniformOutput',false) ; 
                            end
                            
                            % Add the data to the object
                            obj = microstate.functions.process_append(obj,'Imported fieldtrip source data structure',options) ; 
                            obj = obj.add_data(cell2mat(mom)',modality,data.time') ; 
                            
                            
                        case 'average'
                            
                            % Check that dipole moments are column arrays
                            % and remove dipoles outside of the grid
                            data.avg.mom = data.avg.mom(:) ;
                            data.avg.mom = data.avg.mom(data.inside) ; 
                            
                            % Check if scalar or vector source
                            % reconstruction
                            len = cellfun(@(x) size(x,1),data.avg.mom) ; 
                            isscalar = all(len(:) == 1) ; 
                            
                            % Use SVD to reduce to a single dimension
                            if ~isscalar
                                warning('A vector leadfield was used for source reconstruction, meaning source dynamics are >1 dimensional. Using SVD to reduce to 1 dimensional dynamics.')
                                
                                U = cellfun(@(C) svd(C),data.avg.mom,'UniformOutput',false) ; 
                                data.avg.mom = cellfun(@(C,U) real(U(:,1)'*C) , data.avg.mom, U,'UniformOutput',false) ; 
                            end
                            
                            % Add the data to the object
                            obj = microstate.functions.process_append(obj,'Imported fieldtrip source data structure',options) ; 
                            obj = obj.add_data(cell2mat(data.avg.mom)',modality,data.time') ; 
                                
                        otherwise
                            error('Unrecognised source method, must be either rawtrial or average')   
                    end
            end
            
        end
        
        
        
        
        % SPM ---
        function obj = import_spm_meeg(obj,D,varargin)
            % Function to import SPM M/EEG objects. 
	    
            % Check what type of data we have
            ishdr = isa(D,'struct') ; % hdr
            ismeeg = isa(D,'meeg') ; % meeg object
            isfname = isa(D,'char') ; % filename
            
            % Make options array
            options = struct ; 
            
            % Check if source data should be imported
            idx_source = find(strcmp(varargin,'source')) ; 
            if ~isempty(idx_source)
                issource = true ; 
            else
                issource = false ; 
            end
            varargin(idx_source) = [] ; 
            
            % Check if modality was specified
            if length(varargin)~=2*floor(length(varargin)/2) % check even
                error('Additional inputs must be in name-value pairs') ; 
            end
            idx_modality = find(strcmp(varargin,'modality')) ; 
            if ~isempty(idx_modality)
                modality = varargin{idx_modality+1} ; 
                if ~any(strcmp(modality,{'eeg','meg','source','amplitude','timeOnset','Nsamples','Fsample'}))
                    error('modality must be eeg, meg, source, or amplitude')
                end
                hasmodality = true ; 
                options.modality = modality ; 
            else
                hasmodality = false ; 
            end
            
            
            % Check there is only one data type
            datatype = [ishdr,ismeeg,isfname] ; 
            if sum(datatype)==0
                error('Data type is not a SPM header struct, SPM MEEG object, or file name. See SPM manual for information on SPM MEEG data format.')
            elseif sum(datatype)>1
                error('Cannot determine datatype. See SPM manual for information on SPM MEEG data format.')
            end
            datatype = find(datatype) ; 
            
            % If it is a filename, we will import the header and convert
            % type to hdr
            if isfname
                % Split the file
                [path,name,ext] = fileparts(D) ; 
                
                % Check that either a .mat, .dat was supplied, or no
                % extension was given
                if ~isempty(ext) && ~any(strcmp(ext,{'.mat','.dat'}))
                    error('Extension %s not recognised, please specify filename for either a .mat or .dat SPM file',ext) ; 
                end
                
                % Check header file exists
                hdrfile = fullfile(path,[name,'.mat']) ; 
                if ~exist(hdrfile,'file')
                    error('Cannot find SPM header file %s',hdrfile)
                end
                       
                % Read and check header file
                try
                    D = load(hdrfile).D ;
                catch
                    error('Cannot read the SPM structure D from header file %s',hdrfile)
                end
                        
                % Check the header file has the relevant fields
                if ~isstruct(D)
                    error('SPM header file %s should contain an SPM header structure called "D"',hdrfile) ; 
                end
                if ~all(isfield(D,{'channels','data','fname','path','sensors'}))
                    error('SPM header file %s should contain an SPM header structure called "D"',hdrfile) ; 
                end
                
                % Change datatype to hdr
                datatype = 1 ; 
            end
                        
            
            % Switch between data types
            switch datatype
                
                case 1 % header
                    
                    if ~issource
                        % IMPORT SENSOR DATA
                        
                        % Check data file exists
                        datafile = D.fname ; datafile(end-3:end) = '.dat' ; 
                        datafile = fullfile(D.path,datafile) ; 
                        if ~exist(datafile,'file')
                            error('Cannot find SPM data file %s',datafile)
                        end

                        % Read in the data
                        fid = fopen(datafile) ; 
                        data = fread(fid,size(D.data),'float32') ; 
                        fclose(fid) ; 

                        % Make time
                        time = D.timeOnset + (1/D.Fsample)*((1:D.Nsamples)-1) ; 

                        % Get the modality
                        if ~hasmodality
                            chantype = {D.channels.type} ;
                            options.labels = {D.channels(:).label} ; 
                            if all(strcmp(chantype,'EEG')) ; 
                                modality = 'eeg' ; 
                            elseif all(strcmp(chantype,'MEG')) ; 
                                modality = 'meg' ; 
                            else
                                error('Cannot identify data modality, try calling individual.import_spm_meeg(D,modality)')
                            end
                        end
                        
                    else 
                        % IMPORT SOURCE DATA
                        
                        % Check there has been some source data
                        if ~isfield(D,'other') ; 
                            error('No source data saved in the SPM header file')
                        end
                        if ~isfield(D.other,'inv')
                            error('No source data saved in the SPM header file')
                        end
                        
                        % Check if sourceval was specified
                        idx_sourceval = find(strcmp(varargin,'sourceval')) ; 
                        if ~isempty(idx_sourceval)
                            sourceval = varargin{idx_sourceval+1} ; 
                            validateattributes(sourceval,{'numeric'},{'integer','real','nonnegative'},'import_spm','source')
                        else
                            sourceval = length(D.other.inv) ; 
                        end
                        
                        % Get the inverse
                        inv = D.other.inv{sourceval}.inverse ; 
                        
                        % Get number of trials
                        Ntrial = length(inv.trials) ; 
                        if Ntrial>1
                            warning('%d trials in data, concatenating to a single continuous trial.',Ntrial)
                        end
                          
                        % Reconstruct the source data
                        data = [] ; 
                        for i = 1:Ntrial
                            data = inv.J{i}*inv.T' ; 
                        end
                        
                        % Make time
                        time = D.timeOnset + (1/D.Fsample)*((1:D.Nsamples)-1) ; 
                        time = time(1:size(data,2)) ; 

                        % Get the modality
                        if ~hasmodality
                            modality = 'source' ; 
                        end
                        
                        % Source positions
                        options.dipole_position = D.other.inv{sourceval}.mesh.tess_mni.vert ;
                        options.dipole_faces = D.other.inv{sourceval}.mesh.tess_mni.face ;
                        
                    end
                        
                case 2 % meeg object
                    
                    if ~issource
                        % IMPORT SENSOR DATA
                        
                        % Check spm is on the path
                        spm_path = which('spm') ; 
                        addpath(fileparts(spm_path)) ; % just in case not on path

                        % The following should work without needing to
                        % initialize SPM, just having it on the path is
                        % sufficient
                        data = D(:,:,:) ;
                        if size(D,3)>1
                            warning('%d trials in data, concatenating to a single continuous trial.',size(D,3))
                            data = reshape(data,size(D,1),size(D,2)*size(D,3)) ; 
                        end
                        time = D.time ; 

                        % Get the modality
                        chantype = D.chantype ; 
                        options.labels = D.chanlabels ;
                        if all(strcmp(chantype,'EEG')) ; 
                            modality = 'eeg' ; 
                        elseif all(strcmp(chantype,'MEG')) ; 
                            modality = 'meg' ; 
                        else
                            error('Cannot identify data modality, try calling individual.import_spm_meeg(D,modality)')
                        end
                        
                    else
                        % IMPORT SOURCE DATA
                        
                        % Check spm is on the path
                        spm_path = which('spm') ; 
                        addpath(fileparts(spm_path)) ; % just in case not on path
                        
                        
                        % Check there has been some source data
                        if ~isfield(D,'inv')
                            error('No source data saved in the SPM header file')
                        end
                        
                        % Check if sourceval was specified
                        idx_sourceval = find(strcmp(varargin,'sourceval')) ; 
                        if ~isempty(idx_sourceval)
                            sourceval = varargin{idx_sourceval+1} ; 
                            validateattributes(sourceval,{'numeric'},{'integer','real','nonnegative'},'import_spm','source')
                        else
                            sourceval = length(D.inv) ; 
                        end
                        
                        % Get the inverse
                        inv = D.inv{sourceval}.inverse ; 
                        
                        % Get number of trials
                        Ntrial = length(inv.trials) ; 
                        if Ntrial>1
                            warning('%d trials in data, concatenating to a single continuous trial.',Ntrial)
                        end
                          
                        % Reconstruct the source data
                        data = [] ; 
                        for i = 1:Ntrial
                            data = inv.J{i}*inv.T' ; 
                        end
                        
                        % Make time
                        time = D.time(1:size(data,2)) ; 

                        % Get the modality
                        if ~hasmodality
                            modality = 'source' ; 
                        end
                        
                        % Source positions
                        options.dipole_position = D.inv{sourceval}.mesh.tess_mni.vert ;
                        options.dipole_faces = D.inv{sourceval}.mesh.tess_mni.face ;
                    end
                             
            end
            options.modality = modality ;
            if issource
                options.sensorsource = 'source' ; 
            else
                options.sensorsource = 'sensor' ; 
            end
            
            
            % Add the data to the object
            obj = microstate.functions.process_append(obj,'Imported SPM meeg data',options) ; 
            obj = obj.add_data(data',modality,time) ; 
            
        end
        
        % ---- SPM IMAGES ----
        function obj = import_spm_nifti(obj,filename,fsample,zthresh) ; 
	    % Function to import NiFTi images. 
            
            if nargin < 3
                error('Sample rate must be provided to read SPM/NiFTi images. Useage individual.import_spm_nifti(filename,sample_rate)')
            end
            validateattributes(fsample,{'numeric'},{'scalar','real','positive'},'import_spm_nifti','fsample')
            validateattributes(filename,{'char'},{},'import_spm_nifti','fname')
            
            % Initialize output
            options = struct ; 
            
            % Read hdr and img
            hdr = niftiinfo(filename) ; 
            img = niftiread(filename) ; 
            
            % Get threshold
            if nargin<4
                zthresh = 0 ; 
            end
            
            % Get indices of active dipoles
            pow = mean(abs(img),4) ; 
            ind = find(pow>0) ; 
            
            % z-score
            stdz = vecnorm(pow(ind))/sqrt(length(ind)-1) ; 
            z = pow/stdz ; clear pow 
            ind = find(z>zthresh) ; 
            [x,y,z] = ind2sub(size(z),ind) ; 
            
            % Get dipole positions
            P = hdr.Transform.T*[x';y';z';ones(size(x'))] ; 
            options.dipole_positions = P(1:3,:)' ; 
            
            % Make data
            data = zeros(size(img,4),length(ind)) ; 
            for i = 1:length(ind)
                data(:,i) = img(x(i),y(i),z(i),:) ; 
            end
                
            % Make modality 
            modality = 'source' ; % always source
            
            % Add the data to the object
            obj = microstate.functions.process_append(obj,'Imported SPM/NiFTi image',options) ; 
            obj = obj.add_data(data,modality,fsample) ;
        end
            
        
        
        % ---- EEGLAB ----
        function obj = import_eeglab(obj,EEG,varargin) ; 
            % Function to import data from EEGLab.
            
            idxmodality = find(strcmp(varargin,'modality')) ; 
            if ~isempty(idxmodality)
                modality = varargin{idxmodality+1} ; 
            else
                modality = 'eeg'; 
            end
            fsample = EEG.srate;
            
            % Check if multiple trials or single trial, and give
            % warning in the case of multiple trials
            if EEG.trials>1
                warning('%d trials in data, concatenating to a single continuous trial.',EEG.trials)
            end
                    
            for index = 1:EEG.trials
                trial{index}  = EEG.data(:,:,index);
            end
            
            obj = obj.add_data(cell2mat(trial)',modality,fsample) ; 
        end
        
        
        % ---- LORETA (SOURCE ONLY) ----
        function obj = import_loreta(obj,filename,fsample)
            % Function to import source reconstructions from LORETA or sLORETA/eLORETA software.
            % 
            % Much of this function is modified from Fieldtrip function
            % loreta2fieldtrip
            
            % Check inputs are correct
            if nargin < 3
                error('Sample rate must be provided to read LORETA/sLORETA files. Useage individual.import_loreta(filename,sample_rate)')
            end
            validateattributes(fsample,{'numeric'},{'scalar','real','positive'},'import_loreta','fsample')
            validateattributes(filename,{'char'},{},'import_loreta','fname')
            
            % Get the extension of the file name
            [~,~,ext] = fileparts(filename) ; 
            
            % Import different file types
            switch ext
                case '.lorb' 
                    % IMPORT LORETA BINARY FILE
                    
                    % Number of voxels
                    voxnumber = 2394;
                    
                    % Read in positions of dipoles
                    dr = pwd ; 
                    cd(fullfile(microstate.functions.toolbox_path,'+external','LORETA')) ; 
                    options.dipole_position = load('LORETA_dipole_positions','pos') ; 
                    cd(dr)
                    
                    % Mark whether txt or binary
                    ftype = 'bin' ; 
                    
                    % Update options
                    options.software = 'LORETA' ;  
                    
                case {'.txt','.slor'} 
                    % IMPORT SLORETA FILES
                    
                    % Number of voxels
                    voxnumber = 6239;
                    
                    % Read in positions of dipoles
                    dr = pwd ; 
                    cd(fullfile(microstate.functions.toolbox_path,'+external','LORETA')) ; 
                    options.dipole_position = load('sLORETA_dipole_positions','pos') ; 
                    cd(dr)
                    
                    % Mark whether txt or binary
                    if strcmp(ext,'.slor')
                        ftype = 'bin' ;
                    elseif strcmp(ext,'.txt')
                        ftype = 'txt' ; 
                    end
                    
                    % Update options
                    options.software = 'sLORETA/eLORETA' ;  
                 
                    
                % --- Here we deal with the case of 3d current density
                % vectors ---
                case '.lor3'
                    error('LORETA-xyz files (.lor3) are not supported. Compute 1d current density vectors (.lor) using the window EEG/ERP->LORETA, *not* the window EEG/ERP->LORETA-xyz.')
                case '.slor3'
                    error('sLORETA-xyz files (.slor3) are not supported. Compute 1d current density vectors (.slor) using the option "Compute sLORETA" in the "EEG/ERP to sLORETA" window, *not* the option "Compute sLORETA-xyz".')
                % ---    
                    
                
            end % end switch 
            
            % Read in the data
            switch ftype
                case 'bin'
                    
                    % work with binary file
                    fid = fopen(filename,'r', 'ieee-le');
                    % determine the length of the file
                    fseek(fid, 0, 'eof');
                    filesize = ftell(fid);
                    Ntime = filesize/voxnumber/4;
                    % read binary file
                    fseek(fid, 0, 'bof');
                    data = fread(fid, [voxnumber Ntime], 'float=>double');
                    fclose(fid);
                    data = data' ; 
                    
                case 'txt'
                    
                    % read with textfile
                    data = readmatrix(filename); % changed from dlmread in fieldtrip
                    if size(data,1) == voxnumber
                        data = data' ; % turn to column vector
                    elseif size(data,2) == voxnumber
                        % do nothing, this is fine
                    elseif any(size(data) == 3*voxnumber)
                    	error('sLORETA-xyz files are not supported. Compute 1d current density vectors using the option "Compute sLORETA" in the "EEG/ERP to sLORETA" window, *not* the option "Compute sLORETA-xyz".')
                    else
                        error('Number of voxels in .txt file (%d-by-%d) does not match the expected number (%d)',size(data,1),size(data,2),voxnumber) ; 
                    end
                    
            end
            
            % Make time axis
            time = (0:(size(data,1)-1))*(1/fsample) ; 
            
            % Make modality 
            modality = 'source' ; % always source
            
            % Add the data to the object
            obj = microstate.functions.process_append(obj,'Imported LORETA/sLORETA source reconstruction',options) ; 
            obj = obj.add_data(data,modality,time) ;
        end
        
        % ---- BRAINSTORM ----
        function obj = import_brainstorm(obj,bst_data,varargin) ; % bst_source,bst_scout)
            % Function to import brainstorm source data. 
            % 
            % Useage: 
            % EEG/MEG: 
            % obj = obj.import_brainstorm(bst_data,modality)
            % 
            % SOURCE: 
            % obj = obj.import_brainstorm(bst_data,bst_source)
            % 
            % PARCELLATED SOURCE: 
            % obj = obj.import_brainstorm(bst_data,bst_source,bst_scouts)
            % 
            % Inputs: 
            % - bst_data: Right click on the "Data" file containing the
            % EEG/MEG data -> file -> export to matlab. 
            % - modality: 'eeg', 'meg', or 'source'
            % - bst_source: Following source reconstruction, right click on
            % the "Results kernel" file in Brainstorm -> file-> export to
            % Matlab (required for source space)
            % - bst_scouts: In the "Scouts" tab, select the "Scouts"
            % dropdown and export to Matlab (required for source space)
            
            % Check inputs
            if isempty(varargin)
                error('At least two inputs must be specified')
            end
            if length(varargin) > 2
                error('Up to 3 inputs expected')
            end
            if length(varargin) == 2
                if isempty(varargin{2})
                    varargin(2) = [] ; 
                end
            end
            
            % Find out what the extra inputs were
            ismeeg = false ; issource = false ; isparcelsource = false ; 
            bst_source = [] ; bst_scout = [] ; 
            
            if length(varargin) == 2 % Check if parcelsource
                isparcelsource = true ; 
                
                % Check the 2nd input is valid
                bst_scout = varargin{2} ; 
                validateattributes(bst_scout,{'struct'},{},'import_brainstorm','bst_scout',3)
                if ~isfield(bst_scout,'Vertices')
                    error('Input bst_scout must have a field "Vertices"')
                end
                
                % Check the 1st input is valid
                bst_source = varargin{1} ; 
                validateattributes(bst_source,{'struct'},{},'import_brainstorm','bst_source',2)
                if ~isfield(bst_source,'ImagingKernel')
                    error('Input bst_source must have a field "ImagingKernel"')
                end
                
                % Make modality
                modality = 'source' ; 
                    
            elseif isstruct(varargin{1}) % Check if source
                issource = true ; 
                
                % Check the 1st input is valid
                bst_source = varargin{1} ; 
%                 validateattributes(bst_source,{'struct'},{},'import_brainstorm','bst_source',2)
                if ~isfield(bst_source,'ImagingKernel')
                    error('Input bst_source must have a field "ImagingKernel"')
                end
                
                % Make modality
                modality = 'source' ; 
                
            elseif ischar(varargin{1}) % Check if EEG/MEG
                ismeeg = true ; 
                
                % Get the modality
                modality = varargin{1} ; 
                
            else  
                % Output error
                validateattributes(varargin{1},{'struct','char'},{},2)
            end
                
                
            % Get datatype
            datatype = find([ismeeg,issource,isparcelsource]) ; 

            switch datatype
                case 1 % M/EEG data
                    
                    % Add the data to the object
                    obj = microstate.functions.process_append(obj,'Imported Brainstorm sensor data') ; 
                    obj = obj.add_data(bst_data.F',modality,bst_data.Time) ;
                    
                case 2 % Source data
                    
                    % Add the data to the object
                    obj = microstate.functions.process_append(obj,'Imported Brainstorm source data') ; 
                    obj = obj.add_data((bst_source.ImagingKernel*bst_data.F)',modality,bst_data.Time) ;
                    
                case 3 % Parcellated Source data
                    
                    % Parcellate
                    warning('off','stats:pca:ColRankDefX') ; % often rank deficient
                    for i = 1:length(bst_scout)
                        
                        K = bst_source.ImagingKernel(bst_scout(i).Vertices,:) ; 
                        s = K*bst_data.F ; 
                        
                        % Get mean variance
                        mv = mean(var(s')) ; 
                        
                        % PCA
                        [~,s,~,~,varexpl] = pca(s','NumComponents',1) ; 
                        
                        % Rescale to have variance equal to variance
                        % explained by first PC multiplied by average
                        % variance of data
                        data(:,i) = sqrt(varexpl(1)*mv/var(s))*s ; 
                    end
                    warning('off','stats:pca:ColRankDefX') ; % often rank deficient
                         
                    % Add the data to the object
                    obj = microstate.functions.process_append(obj,'Imported Brainstorm source data') ; 
                    obj = obj.add_data(data,modality,bst_data.Time) ;
              
            end
        end
       
        
            
    end

end
        
