classdef cohort
	%% ------ PROPERTIES -------    
    properties
        individual = [] ; % Array of microstate.individual objects representing each individual in the cohort
        condition (:,1) double {mustBeReal, mustBeFinite} = [] ; % Array of integers labelling the condition corresponding to each individual
        conditionlabels = cell(0) ; % Cell array giving each integer value in cohort.condition a label
        globalmaps double % Array of global microstate maps
        stats = struct ; % Structure containing group level microstate statistics
        process = table([],[],'VariableNames',{'Process','Info'}); % List of all processes performed on the data
        
    end
    
    %% ------ METHODS ------
    
    methods
        
        % Initialization --------------------------------------------------
        function obj = cohort(individuals,conditions,datastorage) 
            % Make a cohort object with inputs: individuals, conditions, and datastorage
            
            % 3 argument input
            if nargin < 1; individuals = [] ; end
            if nargin < 2; conditions = []; end
            if nargin < 3; datastorage = 0; end
                
            % add to the object
            obj = obj.add_individuals(individuals,conditions,datastorage) ; 
            
        end
        
        function obj = add_individuals(obj,individuals,conditions,datastorage)
            % Add an individual to the microstate.cohort object.
            
            if nargin < 2
                individuals = [] ; 
            end
            if nargin < 3
                conditions = [] ; 
            end
            if nargin < 4
                datastorage = 0 ; 
            end
            
            % Check individuals input
            if isempty(individuals)
                return % return if no individuals are entered
            end
            if isa(individuals,'microstate.cohort')
                conditions = individuals.conditionlabels(individuals.condition) ; 
                individuals = individuals.individual ; 
            end
            if isstruct(individuals) % if structure array, convert to array of individual objects
                tmpind = individuals(:) ; 
                individuals = [] ; 
                for i = 1:numel(tmpind)
                    individuals(i,1) = microstate.individual(tmpind(i)) ; 
                end
            end
            if ~isa(individuals,'microstate.individual') % check is microstate.individual object
                error('Input individuals should be an array of microstate.individual objects')
            end
            individuals = individuals(:) ; 
            
            % Check datastorage input
            if ~isnumeric(datastorage) || length(datastorage) ~= 1
                error('Input datastorage must be either 0, a float 0-1, or an integer >= 1')
            end
            if datastorage == 0 % do not store data in the cohort structure
                for i = 1:length(individuals)
                    individuals(i).data = [] ; 
                    individuals(i).sample = [] ; 
                end
            elseif (datastorage > 0)&&(datastorage < 1) % store this percentage of the data randomly sampled
                for i = 1:length(individuals)
                    % can't store data if there isn't any
                    if isempty(individuals(i).data)
                        continue
                    end
                    
                    % set gfp to have unit
                    individuals(i) = individuals(i).calculate_gfp() ; 
                    gfp = individuals(i).gfp ; 
                    gfp = gfp-min(gfp) ; 
                    gfp = gfp/nansum(gfp) ; 
                    
                    % check there is a gfp, if not, set it to uniform
                    if isempty(gfp)
                        gfp = ones(1,size(individuals(i).data,1))/size(individuals(i).data,1) ; 
                    end
                    
                    % sample from this distribution
                    numsample = round(datastorage*length(gfp)) ; 
                    y = randsample(length(gfp),numsample,true,gfp) ; 
                    
                    % this has to be done with replacement, but we don't
                    % want replacement, so we need to loop until it is done
                    % without replacement
                    while length(y) ~= length(unique(y))
                        y = unique(y,'stable') ; 
                        
                        % only keep values which haven't been selected
                        remain = 1:length(gfp) ; 
                        remain(y) = [] ; 
                        gfpremain = gfp ; 
                        gfpremain(y) = [] ; 
                        
                        % resample new ones
                        ynew = randsample(length(gfpremain),numsample-length(y),true,gfpremain) ; 
                        ynew = remain(ynew); 
                        
                        % append
                        y = [y(:);ynew(:)] ;  
                    end
                    
                    % keep only these datapoints
                    individuals(i).data = individuals(i).data(y,:) ; 
                    individuals(i).sample = y ; 
                    
                end
                
            elseif datastorage == 1
                
                % do nothing, keep all samples
                
            elseif datastorage > 1 && (round(datastorage) == datastorage) % integer greater than 1
                
                for i = 1:length(individuals)
                    % can't store data if there isn't any
                    if isempty(individuals(i).data)
                        continue
                    end
                    
                    % check datastorage < number of samples
                    if datastorage > size(individuals(i).data,1)
                        if length(individuals) == 1
                            warning('Fewer samples in data than number of samples selected in datastorage, keeping all samples')
                        else
                            warning('Individual %d - Fewer samples in data than number of samples selected in datastorage, keeping all samples')
                        end
                        continue % do nothing, keep all samples
                    end
                    
                    % ensure there is a gfp
                    individuals(i) = individuals(i).calculate_gfp() ; 
                    gfp = individuals(i).gfp ; 
                    
                    % check there is a gfp
                    samplepeaks = true ; % assume sampling peaks unless otherwise
                    if isempty(gfp)
                        if length(individuals) == 1
                            warning('No GFP available (perhaps individual.modality was not specified?), randomly sampling %d points',datastorage) ; 
                        else 
                            warning('Individual %d - No GFP available (perhaps individual.modality was not specified?), randomly sampling %d points',i,datastorage) ; 
                        end
                        gfp = ones(1,size(individuals(i).data,1))/size(individuals(i).data,1) ; 
                        samplepeaks = false ; 
                    end
                    
                    if samplepeaks
                        % find gfp peaks
                        [~,locs] = findpeaks(gfp) ; 
                        
                        % randomly sample peaks
                        if length(locs) >= datastorage
                            y = randperm(length(locs),datastorage) ; 
                            y = locs(y) ; 
                        else
                            y = locs ; 
                            
                            if length(individuals) == 1
                                warning('Too few GFP peaks available (%d available, %d requested), randomly sampling %d points with probability weighted by GFP',length(locs),datastorage,datastorage-length(locs)) ; 
                            else 
                                warning('Individual %d - Too few GFP peaks available (%d available, %d requested), randomly sampling %d points with probability weighted by GFP',i,length(locs),datastorage,datastorage-length(locs)) ; 
                            end
                            
                            samplepeaks = false ; 
                        end
                    else
                        y = [] ; 
                    end
                    
                    if ~samplepeaks % this is a separate if statement and not an else as it might be set to false in the previous statement
                            
                        % set gfp to have unit
                        if var(gfp)>0 % don't do this step if gfp is uniform
                            gfp = individuals(i).gfp-min(individuals(i).gfp) ; 
                        end
                        gfp = individuals(i).gfp/sum(individuals(i).gfp) ; 

                        % sample from this distribution
                        numsample = datastorage-length(y) ; 
                        y = [y(:);randsample(length(gfp),numsample,true,gfp)] ; 

                        % this has to be done with replacement, but we don't
                        % want replacement, so we need to loop until it is done
                        % without replacement
                        while length(y) ~= length(unique(y))
                            y = unique(y,'stable') ; 

                            % only keep values which haven't been selected
                            remain = 1:length(gfp) ; 
                            remain(y) = [] ; 
                            gfpremain = gfp ; 
                            gfpremain(y) = [] ; 

                            % resample new ones
                            ynew = randsample(length(gfpremain),datastorage-length(y),true,gfpremain) ; 
                            ynew = remain(ynew); 

                            % append
                            y = [y(:);ynew(:)] ;  
                        end
                        
                    end

                    % keep only these datapoints
                    individuals(i).data = individuals(i).data(y,:) ; 
                    individuals(i).sample = y ; 
                    
                end
                
            else % error in datastorage
                    
                error('Input datastorage must be either 0, a float 0-1, or an integer >= 1')
                
            end % finished with datastorage input
                
                    
            % Add individuals
            obj.individual = [obj.individual ; individuals] ; 
            
            % Check type of conditions is cell
            if isempty(conditions) % deal with the case of no conditions supplied
                conditions = ones(numel(individuals),1) ; 
            end
            if isnumeric(conditions)
                conditions = mat2cell(conditions(:),ones(numel(conditions),1),1) ;
                conditions = cellfun(@(C) sprintf('Condition %d',C),conditions,'UniformOutput',false) ; 
            elseif isa(conditions,'char') || isa(conditions,'string')
                conditions = {conditions} ; 
            elseif iscell(conditions)
                % do nothing
            else
                error('Cannot identify type of conditions')
            end
            
            % Check all condition labels are char (required for call to
            % function unique later)
            for i = 1:length(conditions)
                conditions{i} = num2str(conditions{i}) ; % converts doubles, strings, char to char
            end
            
            % Relabel conditions
            labels = unique([obj.conditionlabels ; conditions(:)],'stable') ; % using stable ensures that the order of the conditions doesn't change
            tmpconditions = nan(length(conditions),1) ; % overwrite as a double array
            for i = 1:length(labels)
                ind = strcmp(labels{i},conditions) ; 
                tmpconditions(ind) = i ; 
            end
            
            % Add to object
            obj.conditionlabels = labels ; 
            nbefore = length(obj.condition) ; nafter = length(conditions) ; 
            obj.condition = [obj.condition ; tmpconditions] ; 
            
            % update output
            options = struct ; 
            options.datastorage = datastorage ; 
            if exist('locs','var')
                options.numgfppeaks = length(locs) ; 
            end
            obj = microstate.functions.process_append(obj,sprintf('Added individuals %d-%d',nbefore+1,nafter),options) ; 
            
        end % finished with add_individuals
        
        function obj = delete_individuals(obj,ind)
           
            ind = sort(ind,'descend') ; 
            for i = 1:length(ind)
                obj.individual(ind(i)) = [] ; 
                obj.condition(ind(i)) = [] ; 
            end
        end
        
        
        % ---------- INDIVIDUAL LEVEL FUNCTIONS ---------------------------
        function obj = ind_calculate_gfp(obj)
            % Run calculate_gfp for each individual
            for i = 1:length(obj.individual)
                obj.individual(i) = obj.individual(i).calculate_gfp ; 
            end
        end
        
        function obj = ind_preprocess(obj)
            % Run preprocess for each individual
            for i = 1:length(obj.individual)
                obj.individual(i) = obj.individual(i).preprocess ; 
            end
        end
        
        function obj = ind_preprocess_ampenv(obj,varargin)
            % Run preprocess_ampenv for each individual
            for i = 1:length(obj.individual)
                obj.individual(i) = obj.individual(i).preprocess_ampenv(varargin) ; 
            end
        end
        
        function obj = ind_preprocess_orthogonalize(obj,varargin)
            % Run preprocess_ampenv for each individual
            for i = 1:length(obj.individual)
                obj.individual(i) = obj.individual(i).preprocess_orthogonalize(varargin) ; 
            end
        end
        
        function obj = ind_preprocess_filter(obj,flow,fhigh,varargin)
            % Run preprocess_filter for each individual
            switch nargin
                case 1
                    for i = 1:length(obj.individual)
                        obj.individual(i) = obj.individual(i).preprocess_filter ; 
                    end
                case 2
                    for i = 1:length(obj.individual)
                        obj.individual(i) = obj.individual(i).preprocess_filter(flow) ; 
                    end
                case 3
                    for i = 1:length(obj.individual)
                        obj.individual(i) = obj.individual(i).preprocess_filter(flow,fhigh) ; 
                    end
                case 4
                    for i = 1:length(obj.individual)
                        obj.individual(i) = obj.individual(i).preprocess_filter(flow,fhigh,varargin) ; 
                    end
            end
        end
        
        function obj = ind_preprocess_rereference(obj)
            % Run preprocess_rereference for each individual
            for i = 1:length(obj.individual)
                obj.individual(i) = obj.individual(i).preprocess_rereference ; 
            end
        end
        
        function obj = ind_preprocess_resample(obj,fsample_new)
            % Run preprocess_resample for each individual
            switch nargin
                case 1
                    for i = 1:length(obj.individual)
                        obj.individual(i) = obj.individual(i).preprocess_resample ; 
                    end
                case 2
                    for i = 1:length(obj.individual)
                        obj.individual(i) = obj.individual(i).preprocess_resample(fsample_new) ; 
                    end
            end
        end
        
        function obj = ind_stats_all(obj)
           % Run stats_all for each individual
            msg = [] ;
            for i = 1:length(obj.individual)
                fprintf(repmat('\b',1,length(msg))) ; 
                msg = sprintf('Calculating stats for individual %d of %d',i,length(obj.individual)) ; 
                fprintf(msg)
                obj.individual(i) = obj.individual(i).stats_all ; 
            end
        end
        
        function obj = ind_stats_autoinformation(obj,maxlag)
            % Run stats_autoinformation for each individual
            switch nargin 
                case 1
                    for i = 1:length(obj.individual)
                        obj.individual(i) = obj.individual(i).stats_autoinformation ; 
                    end
                case 2
                    for i = 1:length(obj.individual)
                        obj.individual(i) = obj.individual(i).stats_autoinformation(maxlag) ; 
                    end
            end
        end
        
        function obj = ind_stats_complexity(obj,Ntransitions)
            % Run stats_complexity for each individual
            switch nargin 
                case 1
                    for i = 1:length(obj.individual)
                        obj.individual(i) = obj.individual(i).stats_complexity ; 
                    end
                case 2
                    for i = 1:length(obj.individual)
                        obj.individual(i) = obj.individual(i).stats_complexity(Ntransitions) ; 
                    end
            end
        end
        
        function obj = ind_stats_coverage(obj)
            % Run stats_coverage for each individual
            for i = 1:length(obj.individual)
                obj.individual(i) = obj.individual(i).stats_coverage ; 
            end
        end
        
        function obj = ind_stats_duration(obj)
            % Run stats_duration for each individual
            for i = 1:length(obj.individual)
                obj.individual(i) = obj.individual(i).stats_duration ; 
            end
        end
        
        function obj = ind_stats_gfp_peaksfreq(obj)
            % Run stats_gfp_peaksfrq for each individual
            for i = 1:length(obj.individual)
                obj.individual(i) = obj.individual(i).stats_gfp_peaksfreq ; 
            end
        end
        
        function obj = ind_stats_hurst(obj)
            % Run stats_hurst for each individual
            for i = 1:length(obj.individual)
                obj.individual(i) = obj.individual(i).stats_hurst ; 
            end
        end
        
        function obj = ind_stats_markov(obj)
            % Run stats_markov for each individual
            for i = 1:length(obj.individual)
                obj.individual(i) = obj.individual(i).stats_markov ; 
            end
        end
        
        function obj = ind_stats_occurrence(obj)
            % Run stats_occurrence for each individual
            for i = 1:length(obj.individual)
                obj.individual(i) = obj.individual(i).stats_occurrence; 
            end
        end
        
        function obj = ind_stats_syntax(obj)
            % Run stats_syntax for each individual
            for i = 1:length(obj.individual)
                obj.individual(i) = obj.individual(i).stats_syntax ; 
            end
        end
        
        function nets = ind_networks_wpli(obj,frq,keepstates,epochlength)
            % Defaults
            if nargin < 4
                epochlength = 1280; 
            end
            if nargin < 3
                keepstates = false ; 
            end
            if nargin < 2
                frq = [] ;
            end
            % Run networks for each individual
            for i = 1:length(obj.individual)
                nets{i} = obj.individual(i).networks_wpli(frq,keepstates,epochlength) ; 
            end
        end
            
        
        % -------- COLLECT DATA FROM INDIVIDUALS TO GROUP LEVEL ----------
        function obj = cohort_stats(obj)
            % collect microstate statistics into a single structure
            
            % Get the number of conditions
            numconditions = max(unique(obj.condition)) ; 

               
            % Loop over conditions
            for cond = 1:numconditions
                ind = find(obj.condition == cond) ; 
                statscond = struct ; 

                for i = 1:length(ind)
                    % Append the stats
%                     try
                        statscond = microstate.functions.structappend(statscond,obj.individual(ind(i)).stats) ; 
%                     catch
%                         error('Could not append stats for individual %d. Ensure all individuals have the same statistics calculated.',ind(i))
%                     end

                    % Check all stats are present for all participants
                    checkstructlength(statscond,i) ;
                end

                statsi(cond) = statscond ; 
            end  
            obj.stats = statsi ; 
            
            obj = microstate.functions.process_append(obj,'Collected group statistics') ; 
            
            function checkstructlength(x,i)
                % check length of all fields in structure are i
                f = fieldnames(x) ; 
                for j = 1:length(f)
                    if ~isstruct(x.(f{j}))
                        if size(x.(f{j}),1) ~= i
                            error('Some stats are not present for all individuals (individual %d: %s). Ensure all individuals have the same statistics calculated.',i,f{j})
                        end
                    else
                        checkstructlength(x.(f{j}),i) ; 
                    end
                end
            end
            
            
        end
        
   
    end
    
    
end
        