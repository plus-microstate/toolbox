function obj = simulate_seq_randomwalk(obj,varargin)
% Simulate a random-walk decision-tree microstate sequence

    % check inputs
    options = microstate.functions.make_options(varargin) ; 
    
    % default options
    defaults = {'Nsample',[] ; 
                'fsample',[] ; 
                'Nstates',4 ; 
                'Hurnst',0.71 ; 
                'meandur',[] ; % will be set to 50 later if filterord wasn't specified
                'filterord',[]} ; 
    options = microstate.functions.add_options(options,defaults) ; clear defaults
    
    
    % check object and inputs match
    hastime = ~isempty(obj.time) ; 
    hasNsample = ~isempty(options.Nsample) ; 
    hasfsample = ~isempty(options.fsample) ;
    
    % - No inputs -
    if ~hastime && ~hasNsample && ~hasfsample % doesn't have sufficient inputs
        error('Either the microstate object should have a time axis, or inputs Nsample and fsample should be specified')
    % - Good inputs -
    elseif hastime && ~hasNsample && ~hasfsample % get Nsample and fsample from the time axis
        options.Nsample = length(obj.time) ; 
        options.fsample = 1/mean(diff(obj.time)) ; 
    elseif ~hastime && hasNsample && hasfsample % make time axis from Nsample and fsample
        obj.time = (1/options.fsample)*(0:(Nsample-1)) ; 
    % - Deal with inputs we can work with -
    elseif hastime && hasNsample && ~hasfsample % Both Nsample and time specified, no fsample
        if options.Nsample ~= length(obj.time)
            warning('Both a time axis and Nsample were specified, but do not match. Using length(microstate.time) instead of Nsample for consistency.')
            options.Nsample = length(obj.time) ; 
        end
        options.fsample = 1/mean(diff(obj.time)) ; 
    elseif hastime && ~hasNsample && hasfsample % has fsample and time axis, but no Nsample
        options.Nsample = length(obj.time) ; 
        if options.fsample ~= 1/mean(diff(obj.time))
            warning('Both a time axis and fsample were specified, but do not match. Extracting fsample from time axis for consistency.')
            options.fsample = 1/mean(diff(obj.time)) ; 
        end
     elseif hastime && hasNsample && hasfsample
        if options.Nsample ~= length(obj.time)
            warning('Both a time axis and Nsample were specified, but do not match. Using length(microstate.time) instead of Nsample for consistency.')
            options.Nsample = length(obj.time) ; 
        end
        if options.fsample ~= 1/mean(diff(obj.time))
            warning('Both a time axis and fsample were specified, but do not match. Extracting fsample from time axis for consistency.')
            options.fsample = 1/mean(diff(obj.time)) ; 
        end
    % - Deal with inputs without sufficient information
    elseif ~hastime && hasNsample && ~hasfsample % Nsample specified, fsample not specified
        error('If Nsample is specified, fsample must also be specified') ; 
    elseif ~hastime && ~hasNsample && hasfsample % fsample specified, Nsample not specified
        error('If fsample is specified, Nsample must also be specified') ; 
    end
    
    % check whether both filterord and meandur were specified
    if isempty(options.meandur) && isempty(options.filterord) 
        options.meandur = 50e-3 ; % set default
    elseif ~isempty(options.meandur) && ~isempty(options.filterord)
        warning('Ignoring input meandur, using specified filterord')
        options.meandur = [] ; 
    end
    
    
    % build a tree for assigning randomwalks to microstate sequences
    numgroups = options.Nstates ; numlevels = 1 ; 
    for i = 1:options.Nstates; groups{1}{i} = i ; end
    while numgroups>1
        numwalksi = floor(numgroups/2) ;
        for i = 1:numwalksi
            groups{numlevels+1}{i} = [groups{numlevels}{2*i-1},groups{numlevels}{2*i}] ; 
        end
        for i = (2*floor(numgroups/2)+1):numgroups
            groups{numlevels+1}{end+1} = groups{numlevels}{i} ;
        end
        numlevels = numlevels+1 ; 
        numgroups = length(groups{numlevels}) ; 
    end
      
    for i = 1:options.Nstates-1 % loop over all random walks (there are always Nstates-1 randomwalks)
        % Use the wavelet toolbox to generate fractal brownian motion with Hurnst exponent 0.71 [1] 
        vi(i,:) = wfbm(options.Hurnst,options.Nsample+1) ; 
    end

    % find filter order
    if isempty(options.filterord)
        for filtord0 = 1:200
            [err(filtord0),md(filtord0),v(filtord0,:)] = makev(filtord0,vi) ;  
        end
        err = smooth(err,20) ; 
        [~,ind] = min(abs(err)) ; 
        v = v(ind,:) ; 
        options.filterord = ind ; 
    else
        v = makev(filterord,vi) ; 
    end
    
    obj.label = v ; 
    
    obj = microstate.functions.process_append(obj,'Added simulated (randomwalk) microstate sequence',options) ; 
        
    % Function to make + filter sequence
    function [err,md,v] = makev(filtord,vi0)

        vi1 = movmean(vi0',filtord)'; 
        x = diff(vi1,1,2)>0 ; 


        % Extract the microstate sequence from the random walk data
        v = zeros(1,options.Nsample) ; % will be the microstate sequence
        V = true(options.Nstates,options.Nsample) ; 
        groups0 = groups(1:end-1) ; 
        numwalksi = 0 ; 
        for i= length(groups0):-1:1
            Ni = length(groups0{i})/2 ; 
            for j = 1:Ni
                numwalksi = numwalksi+1 ; 
                ind1 = groups0{i}{2*j-1} ; 
                V(ind1,:) = V(ind1,:) & repmat(x(numwalksi,:),length(ind1),1) ;

                % remove any previous reference to this group as it is not a split
                for oi = 1:i-1
                    idx = false(1,length(groups0{oi})) ; 
                    for oj = 1:length(groups0{oi})
                        idx(oj) = isequal(groups0{oi}{oj},ind1) ; 
                    end
                    groups0{oi}(idx) = [] ; 
                end

                ind0 = groups0{i}{2*j} ; 
                V(ind0,:) = V(ind0,:) & repmat(~x(numwalksi,:),length(ind0),1) ; 

                % remove any previous reference to this group as it is not a split
                for oi = 1:i-1
                    idx = false(1,length(groups0{oi})) ; 
                    for oj = 1:length(groups0{oi})
                        idx(oj) = isequal(groups0{oi}{oj},ind0) ; 
                    end
                    groups0{oi}(idx) = [] ; 
                end
            end
        end

        [row,col] = find(V) ; 
        v(col) = row ; 

        % get mean duration
        md = mean(diff(find(diff(v)))/options.fsample) ;
        err = md-options.meandur ; 

        vi1 = [] ; 

    end


end

    