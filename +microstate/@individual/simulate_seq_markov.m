function obj = simulate_seq_markov(obj,varargin)
% Simulate a Markov chain microstate sequence
    % check inputs
    options = microstate.functions.make_options(varargin) ; 
    
    % default options
    defaults = {'Nsample',[] ; 
                'PriorProb',[]} ; 
    options = microstate.functions.add_options(options,defaults) ; clear defaults
    
    
    % check object and inputs match
    hastime = ~isempty(obj.time) || ~isempty(obj.label) ; 
    hasNsample = ~isempty(options.Nsample) ; 
    
    if ~hastime && ~hasNsample % doesn't have sufficient inputs
        error('Either the microstate object should have a time axis or input Nsample should be specified')
    elseif hastime && ~hasNsample % get Nsample and fsample from the time axis
        if ~isempty(obj.time)
            options.Nsample = length(obj.time) ; 
        else
            options.Nsample = length(obj.label) ; 
        end
    elseif hasNsample && ~hastime
        % do nothing, this is ok
    elseif hastime && hasNsample % Both Nsample and time specified
        if ~isempty(obj.time)
            if options.Nsample ~= length(obj.time)
                warning('Both a time axis and Nsample were specified, but do not match. Using length(microstate.time) instead of Nsample for consistency.')
                options.Nsample = length(obj.time) ; 
            end
        else
            if options.Nsample ~= length(obj.label)
                warning('Both a time axis and Nsample were specified, but do not match. Using length(microstate.time) instead of Nsample for consistency.')
                options.Nsample = length(obj.label) ; 
            end
        end
    end
    
    % check the microstate object has a Markov matrix or a sequence
    hasT = false ; 
    if isfield(obj.stats,'markov')
        if ~isempty(obj.stats.markov.matrix)
            hasT = true ; 
        end
    end
    haslabel = ~isempty(obj.label) ; 
    if ~hasT && ~haslabel
        error('The microstate object should either have a markov transitioning matrix or a microstate sequence (label) included')
    elseif hasT && haslabel
        warning('Overwriting current labels with new simulated labels')
    elseif ~hasT && haslabel
        obj = obj.stats_markov() ; % get Markov matrix
    end
    T = obj.stats.markov.matrix ; 
    
    % get number of classes
    Ns = size(T,1) ; 

    % set P0 if not specified
    if isempty(options.PriorProb)
        options.PriorProb = repmat(1/Ns,Ns,1) ; 
    else
        options.PriorProb = options.PriorProb/sum(options.PriorProb) ; % ensure it sums to 1
    end

    % initialize
    v = zeros(1,options.Nsample) ; % this will be our Markov generated microstate sequence
    x0 = find(rand < cumsum(options.PriorProb)) ; v(1) = x0(1) ; % draw the first microstate from the P0 distribution 

    % generate Markov chain
    for i = 1:(options.Nsample-1)
        x0 = v(i) ; % current state
        Ti = T(x0,:) ; % probability of transition
        x1 = find(rand < cumsum(Ti)) ; x1 = x1(1) ; % draw next state
        v(i+1) = x1; % assign to time series
    end
    
    obj.label = v ; 
    obj = microstate.functions.process_append(obj,'Added simulated (Markov) microstate sequence',options) ; 


end

    