function obj = simulate_seq_randomwalk(obj,varargin)
% Simulate a random-walk decision-tree microstate sequence

    % check inputs
    options = microstate.functions.make_options(varargin) ; 
    
    % default options
    defaults = {'Nsample',[] ; 
                'fsample',[] ; 
                'Nstates',4 ; 
                'Hurst',0.71 ; 
                'meandur',[] ; % will be set to 50 later if filterord wasn't specified
                'filterord',[] ; 
                'selectord','old' ; 
                'mindur',0 ; 
                'smoothpadding',1000} ; 
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
        obj.time = (1/options.fsample)*(0:(options.Nsample-1)) ; 
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
        % Use the wavelet toolbox to generate fractal brownian motion with Hurst exponent 0.71 [1] 
        vi(i,:) = wfbm(options.Hurst,options.Nsample+2*options.smoothpadding+1) ; 
    end

    % find filter order
    if isempty(options.filterord)
        if strcmpi(options.selectord,'ga')
            optimopts = optimoptions('ga','MaxGenerations',5,'PlotFcn','gaplotscores','PopulationSize',20) ; 
            filtord0 = ga(@(ord) makev(ord,vi),1,[],[],[],[],0,round(2*options.meandur*options.fsample),[],1,optimopts) ; 
            [err,md,v] = makev(filtord0,vi) ;
            options.filterord = filtord0 ; 
        elseif strcmpi(options.selectord,'steps')
            disp('Staircase to optimize mean microstate duration')
            disp('Iteration     |     Filter order     |     Duration')
            filtord0 = 2*round(options.meandur*options.fsample) ; 
            orderstried = [] ; 
            maxiter = 500 ; 
            iter = 0 ; 
            while (iter<maxiter) & ~any(orderstried==filtord0) ; 
                iter = iter+1 ; 
                orderstried(iter) = filtord0 ; 
                msg = sprintf('%d',iter) ; 
                msg = [msg,repmat(' ',1,14-length(msg)),'|     '] ; 
                msg = sprintf('%s%d',msg,filtord0) ; 
                msg = [msg,repmat(' ',1,37-length(msg)),'|     '] ; 
                [err,md,v] = makev(filtord0,vi) ;
                msg = sprintf('%s%.05f\n',msg,md) ; 
                fprintf(msg); 
                err = (md-options.meandur)*options.fsample ; 
                step = -ceil(1e-1*err) ; 
                filtord0 = max(filtord0+step , 1) ; 
                
            end
            options.filterord = filtord0 ; 
                
        else
            for filtord0 = 1:200
                [err(filtord0),md(filtord0),v(filtord0,:)] = makev(filtord0,vi) ;  
            end
            err = smooth(err,20) ; 
            [~,ind] = min(abs(err)) ; 
            v = v(ind,:) ; 
            options.filterord = ind ;
        end
    else
        v = makev(filterord,vi) ; 
    end
    
    obj.label = v ; 
    
    obj = microstate.functions.process_append(obj,'Added simulated (randomwalk) microstate sequence',options) ; 
        
    % Function to make + filter sequence
    function [err,md,v] = makev(filtord,vi0)
        % disp(num2str(filtord))
        try
            if filtord>0
                vi1 = movmean(vi0',filtord)';
            end
        vi1 = vi1(:,(options.smoothpadding+1):(end-options.smoothpadding)) ; 
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
        
        % get durations
        change = find(diff(v)) ; 
        change = [1,change+1 ; change , length(v)]; 
        dur = diff(change)/options.fsample ; 
        % ind = find(dur<options.mindur) ;
        [dur,ind2] = sort(dur) ; 
        % ind = ind(ind2) ; 
        while dur(1)<options.mindur
            
            try
                selind = randi(2) ; 
                sels = [v(change(2,ind2(1))+1) ; v(change(1,ind2(1))-1)];
                v(change(1,ind2(1)):change(2,ind2(1))) = sels(selind) ; 
            catch   
            try
                v(change(1,ind2(1)):change(2,ind2(1))) = v(change(2,ind2(1))+1) ;
            catch
                v(change(1,ind2(1)):change(2,ind2(1))) = v(change(1,ind2(1))-1) ;
            end
            end
            
            
            change = find(diff(v)) ; 
            change = [1,change+1 ; change , length(v)]; 
            dur = diff(change)/options.fsample ; 
            % ind = find(dur<options.mindur) ;
            [dur,ind2] = sort(dur) ; 
            % ind = ind(ind2) ; 
        end
            

        % get mean duration
        change = find(diff(v)) ;
        change = [1,change+1 ; change , length(v)]; 
        dur = diff(change)/options.fsample ;
        md = mean(dur) ;
        err = abs(md-options.meandur) ; 

        vi1 = [] ; 

        catch
            err = inf ; md = [] ; v = [] ; 
        end
    end


end

    