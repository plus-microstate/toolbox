function obj = stats_gev(obj,addprocess) ; 
% Calculate GEV of microstates
    % check we have GFP, maps, and labels
    if isempty(obj.gfp)
        obj = obj.calculate_gfp ; 
    end
    if isempty(obj.maps)
        error('Maps are required to calculate GEV')
    end
    if isempty(obj.label)
        obj = obj.cluster_alignmaps ; 
    end
    
    % Get data
    x = obj.data ; 
    c = obj.maps ; 
    
    % get 'correlation' functions
    switch obj.modality
        case 'eeg'
            x = x-mean(x,2) ; % re-reference to average
            corrfun = @(x,c) abs( (x./vecnorm(x,2,2)) * ((c-mean(c,2))./vecnorm((c-mean(c,2)),2,2))' )' ; 
        case 'meg'
            % do nothing
            corrfun = @(x,c) abs( (x./vecnorm(x,2,2)) * (c./vecnorm(c,2,2))' )' ; 
        case 'source'
            x = abs(x) ; % .^2 ; % set to magnitude
            c = abs(c) ; % .^2 ; 
            corrfun = @(x,c) ( (x./vecnorm(x,2,2)) * (c./vecnorm(c,2,2))' )' ;
        case 'ampenv'
            % do nothing
            corrfun = @(x,c) ( (x./vecnorm(x,2,2)) * (c./vecnorm(c,2,2))' )' ; 
    end
    
    % calculate GEV
    denom = nansum(obj.gfp.^2) ; 
    gev = 0 ; 
    for i = 1:size(c,2)
        ind = find(obj.label == i) ; 
        if isempty(ind)
            continue
        end
        xi = x(ind,:) ; 
        gev = gev + nansum( (obj.gfp(ind).^2).*(corrfun(xi,c(:,i)').^2) ) ; 
    end
    gev = gev/denom ; 
    
    
    % Output options
    obj.stats.gev = gev ; 
    if nargin < 2
        addprocess = false ; 
    end
    if addprocess
        obj = microstate.functions.process_append(obj,'Calculated GEV') ;
    end

end