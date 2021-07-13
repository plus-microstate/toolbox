function [obj,mapsim] = cluster_alignmaps(obj,usepeaks)
% Assign a microstate class to each time point, building a discrete microstate time series.
% 
% Given a set of template maps c, a set of maps at global field potential
% peaks X, and some way to map the peak maps onto the data v, we can
% calculate a time series which classifies each data point into a
% microstate. This is done by finding the template map each peak map is
% closest to by taking the absolute value of correlation [1]. 

% Luke Tait 15/05/2018
%%  Check inputs

    if nargin < 2
        usepeaks = true ; 
    end

    % check there are maps
    if isempty(obj.maps)
        error('No microstate maps included in individual object. Use function cluster_estimatemaps to estimate the maps')
    end
    
    % check there is data
    if isempty(obj.data)
        error('Data is required')
    end
    
    % Get dimensionality
    [n,p] = size(obj.data) ; % n: number of points, p: dimensionality of data
    [pc,k] = size(obj.maps) ; % k: number of maps, pc: dimensionality of data 
    if pc ~= p
        error('Maps should be of the same dimensionality as the data')
    end

    if usepeaks
        % calculate gfp
        if isempty(obj.gfp)
            obj = obj.calculate_gfp ; 
        end

        % Assign classes
        [xpeak,~,peak2sample] = obj.cluster_get_gfppeaks ; % calculate GFP peaks
        n = size(xpeak,1) ; % new value of n
    else
        xpeak = obj.data ; 
        peak2sample = 1:n ; 
    end
    
    % get maps
    maps = obj.maps ; 
    

    
    % Find class of each map by calculating to which template the map has
    % minimum distance
    
    % get the distance function
    switch obj.modality
        case 'eeg'
            xpeak = xpeak-mean(xpeak,2) ; % re-reference to average
            distfun = @(x,c) 1-abs( (x./vecnorm(x,2,2)) * ((c-mean(c,2))./vecnorm((c-mean(c,2)),2,2))' )' ; 
        case 'meg'
            % do nothing
            distfun = @(x,c) 1-abs( (x./vecnorm(x,2,2)) * (c./vecnorm(c,2,2))' )' ; 
        case 'source'
            xpeak = abs(xpeak) ; % set to magnitude
            maps = abs(maps) ; 
            distfun = @(x,c) 1-( (x./vecnorm(x,2,2)) * (c./vecnorm(c,2,2))' )' ;
        case 'ampenv'
            % do nothing
            maps = abs(maps) ; % set to amplitude
            distfun = @(x,c) 1-( (x./vecnorm(x,2,2)) * (c./vecnorm(c,2,2))' )' ; 
    end
    d = distfun(xpeak,maps') ; 
    mapsim = 1-d ; 
    [~,idx] = min(d,[],1) ; 

    % interpolate over all samples, not just peaks, using peak2sample
    obj.label = zeros(size(peak2sample)) ; % initialize output V
    for i = 1:n % loop over each gfp peak
        ind = peak2sample == i ; % find time points corresponding to peak
        obj.label(ind) = idx(i) ; % assign each of these time points the class corresponding to the peak
    end
    
    
    
    options = struct ; 
    options.usepeaks = usepeaks ; 
    obj = microstate.functions.process_append(obj,'Backfit maps to data',options) ;
    
end

