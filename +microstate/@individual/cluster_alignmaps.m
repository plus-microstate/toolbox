function [obj,mapsim] = cluster_alignmaps(obj,usepeaks,varargin) ; 
% Assign a microstate class to each time point, building a discrete microstate time series.
% 
% Given a set of template maps c, a set of maps at global field potential
% peaks X, and some way to map the peak maps onto the data v, we can
% calculate a time series which classifies each data point into a
% microstate. This is done by finding the template map each peak map is
% closest to by taking the absolute value of correlation [1]. 

% Luke Tait 15/05/2018
%%  Check inputs

    if nargin < 2 || isempty(usepeaks)
        usepeaks = true ; 
    end
    
    options = microstate.functions.make_options(varargin) ; 

    % default options
    defaults = {'minmslength',0; 
                'smooth',0 ; 
                'corrthresh',0 ; 
                'keep_polarity',false ;
                } ;  
    options = microstate.functions.add_options(options,defaults) ; clear defaults
    

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
            if options.keep_polarity
                distfun = @(x,c) 1-( (x./vecnorm(x,2,2)) * ((c-mean(c,2))./vecnorm((c-mean(c,2)),2,2))' )' ;
            else
                distfun = @(x,c) 1-abs( (x./vecnorm(x,2,2)) * ((c-mean(c,2))./vecnorm((c-mean(c,2)),2,2))' )' ; 
            end
        case 'meg'
            % do nothing
            if options.keep_polarity
                distfun = @(x,c) 1-( (x./vecnorm(x,2,2)) * (c./vecnorm(c,2,2))' )' ; 
            else
                distfun = @(x,c) 1-abs( (x./vecnorm(x,2,2)) * (c./vecnorm(c,2,2))' )' ; 
            end
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
    
    % interpolate over all samples, not just peaks, using peak2sample
    obj.label = zeros(size(peak2sample)) ; % initialize output V
    msprob = zeros(k,length(peak2sample)) ; 
    for i = 1:n % loop over each gfp peak
        ind = peak2sample == i ; % find time points corresponding to peak
        msprob(:,ind) = repmat(mapsim(:,i),1,sum(ind)) ; 
        % obj.label(ind) = idx(i) ; % assign each of these time points the class corresponding to the peak
    end
    
    % options.smooth
    if options.smooth>0
        msprob = smoothdata(msprob','gaussian',options.smooth)' ; 
    end
    
    % get sequence
    [msprob,lbl] = max(msprob,[],1) ; 
    
    % remove short microstates
    if options.minmslength>0
        transition = [0,find(diff(lbl)~=0),length(lbl)] ; 
        difflbl = diff(transition) ; 
        shortdur = find(difflbl<options.minmslength) ;
        inds = [transition(1:end-1)+1 ; transition(2:end)] ; 
        for i = 1:length(shortdur)
            lbl(inds(1,shortdur(i)):inds(2,shortdur(i))) = nan ; 
        end
        
        % interpolate the short states
        lbl = interp1(find(~isnan(lbl)),lbl(~isnan(lbl)),1:length(lbl),'nearest') ; 
    
    end
    
    % threshold
    if options.corrthresh > 0 
        lbl(msprob<options.corrthresh) = nan ; 
    end
    
    
    obj.label = lbl ; 
    
    
    
    options = struct ; 
    options.usepeaks = usepeaks ; 
    obj = microstate.functions.process_append(obj,'Backfit maps to data',options) ;
    
end

