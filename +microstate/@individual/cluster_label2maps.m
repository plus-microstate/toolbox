function [obj,mapsim] = cluster_label2maps(obj,usepeaks,keep_polarity)
% Calculate centroid maps given data and microstate labels

%%  Check inputs

    if nargin < 3
        keep_polarity = false ; 
    end
    if nargin < 2
        usepeaks = true ; 
    end

    % check there are maps
    if isempty(obj.label)
        error('No microstate labels provided')
    end
    
    % check there is data
    if isempty(obj.data)
        error('Data is required')
    end
    
    % Get dimensionality
    [n,p] = size(obj.data) ; % n: number of points, p: dimensionality of data
    k = max(obj.label) ; % k: number of maps, pc: dimensionality of data 

    if usepeaks
        % calculate gfp
        if isempty(obj.gfp)
            obj = obj.calculate_gfp ; 
        end

        % Assign classes
        [xpeak,peaks] = obj.cluster_get_gfppeaks ; % calculate GFP peaks
        n = size(xpeak,1) ; % new value of n
    else
        xpeak = obj.data ; 
        peaks = 1:n ; 
    end
    
    % get maps
    label = obj.label(peaks) ; 
    
    
    % get the distance function
    switch obj.modality
        case 'eeg'
            xpeak = xpeak-mean(xpeak,2) ; % re-reference to average
            centfun = @(v1) v1 ; 
        case 'meg'
            % do nothing
            centfun = @(v1) v1 ; 
        case 'source'
            xpeak = abs(xpeak) ; % set to magnitude
            centfun = @(v1) abs(v1) ; 
        case 'ampenv'
            % do nothing
            centfun = @(v1) abs(v1) ; 
    end
    
    % loop over maps
    for i = 1:k
        
        Xj = xpeak(label == i,:) ; % Maps previously assigned to cluster j
        if isempty(Xj) % Not found any in cluster
            v1 = nan(p,1) ; 
        else
            [V,D] = eig(Xj'*Xj) ; 
            [~,indmax] = max(diag(D)) ; % should just be the biggest value, but put this in in case
            v1 = V(:,indmax) ; % get first eigenvector
            
            if keep_polarity
                C = v1'*Xj' ; 
                v1 = sign(mean(C))*v1 ; 
            end
        end
        c(i,:) = centfun(v1) ; % assign as centroid
        
    end
    
    obj.maps = c' ; 
    
    options = struct ; 
    options.usepeaks = usepeaks ; 
    obj = microstate.functions.process_append(obj,'Calculated maps from labels',options) ;
    

end

