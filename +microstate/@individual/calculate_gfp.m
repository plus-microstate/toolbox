function obj = calculate_gfp(obj) ; 
% Calculate GFP of time series
    if isempty(obj.modality)
        return
    end
    switch obj.modality
        case 'eeg'
            gfp = std(obj.data,[],2) ;
            method = 'std' ; 
        case {'meg','source','ampenv'}
            gfp = vecnorm(obj.data,2,2) ; 
            method = 'vecnorm' ; 
    end
    gfp(obj.bad_samples) = nan ;
    
    options = struct ; 
    options.method = method ; 

    obj.gfp = gfp ; 
    obj = microstate.functions.process_append(obj,'Calculated GFP',options) ;

end