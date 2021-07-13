function obj = preprocess_resample(obj,fsample_new)
% Resample data
    % get size of data
    [m,n] = size(obj.data) ;
    
    % get sampling frequency
    fsample_old = 1/mean(diff(obj.time)) ; 
    if fsample_old == fsample_new
        return
    end

    % resample
    obj.data = resample(obj.data,round(fsample_new),round(fsample_old)) ; 
    tnew = (1/fsample_new)*(0:(size(obj.data,1)-1)) + obj.time(1) ;
    
    % recalculate microstate labels
    if ~isempty(obj.label)
        obj.label = interp1(obj.time',obj.label',tnew')' ; 
    end
    obj.time = tnew ; 
    
    % update the process
    options = struct ; 
    options.fsample_old = fsample_old ; options.fsample_new = fsample_new ; 
    obj = microstate.functions.process_append(obj,'Resampled data',options) ; 
    
    % update the gfp
    obj = obj.calculate_gfp() ; 
end