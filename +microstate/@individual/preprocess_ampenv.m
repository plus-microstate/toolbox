function obj = preprocess_ampenv(obj,varargin)
% Calculate amplitude envelope of data

    % check inputs
    options = microstate.functions.make_options(varargin) ; 
    
    % default options
    defaults = {'bpfilter',[] ; 
                'resample',[]} ; 
    options = microstate.functions.add_options(options,defaults) ; clear defaults

    % get size of data
    [m,n] = size(obj.data) ;
    
    % bandpass filter
    if ~isempty(options.bpfilter)
        validateattributes(options.bpfilter,{'double'},{'size',[1,2]},'preprocess_convert2ampenv','bpfilt')
        obj = obj.preprocess_filter(options.bpfilter(1),options.bpfilter(2)) ;
    end
    
    % Check fieldtrip external hilbert is not on path
    fpath = fileparts(which('hilbert')) ; 
    onpath = contains(fpath,'fieldtrip') ; 
    if onpath
        rmpath(fpath) ; 
    end
        
    % get amplitude envelope
    obj.data = abs(hilbert(obj.data)) ; 
    obj.modality = 'ampenv' ; 
    
    % add path back if needs be
    if onpath
        addpath(fpath) ; 
    end
    
    % update the process
    obj = microstate.functions.process_append(obj,'Calculated amplitude envelope') ; % options are included in the subfunctions
    
    % resample
    if ~isempty(options.resample)
        validateattributes(options.resample,{'double'},{'scalar','finite','real','positive'},'preprocess_convert2ampenv','resample')
        obj = obj.preprocess_resample(options.resample) ; 
    end
    
    % update the gfp
    obj = obj.calculate_gfp() ; 
    
end