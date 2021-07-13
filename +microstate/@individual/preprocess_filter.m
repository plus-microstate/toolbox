function obj = preprocess_filter(obj,flow,fhigh,varargin)
% Bandpass or bandstop filter data

    % check inputs
    options = microstate.functions.make_options(varargin) ; 
    
    % default options
    defaults = {'type','pass' ; % pass or stop
                'order',2 ; % 4th order
                } ; 
    options = microstate.functions.add_options(options,defaults) ; clear defaults

    % get size of data
    [m,~] = size(obj.data) ;
    
    % get sampling frequency
    fsample = 1/mean(diff(obj.time)) ; 

    % nyquist frequency
    fn = fsample/2;
    
    % compute filter coefficients
    if strcmp(options.type,'pass')
        if flow == 0 % lowpass filter
            [b,a] = butter(options.order,fhigh/fn) ; 
        elseif fhigh > fn % highpass filter
            [b,a] = butter(options.order,flow/fn,'high') ; 
        else % bandpass filter
            [b,a] = butter(options.order,[flow/fn,fhigh/fn]);
        end
    elseif strcmp(options.type,'stop')
        [b,a] = butter(options.order,[flow/fn,fhigh/fn],'stop') ; 
    end

    % remove mean (to be added after)
    mx = mean(obj.data,1) ; 
    obj.data = obj.data-repmat(mx,m,1) ; 

    % Check fieldtrip external hilbert is not on path
    fpath = fileparts(which('filtfilt')) ; 
    onpath = contains(fpath,'fieldtrip') ; 
    if onpath
        rmpath(fpath) ; 
    end
        
    % filter
    x_filt = filtfilt(b,a,obj.data) ; 
    
    % add means back 
    obj.data = x_filt+repmat(mx,m,1) ; 
    
    % add path back if needs be
    if onpath
        addpath(fpath) ; 
    end
    
    % update the process
    options.flow = flow ; options.fhigh = fhigh ; options.fsample = fsample ; 
    obj = microstate.functions.process_append(obj,'Filtered data',options) ; 
    
    % update the gfp
    obj = obj.calculate_gfp() ; 
end