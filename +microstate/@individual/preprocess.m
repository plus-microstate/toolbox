function obj = preprocess(obj)
% Wrapper function for a default preprocessing pipeline. 
% The pipeline will filter 1-30 Hz, re-reference to average (EEG only), and
% resample at 256 Hz (unless sampling rate is lower than this value).

    obj = obj.preprocess_filter(1,30) ; % filter 1-30 Hz
    obj = obj.preprocess_rereference ; % if EEG, re-reference to average
    obj = obj.preprocess_resample(min(256,1/mean(diff(obj.time)))) ; % resample at 256 if current sampling rate is higher
    
end