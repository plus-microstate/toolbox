function obj = preprocess_demean(obj) ; 
% Re-reference data to average
    if isempty(obj.modality)
        return
    end
    obj.data = obj.data-mean(obj.data,1) ; 

    options.msg = msg ; 
    obj = microstate.functions.process_append(obj,'Applied re-reference',options) ;

end