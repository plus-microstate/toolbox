function obj = preprocess_rereference(obj) ; 
% Re-reference data to average
    if isempty(obj.modality)
        return
    end
    switch obj.modality
        case 'eeg'
            obj.data = obj.data-mean(obj.data,2) ; 
            msg = 'Re-referenced to average' ; 
        case {'meg','source','ampenv'}
            msg = sprintf('Modality %s has no reference',obj.modality) ;
    end

    options.msg = msg ; 
    obj = microstate.functions.process_append(obj,'Applied re-reference',options) ;

end