function obj = preprocess_demean(obj,polynomialord) ; 
% Re-reference data to average
    if isempty(obj.modality)
        return
    end
    if nargin<2 || isempty(polynomialord)
        polynomialord=1 ; 
    end
    obj.data = detrend(obj.data,polynomialord) ; 

    obj = microstate.functions.process_append(obj,'Detrended',struct('polynomialord',polynomialord)) ;

end