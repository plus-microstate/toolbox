function obj = select_time(obj,time)

switch class(obj)
    
    case 'microstate.cohort'
        for i = 1:length(obj.individual) 
            obj.individual(i) = microstate.functions.select_time(obj.individual(i),time) ; 
        end
        
    case 'microstate.individual' 
        
        idx = (obj.time>=time(1)) ... 
                & (obj.time<=time(2)) ;
        if ~isempty(obj.data)
            obj.data = obj.data(idx,:) ; 
        end
        if ~isempty(obj.time)
            obj.time = obj.time(idx) ;
        end
        if ~isempty(obj.label)
            obj.label = obj.label(idx) ; 
        end
        if ~isempty(obj.sample)
            obj.sample = obj.sample(idx) ;
        end
        if ~isempty(obj.gfp)
            obj.gfp = obj.gfp(idx) ;
        end
end 
    