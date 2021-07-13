function s = exportstruct(obj)
    f = fieldnames(obj) ;
    s = struct ; 
    for i = 1:length(f)
        s.(f{i}) = obj.(f{i}) ; 
    end
end