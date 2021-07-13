% Import struct
function obj = importstruct(obj,s)
    f = fieldnames(s) ;
    for i = 1:length(f)
        obj.(f{i}) = s.(f{i}) ; 
    end
end