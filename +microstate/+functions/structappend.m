function s = structappend(s0,s1) ;

    f1 = fieldnames(s1) ; 
    s = s0 ; 

    for i = 1:length(f1)

        % deal with matrices
        if strcmp(class(s1.(f1{i})),'double') ; 
            dim = size(s1.(f1{i})) ; 
            if (length(dim) == 2) && (dim(2) == 1) % column vector
                s1.(f1{i}) = s1.(f1{i})' ; 
            elseif (length(dim) == 2) && (dim(1) == 1)
                % do nothing
            elseif length(dim)>=2 % matrices
                s1.(f1{i}) = permute(s1.(f1{i}),[length(dim)+1,1:length(dim)]) ; 
            end
        end

        % deal with the case where there is nothing to be updated
        if ~isfield(s0,f1{i})
            if ~isstruct(s1.(f1{i}))
                s.(f1{i}) = s1.(f1{i}) ; 
                continue
            else
                s.(f1{i}) = struct ; 
                s.(f1{i}) = microstate.functions.structappend(s.(f1{i}),s1.(f1{i})) ; 
                continue
            end
        end

        % now, check that both structures think this field is the same type
        % (e.g. you cannot append a cell onto a matrix)
        if ~strcmp(class(s0.(f1{i})),class(s1.(f1{i})))
            error('Cannot append field %s since it is type %s in the original structure and type %s in the structure to be appended',f1{i},class(s0.(f1{i})),class(s1.(f1{i}))) ; 
        end

        % append to the end
        switch class(s1.(f1{i}))
            case 'struct'
                % this is the case of substructures - let us recursively call
                % the function to append substructures
                s.(f1{i}) = microstate.functions.structappend(s0.(f1{i}),s1.(f1{i})) ; 

            case {'double','cell'}
                % ensure cells being concatenated are the same length
                if iscell(s1.(f1{i}))
                    l1 = size(s1.(f1{i}),2) ; 
                    l0 = size(s0.(f1{i}),2) ; 
                    if l0>l1
                        for j = 1:size(s1.(f1{i}),1)
                            for k = (l1+1):l0
                                s1.(f1{i}){j,k} = [] ; 
                            end
                        end
                    elseif l1>l0
                        for j = 1:size(s0.(f1{i}),1)
                            for k = (l0+1):l1
                                s0.(f1{i}){j,k} = [] ; 
                            end
                        end
                    end
                end

                s.(f1{i}) = [s0.(f1{i}) ; s1.(f1{i})] ; 

            otherwise
                try
                    s.(f1{i}) = [s0.(f1{i}) ; s1.(f1{i})] ; 
                catch
                    warning('Field %s not appended as classes %s are currently not supported in the structappend function',f1{i},class(s0.(f1{i}))) ; 
                end
        end

    end 
end % structappend