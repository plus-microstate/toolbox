function options = make_options(opts_in)

    % Initialize options structure
    options = struct ; 

    % Deal with case of no options in
    if isempty(opts_in) 
        return
    end

    % Deal with structure options
    if length(opts_in) == 1 && isstruct(opts_in{1})
        opts_in = opts_in{1} ; 
    end

    % Options as a cell
    if iscell(opts_in)
        if length(opts_in) == 1 && isempty(opts_in{1})
            return
        end
        
        % Check options are in name value pairs
        if mod(length(opts_in),2) % check its even numbered
            error('Options must be name-value pairs')
        end
        opts_in = reshape(opts_in',2,length(opts_in)/2)' ; % make a column vector
        for i = 1:size(opts_in,1)
            if ~ischar(opts_in{i,1})
                error('Options must be name-value pairs') 
            end
        end

        % input options to options structure
        for i = 1:size(opts_in,1)
            options.(opts_in{i,1}) = opts_in{i,2} ; 
        end
        return
    end

    % Options as a structure
    if isstruct(opts_in)
        fields = fieldnames(opts_in) ; 
        for i = 1:length(fields)
            options.(fields{i}) = opts_in.(fields{i}) ; 
        end
    end

end % make_options