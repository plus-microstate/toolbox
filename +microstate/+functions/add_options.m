% Add options
function old_options = add_options(old_options,new_options,overwrite)

    % Check whether to overwrite
    if nargin < 3
        overwrite = false ; 
    end
    validateattributes(overwrite,{'logical'},{'scalar'})

    % new options can be in one of three formats
    if isstruct(new_options)
        % do nothing
    elseif iscell(new_options) && size(new_options,1) == 1
        new_options = microstate.functions.make_options(new_options) ; 
    elseif iscell(new_options) && size(new_options,2) == 2
        new_options = new_options' ; 
        new_options = microstate.functions.make_options(new_options(:)') ; 
    else
        error('incorrect format for new_options')
    end

    % Check both options inputs are structures
    validateattributes(old_options,{'struct'},{'scalar'})
    validateattributes(new_options,{'struct'},{'scalar'})

    % Get fields in new options
    fields = fieldnames(new_options) ; 

    % Update new options
    if overwrite % overwrite old options with new options
        for i = 1:length(fields)
            old_options.(fields{i}) = new_options.(fields{i}) ; 
        end
        return
    else % do keep old options (e.g. if new options are defaults)
        for i = 1:length(fields)
            if ~isfield(old_options,fields{i}) % only add the new option if it isn't already in old_options
                old_options.(fields{i}) = new_options.(fields{i}) ;
            end
        end
    end

end % add_options