function obj = stats_coverage(obj)
% Calculate coverage of microstates
    % Check inputs
    if isempty(obj.label)
        error('To calculation microstate coverage, property label is required')
    end
    
    % Number of states
    if ~isempty(obj.maps)
        Nstates = size(obj.maps,2) ; 
    else
        Nstates = max(obj.label) ; 
    end
    
    % Get coverage of each microstate
    coverage = zeros(1,Nstates) ; 
    for i = 1:Nstates
        coverage(i) = sum(obj.label == i)/length(obj.label) ; 
    end
    
    % save to microstate object
    obj.stats.coverage = coverage ; 
    
    % Append to process
    obj = microstate.functions.process_append(obj,'Calculated statistic: coverage') ; 
    
end