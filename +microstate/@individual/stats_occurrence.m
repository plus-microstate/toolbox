function obj = stats_occurrence(obj)
% Calculate number of occurrences of a microstate per second
    % Check inputs
    if isempty(obj.time) || isempty(obj.label)
        error('To calculation microstate occurrence, properties time and label are required')
    end
    if length(obj.time) ~= length(obj.label)
        error('microstate.label must be the same length as microstate.time')
    end
    
    % Number of states
    if ~isempty(obj.maps)
        Nstates = size(obj.maps,2) ; 
    else
        Nstates = max(obj.label) ; 
    end
    
    % Get microstate transitioning sequence
    transition = diff(obj.label) ~= 0 ; 
    label = [obj.label(transition),obj.label(end)] ; 
    
    % get occurrences of each state
    numstate = zeros(1,Nstates) ; 
    for i = 1:max(Nstates)
        numstate(i) = sum(label == i) ; 
    end
    occ = numstate/(obj.time(end)-obj.time(1)) ; 
    
    % save to microstate object
    obj.stats.occurrence = occ ; 
    
    % Append process
    obj = microstate.functions.process_append(obj,'Calculated statistic: occurrence') ; 
    
end