function obj = stats_duration(obj)
% Calculate durations of microstates
    % Check inputs
    if isempty(obj.time) || isempty(obj.label)
        error('To calculation microstate duration, properties time and label are required')
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
    
    % Get durations and their labels
    transition = find(diff(obj.label)) ; 
    dv = diff(obj.time(transition)) ; % duration of each state
    label = obj.label(transition) ; label = label(2:end) ;  % label of each state
    % note, first and final states ignored in this code. this is good,
    % since we don't know when they started/ended respectively. 
    
    % Get mean duration
    md = mean(dv) ; 
    
    % get durations of each state
    dur = zeros(1,Nstates) ; 
    for i = 1:Nstates
        dur(i) = mean(dv(label == i)) ; 
    end
    
    % save to microstate object
    obj.stats.mean_duration = md ; 
    obj.stats.duration = dur ; 
    
    % Append process
    obj = microstate.functions.process_append(obj,'Calculated statistic: durations') ; 
    
end