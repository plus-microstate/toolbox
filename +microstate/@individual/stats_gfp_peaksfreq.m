function obj = stats_gfp_peaksfreq(obj)
% Calculate number of GFP peaks per second
    if isempty(obj.gfp)
        obj = obj.calculate_gfp ; 
    end

    % Check inputs
    if isempty(obj.time)
        error('To calculation gfp peaks per second, time is required')
    end
    if length(obj.time) ~= length(obj.gfp)
        error('microstate.gfp must be the same length as microstate.time')
    end
    
    % Find peaks of gfp
    [~,pks] = findpeaks(obj.gfp) ; % locations of peaks
    pks = length(pks) ; % number of peaks
    pks = pks/(obj.time(end)-obj.time(1)) ; % number of peaks per second
    
    % update output
    obj.stats.gfp_peaksfreq = pks ; 
    
    % Append process
    obj = microstate.functions.process_append(obj,'Calculated statistic: gfp peaks per second') ; 
 
end

    