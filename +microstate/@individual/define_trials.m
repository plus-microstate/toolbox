function trls = define_trials(obj,varargin) ; 

% Useage: 
% ms = ms.define_trials(epoch_length) ; % splits into trials of length
% epoch_length with no overlap
% 
% ms = ms.define_trials(epoch_length,overlap) ; % splits into trials of
% length epoch_length with overlap of overlap seconds.
% 
% ms = ms.define_trials(latency,stim_times) ; % splits into trials from A
% to B seconds following each stimulus, where latency = [A,B]. E.g. for
% -100 to 350 ms around the stimulus, specify latency = [-0.1,0.35]. 
% 
% ms = ms.define_trials(time_vec,stim_times) ; % splits into trials along
% the specified vector of times, where data/microstate labels are
% interpolated along these times. Times are relative to each stimulus, e.g.
% time_vec = -0.1:0.01:0.35 will give -100 ms to 350 ms around the stimulus
% in steps of 10ms. 
% 
% ms = ms.define_trials(_,stim_times,conditions) ; % as in the previous
% two, but labels for the conditions are specified

% Find out the inputs
switch nargin
    case 1
        error('Insufficient inputs to define trials, see useage description')
        
    case 2 % epoch_length
        method = 'fixed_epoch' ; 
        epoch_length = varargin{1} ; 
        overlap = 0 ; 
        
    case 3 % either epoch_length,overlap or latency,stim_times or t_vec,stim_times
        if isscalar(varargin{1}) && isscalar(varargin{2}) % epoch_length,overlap
            method = 'fixed_epoch' ; 
            epoch_length = varargin{1} ; 
            overlap = varargin{2} ; 
            
        elseif length(varargin{1})==2 % latency,stim_times
            method = 'latency' ; 
            latency = varargin{1}(:) ; 
            stim_times = varargin{2}(:) ; 
            conditions = ones(size(stim_times)) ; 
            
        else % t_vec,stim_times
            method = 'time_vec' ; 
            time_vec = varargin{1}(:) ; 
            stim_times = varargin{2}(:) ; 
            conditions = ones(size(stim_times)) ; 
        end
        
    case 4 % either latency,stim_times or t_vec,stim_times
        if length(varargin{1})==2 % latency,stim_times
            method = 'latency' ; 
            latency = varargin{1}(:) ; 
            stim_times = varargin{2}(:) ; 
            conditions = varargin{3}(:) ; 
            if length(conditions) ~= length(stim_times)
                error('length(conditions) should equal length(stim_times)')
            end
            
        else % t_vec,stim_times
            method = 'time_vec' ; 
            time_vec = varargin{1}(:) ; 
            stim_times = varargin{2}(:) ; 
            conditions = varargin{3}(:) ; 
            if length(conditions) ~= length(stim_times)
                error('length(conditions) should equal length(stim_times)')
            end
        end
end


% Make the trials
trls = microstate.cohort ; 
switch method
    case 'fixed_epoch'
        t_start = obj.time(1) ; 
        t_end = obj.time(end) ; 
        
        epoch(:,1) = t_start:(epoch_length-overlap):(t_end-epoch_length) ; 
        epoch(:,2) = epoch(:,1)+epoch_length ; 
        
        for i = 1:size(epoch,1)
            ms = microstate.functions.select_time(obj,epoch(i,:)) ;
            ms.time = ms.time-ms.time(1) ; 
            trls = trls.add_individuals(ms,[],1) ; 
        end
        
    case 'latency' 
            
        for i = 1:length(stim_times) ; 
            if (obj.time(1)>stim_times(i)+latency(1)) || (obj.time(end)<stim_times(i)+latency(2))
                continue
            end
            ms = microstate.functions.select_time(obj,stim_times(i)+latency) ; 
            ms.time = ms.time-stim_times(i) ; 
            trls = trls.add_individuals(ms,conditions(i),1) ; 
        end
        
    case 'time_vec'
        
        for i = 1:length(stim_times)
            if (obj.time(1)>stim_times(i)+time_vec(1)) || (obj.time(end)<stim_times(i)+time_vec(end))
                continue
            end
            
            x = [] ; lbl = [] ;
            if ~isempty(obj.data)
                x = interp1(obj.time,obj.data,stim_times(i)+time_vec); 
            end
            if ~isempty(obj.label)
                lbl = interp1(obj.time,obj.label,stim_times(i)+time_vec,'next') ; 
            end
            
            ms = microstate.individual(x,obj.modality,time_vec) ; 
            ms.label = lbl ; 
            trls = trls.add_individuals(ms,conditions(i),1) ; 
        end
        
end
        
        

end