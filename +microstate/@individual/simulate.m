function obj = simulate(obj,varargin) ; 
% Wrapper for default simulation pipelines. 
% By default, this method will: 
% - Generate time axis. If one is a property of the microstate object, it
%   will use this. Otherwise, 2 seconds of data at a sampling rate of 256
%   will be simulated. 
% - Get number of states and number of ROIs. If microstate.maps is
%   specified, then these will directly be taken from this property.
%   Otherwise, the following properties will be checked (in the order
%   listed) to try to choose the number of states: label, markov_matrix,
%   syntax_matrix. If each of these are empty, by default 4 states are
%   used. If maps aren't specified, then maps will be set to the identity
%   matrix and the same number of ROIs as states will be used. 
% - If no microstate sequence (microstate.label) is supplied, generate one.
%   If the Markov matrix is supplied, this will be used to generate the
%   labels, otherwise a randomwalk will be. 
% - Simulate the data using the default parameters.

    % Check for maximum time
    idx = find(strcmp(varargin,'tstart')) ; 
    if ~isempty(idx)
        tstart = varargin{idx+1} ; 
    else
        tstart = 0 ; 
    end
    
    % Check for start time
    idx = find(strcmp(varargin,'tend')) ; 
    if ~isempty(idx)
        tend = varargin{idx+1} ; 
    else
        tend = 60 ; 
    end
    
    % Check for sample rate
    idx = find(strcmp(varargin,'fsample')) ; 
    if ~isempty(idx)
        fsample = varargin{idx+1} ; 
    else
        fsample = 256 ; 
    end
    
    % Get time statistics
    if isempty(obj.time)
        obj.time = tstart:(1/fsample):tend ; % 2 seconds by default
    end
    fsample = 1/mean(diff(obj.time)) ; 
    Nt = length(obj.time) ; 
    
    % Get number of states
    if ~isempty(obj.maps)
        [~,Ns] = size(obj.maps) ; 
    elseif isempty(obj.maps) && ~isempty(obj.label)
        Ns = max(obj.label) ; 
        obj.maps = eye(Ns) ; 
    elseif isempty(obj.maps) && ~isempty(obj.markov_matrix)
        Ns = size(obj.markov_matrix,1) ; 
        obj.maps = eye(Ns) ; 
    elseif isempty(obj.maps) && ~isempty(obj.syntax_matrix)
        Ns = size(obj.syntax_matrix,1) ; 
        obj.maps = eye(Ns) ; 
    else
        Ns = 4 ; % by default, use 4 states
        obj.maps = eye(Ns) ; 
    end
        
    % Simulate microstate sequence - if Markov matrix supplied use this,
    % otherwise use randomwalk method
    if isempty(obj.label)
        if isfield(obj.stats,'markov')
            obj = obj.simulate_seq_markov(varargin{:}) ; % generate markov sequence
        else
            obj = obj.simulate_seq_randomwalk('Nstates',Ns,varargin{:}) ; % generate randomwalk
        end
    end
    
    % Simulate data
    obj = obj.simulate_data(varargin{:}) ; 
    
end
    