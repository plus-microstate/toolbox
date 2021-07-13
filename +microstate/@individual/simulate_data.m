function obj = simulate_data(obj,varargin)
% Simulate time series data
    % check inputs
    options = microstate.functions.make_options(varargin) ; 
    
    % default options
    defaults = {'params',struct ; 
                'fsample',[]} ; 
    options = microstate.functions.add_options(options,defaults) ; clear defaults
    
    % check obj has the right values
    if isempty(obj.label)
        error('The microstate object must have labels corresponding to which state is active')
    end
    
    % Get number of states to simulates
    Ns = max(obj.label) ; 
    
    % check for either time axis or fsample
    hastime = ~isempty(obj.time) ; 
    hasfsample = ~isempty(options.fsample) ; 
    
    if hastime && hasfsample
        if options.fsample ~= 1/mean(diff(obj.time))
            warning('Both a time axis and fsample were specified, but do not match. Extracting fsample from time axis for consistency.')
            options.fsample = 1/mean(diff(obj.time)) ; 
        end
    elseif hastime && ~hasfsample
        options.fsample = 1/mean(diff(obj.time)) ; 
    elseif ~hastime && hasfsample
        obj.time = (1/fsample)*(1:length(obj.label)) ; 
    elseif ~hastime && ~hasfsample
        error('Either a time axis or a sampling frequency is required')
    end
    
    % check for maps
    if isempty(obj.maps)
        obj.maps = eye(Ns) ; 
    end
        
    % make parameters
    options.tend = obj.time(end) ;
    options.tstart = obj.time(1)-2 ; % 2 seconds initialization
    options.W = obj.maps ; 
    options = wc_simulation_parameters(options) ; 
    obj.data = simulate_microstates(options,[obj.label ; obj.time]) ; 
    
    % update process
    obj = microstate.functions.process_append(obj,'Simulated data',options) ; 
    
    % re-calculate gfp
    obj = obj.calculate_gfp() ; 
    
    

    
    
    % Wilson-Cowan parameters
    function params = wc_simulation_parameters(varargin) ; 
 
        params_in = checkoptions(varargin) ; 

        % make default parameters 
        params.tauE = 10 ; % excitatory time constant, ms
        params.tauI = 20 ; % inverse excitatory synaptic time constant, s^-1
        params.a = 0.1 ; % rewrite to a[w+*E - I + P - b]
        params.b = 40 ; % inverse inhibitory synaptic time constant, s^-1
        params.c = 100 ; % connectivity between populations, a.u. 
        params.wI = 1.5 ;
        params.wplus = 1.4 ; 
        params.Poff = 10 ; 
        params.Pon = 13 ; % mean firing rate of thalamus, Hz
        params.sigma = sqrt(0.01) ; 
        params.SNR = db2pow(0) ; 
        params.W = 1 ; 

        params.tstart = -2 ; % time to start. negative values are to eliminate transients
        params.tend = 60 ; % time to end simulation
        params.dt = 0.1e-3 ; % maximum time step for any integration
        params.fsample = 256 ; 

        % add custom parameters to the params structure
        vars = fieldnames(params_in);  
        for i = 1:length(vars)
            params.(vars{i}) = params_in.(vars{i}) ; 
        end

        % add parameters which are derived from other parameters - typically
        % these are not changed, but the option to change them is available
        if ~isfield(params,'f')
            params.f = @(v) params.c./(1+exp(-params.a*(v-params.b))) ; 
        end
        params.Ns = size(params.W,2) ; 
        params.Nr = size(params.W,1) ; 
        params.P0 = repmat(1/params.Ns,params.Ns,1) ; 

        function options = checkoptions(optionsin)
            errormsg = 'inputs must be specified as a structure or in name-value pairs' ; 
            % make options structure
            if isempty(optionsin)
                options = struct ; 
            elseif isstruct(optionsin)
                options = optionsin ; 
            elseif iscell(optionsin)
                if length(optionsin)==1
                    if isstruct(optionsin{1})
                        options = optionsin{1} ; 
                    else
                        error(errormsg)
                    end
                elseif length(optionsin)/2==floor(length(optionsin)/2) ; % check even number of pairs
                    optionsin = optionsin(:) ; 
                    optionsin = reshape(optionsin,2,length(optionsin)/2) ; 
                    options = struct ; 
                    for i = 1:size(optionsin,2)
                        options.(optionsin{1,i}) = optionsin{2,i} ; 
                    end
                else
                    error(errormsg)
                end
            else
                error(errormsg)
            end
        end
    end


    % Simulation function --------------
    function sim = simulate_microstates(params,label)


        % Generate time axis
        tpoints = 1000*(params.tstart:params.dt:params.tend) ; % simulation time axis
        Nt = length(tpoints) ; % number of simulation time point
        
        % interpolate label onto this new time axis
        tlabel = label(2,:) ; 
        label = label(1,:) ; 
        simlabel = interp1(1000*tlabel,label,tpoints,'nearest') ; 
        simlabel(tpoints<1000*tlabel(1)) = label(1) ; 

        % update input
        params.P = params.Poff*ones(params.Ns,1) ; params.P(simlabel(1))=params.Pon ; 

        % Initialize time series
        % x0 = fsolve(@(x) wc_ode(0,x,pars),rand(2*Ns,1)) ; 
        x = zeros(2*params.Ns,length(tpoints)) ; x(:,1) = 300*rand(2*params.Ns,1) ; % x0 ; 
%         n = zeros(params.Nr,length(tpoints)) ; n(:,1) = rand(params.Nr,1) ; 


        % Simulate
        msg = [] ;
        h = 1000*params.dt ; 
        for i = 1:Nt-1

            if ~mod(i,5000)
                fprintf(repmat('\b',1,length(msg)))
                msg = sprintf('Simulating: %.2f percent',100*i/(Nt-1)) ; 
                fprintf(msg)
            end

            % check if the state has transitioned
            params.P = params.Poff*ones(params.Ns,1) ; params.P(simlabel(i))=params.Pon ;  

            % Calculate next step
            eta = randn(2*params.Ns,1) ; 
            f0 = wc_ode(tpoints(i),x(:,i),params) ; 
            k1 = x(:,i) + f0*h + sqrt(h)*params.sigma.*eta ;
            f1 = wc_ode(tpoints(i)+h,k1,params) ; 
            x(:,i+1) = x(:,i) + (h/2)*(f0+f1) + sqrt(h)*params.sigma.*eta ;

%             % Calculate noise step
%             eta = randn(params.Nr,1) ; 
%             f0 = wc_ode_E(tpoints(i),n(:,i),params) ; 
%             k1 = n(:,i) + f0*h + sqrt(h)*params.sigma.*eta ;
%             f1 = wc_ode_E(tpoints(i)+h,k1,params) ; 
%             n(:,i+1) = n(:,i) + (h/2)*(f0+f1) + sqrt(h)*params.sigma.*eta ;

            if sum(isnan(x(:,i+1)))
                return
            end

        end
        fprintf('\n')
        
        % Interpolate back onto original time axis
        x = interp1(tpoints'/1000,x',tlabel')' ; 

        % make LFP
        lfp = params.wplus*x(1:params.Ns,:) - x((params.Ns+1):(2*params.Ns),:) ; 
        lfp = lfp-mean(lfp,2) ; 

        % Make ROI time series
        wlfp = params.W*lfp ; % map components to ROIs/sensors
        nlfp = microstate.functions.pinknoise(size(wlfp'))' ; nlfp = nlfp-mean(nlfp,2) ; 
        lam = sqrt( sum(var(wlfp')) / ( sum(var(nlfp'))*params.SNR^2 ) ) ; % calculate lambda
        sim = wlfp+lam*nlfp ; sim = sim' ; % add noise
        
        


    end


    % ODE functions ---------
    function xdot = wc_ode(t,x,params)

        E = x(1:params.Ns) ; 
        I = x((params.Ns+1):(2*params.Ns)) ; 

        xdot = [(1./params.tauE).*(-E + params.f(params.P + params.wplus*E - I)) ; 
                (1./params.tauI).*(-I + params.f(params.wI*E))] ; 
    end

    function xdot = wc_ode_E(t,E,params)

        I0 = params.f(params.wI*E) ; 
        xdot = (1/mean(params.tauE))*(-E + params.f(params.Poff + params.wplus*E - I0)) ; 

    end

    
    
    

end