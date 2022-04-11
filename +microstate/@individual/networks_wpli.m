function [p,confusionmat,nets] = networks_wpli(obj,frq,keepstates,epochlength,staticflag) ; 
% Calculate microstate segmented functional connectivity

    % check we have GFP, maps, and labels
    if isempty(obj.data)
        error('Data is required to calculate microstate networks')
    end
    if isempty(obj.label)
        obj = obj.cluster_alignmaps ; 
    end
    
    % Extra inputs
    if nargin < 2
        frq = []; 
    end
    if nargin < 3 || isempty(keepstates)
        keepstates = false ; 
    end
    if nargin < 4 || isempty(epochlength)
        epochlength = 1280 ; % 5 seconds at 256 Hz
    end
    if nargin < 5 || isempty(staticflag)
        staticflag = false ; 
    end
    
    
    % --- SIMPLE CASE - STATIC NETWORK ---
    if staticflag
        % Calculate phase
        if ~isempty(frq)
            fobj = obj.preprocess_filter(frq(1),frq(2)) ; 
        else
            fobj = obj ; 
        end
        ph = angle(hilbert(fobj.data)) ; 
        z = abs(hilbert(fobj.data)) ; 
        clear fobj

        % Epoch the data
        window = 1:epochlength:size(ph,1) ; 
        window = [window(1:end-1)' , window(2:end)'-1] ; 

        % Loop over epochs 
        nets = [] ; 
        for j = 1:size(window,1)
            ind = window(j,1):window(j,2) ; 
            nets = cat(3,nets,wpli(ph(ind,:),z(ind,:))) ;
        end
        p = [] ; confusionmat = [] ; 
        nets = mean(nets,3) ; 
        return
    end
    
    % --- MICROSTATE SEGMENTED FC ---
    
    % Number of states
    if ~isempty(obj.maps)
        Nstates = size(obj.maps,2) ; 
    else
        Nstates = max(obj.label) ; 
    end
    
    % Get durations and their labels
    transition = find(diff(obj.label)) ; 
    transition = [0,transition,length(obj.label)] ; 
    states = [transition(1:end-1)'+1 , transition(2:end)'] ;
    dur = diff(states,[],2)+1 ; 
    
    % Remove states containing bad samples
    bad_samples = obj.bad_samples ;
    Nloop = 0 ; 
    while ~isempty(bad_samples)
        isbad = find((bad_samples(1) >= states(:,1)) & (bad_samples(1) <= states(:,2))) ; 

        indstate = find((bad_samples>=states(isbad,1))&(bad_samples<=states(isbad,2))) ; 
        bad_samples(indstate)=[] ;  
        
        states(isbad,:) = [] ; 
        dur(isbad,:) = [] ; 
        
        Nloop = Nloop+1 ; 
    end
    

    % Get labels of states
    lbl = [] ; 
    for i = 1:size(states,1)
        lbl(i,1) = obj.label(states(i,1)) ; 
    end
    
    % Calculate phase
    if ~isempty(frq)
        fobj = obj.preprocess_filter(frq(1),frq(2)) ; 
    else
        fobj = obj ; 
    end
    ph = angle(hilbert(fobj.data)) ; 
    z = abs(hilbert(fobj.data)) ; 
    clear fobj
    
    % Loop over microstates
    nets = cell(1,Nstates) ; 
    for i = 1:Nstates
        
        % Find which states are microstate i
        ind = find(lbl == i) ; 
        
        % Loop over these states and concatenate
        phi = [] ; zi = [] ; 
        for j = 1:length(ind) 
            state = states(ind(j),:) ; 
            phi = [phi ; ph(state(1):state(2),:)] ; 
            zi = [zi ; z(state(1):state(2),:)] ; 
%             nets{i} = cat(3,nets{i},wpli(ph(state(1):state(2),:),z(state(1):state(2),:))) ; 
        end
        
        % Epoch the data
        window = 1:epochlength:size(phi,1) ; 
        window = [window(1:end-1)' , window(2:end)'-1] ; 
        
        % Loop over epochs
        msg = [] ; 
        for j = 1:size(window,1)
            fprintf(repmat('\b',1,length(msg)))
            msg = sprintf('Calculating wPLI for microstate %d of %d, epoch %d of %d',i,Nstates,j,size(window,1)) ; 
            fprintf(msg) ; 
            ind = window(j,1):window(j,2) ; 
            nets{i} = cat(3,nets{i},wpli(phi(ind,:),zi(ind,:))) ;
        end
        fprintf(repmat('\b',1,length(msg)))
        
    end
    
    % MVPA
    p = [] ; confusionmat = [] ;
%     Nperms = 1000 ; 
%     [p,confusionmat] = microstate.functions.networks_mvpa(nets,Nperms) ; 
%     
    % If just wanting the average, take this average
    for i = 1:Nstates
        if ~keepstates
            nets{i} = mean(nets{i},3) ; 
        end
    end
    
end

function C = wpli(ph,z)
    % Vinck et al. (2011) NeuroImage 55:1547-1565
    % Colclough et al. (2016) NeuroImage 138:284-293
                     
    % calculate wPLI
    n = size(ph,2) ; 
    C = zeros(n) ; 
    for i = 1:n
        dphi = ph-ph(:,i) ; % repmat(x(:,i),1,n) ;
        z1z2 = (z.*z(:,i)).*sin(dphi) ; % repmat(z(:,i),1,n) ;
        C(i,:) = abs(mean(z1z2))./mean(abs(z1z2)) ; 
    end

end % wPLI