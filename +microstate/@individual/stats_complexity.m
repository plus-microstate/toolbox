function obj = stats_complexity(obj,Ntransitions)
% Calculate microstate transition complexity
    % 2nd input (Ntransitions): 
    % -1: use full sequence instead of transitioning sequence
    %  0: use all transitions
    % integer > 0: use this number of transitions
    
    % Check 2nd input
    if nargin < 2
        Ntransitions = 0 ; % use all transitions
    end
    if Ntransitions < -1 || round(Ntransitions) ~= Ntransitions % ensure it is an integer >= -1
        error('The 2nd input Ntransitions must be an integer >= -1')
    end


    % Check inputs
    if isempty(obj.label)
        error('To calculation microstate complexity, property label is required')
    end
    
    % Get microstate transitioning sequence
    if Ntransitions >= 0 % use transitioning sequence
        transition = diff(obj.label) ~= 0 ; 
        label = [obj.label(transition),obj.label(end)] ; 
        if Ntransitions > 0
            label = label(1:Ntransitions) ; 
        end
    else % -1: use full sequence
        label = obj.label ; 
    end
    
    % Calculate lzc
    C = lzc(label); 
    
    % save to microstate object
    obj.stats.complexity.complexity_raw = C ; 
    
    % Calculate normalized complexity
    a = length(unique(label)) ; 
    n = length(label) ; 
    upperbound = n/loga(n,a) ; 
    Cn = C/upperbound ; 
    obj.stats.complexity.complexity = Cn ; 
    
%     if surrflag
%         % Calculate lzc of random sequences
%         nsurrogate = 500 ; 
%         msg = [] ; 
%         pct_old = 0 ; 
%         for i = 1:nsurrogate
%             pct_new = floor(100*i/nsurrogate) ; 
%             if pct_new>pct_old
%                 fprintf(repmat('\b',1,length(msg))) ; 
%                 msg = sprintf('LZC surrogate distribution: %d%%%%',pct_new) ; 
%                 fprintf(msg)
%                 msg = msg(1:end-1) ;  
%                 pct_old = pct_new ; 
%             end
%             lbl2 = label(randperm(length(label))) ; 
%             Cs(i) = lzc(lbl2) ; 
%         end
%         fprintf(repmat('\b',1,length(msg))) ;
%         % z-score
%         Z = (C-mean(Cs))/std(Cs) ; 
%         p = 1-2*abs((sum(C>Cs)-(nsurrogate/2))/nsurrogate) ; % two-tailed test
%         
%         % save to microstate object
%         obj.stats.complexity.zscore = Z ; 
%         obj.stats.complexity.p = p ; 
%     end
%     
    
    % Append to process
    options = struct ; 
    options.sequence_length = length(label) ; 
    obj = microstate.functions.process_append(obj,'Calculated statistic: microstate complexity',options) ; 
    
end

function C = lzc(label)

    % Calculate LZC
    pointer = 1 ; % pointer (start of prefix)
    prefixlen = 1 ; % length of prefix
    i = 1 ; % letter in the sequence
    C = 1 ; % Lempel-Ziv complexity
    imax = i ; % how far through prefix are we when something changes
    while (prefixlen+i) <= length(label) % this ensures we keep going until the whole sequence is studied
        
        % check if the i'th value of the prefix matches the i'th value
        % after the prefix
        
        if label(pointer-1+i) == label(prefixlen+i) % if they match...
            i = i+1 ; % move onto the next letter in prefix
            
            
        else % if they don't match...
            imax = max(i,imax) ; % update imax, which says how many letters it matched for
            pointer = pointer+1 ; % shift the pointer across one
            
            if pointer == prefixlen+1 % if the new pointer is after the prefix...
                C = C+1 ; % increase complexity, as a new sequence is found
                prefixlen = prefixlen+imax ; % increase the length of the prefix to the new point
                i = 1 ; % reset pointer length
                pointer = 1 ; % reset pointer
                imax = i ; % reset how far through prefix we are, as we have a new prefix
                
            else % if the new pointer is in the prefix
                i = 1 ; % reset pointer length, and start from this pointer
            end
        end
    end
    if i ~= 1 % the final string used hasn't been counted, so we need to add that
        C = C+1 ; 
    end
    
end

function loga_n = loga(n,a)

    loga_n = log(n)/log(a) ; 
    
end