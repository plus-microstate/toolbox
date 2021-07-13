function obj = stats_syntax(obj)
% Calculate syntax matrix
    % Check inputs
    if isempty(obj.label)
        error('To calculation microstate syntax properties, property label is required')
    end
    
    % Get microstate transitioning sequence
    transition = diff(obj.label) ~= 0 ; 
    label = [obj.label(transition),obj.label(end)] ; 
    
    % Number of states
    if ~isempty(obj.maps)
        Ns = size(obj.maps,2) ; 
    else
        Ns = max(obj.label) ; 
    end
    
    % Total number of microstates observed
    N = length(label) ; 
    
    % Find probability of each state
    for i = 1:Ns
        Px(i) = sum(label == i)/N ; 
    end
   
    % Save syntax matrix
    T = zeros(Ns) ; 
    v1 = obj.label(2:end) ; v0 = obj.label(1:end-1) ; 
    for i = 1:Ns
        for j = 1:Ns
            T(i,j) = sum((v0==i)&(v1==j)) ; 
        end
    end
    S = T-diag(diag(T)) ; % ignore self transitions in syntax matrix
    T = T./sum(T,2) ; % normalize rows to unity
    S = S./sum(S,2) ; % normalize rows to unity
    
    obj.stats.syntax.matrix = S ; 
    
    % Deal with case of zeros
    idx0 = find(Px == 0) ; 
    orgnum = 1:Ns ; orgnum(idx0) = []; 
    
    Ns = Ns-length(idx0) ; 
    Px(idx0) = [] ; 
    
    oldlabel = label ; 
    for i = 1:Ns
        label(oldlabel==orgnum(i)) = i ; 
    end
    
    % Make randomly distributed syntax matrix (P) and true syntax matrix
    % (S) - NOTE: we do not use the syntax matrix from
    % microstate.calculate_syntax, as this normalizes rows to unity, while
    % here we normalize the total matrix to unity
    clear P S
    v0 = label(1:end-1) ; v1 = label(2:end) ; 
    for i = 1:Ns
        for j = 1:Ns
            P(i,j) = Px(i)*Px(j)/(1-Px(i)) ; 
            S(i,j) = sum((v0==i)&(v1==j))/(N-1) ; 
        end
    end
    P = P-diag(diag(P)) ; % set diagonals to zero
    
    % Vectorize both matrices
    for i = 1:Ns ; P(i,i) = nan ; S(i,i) = nan ; end
    isgood = ~isnan(P) & ~isnan(S) ; 
    P = P(isgood) ; S = S(isgood) ; 
    
    % if either of these values are zero, they will give an error in
    % permutation testing, so remove zero values
    iszero = ~all([P,S],2); 
    P(iszero) = [] ; S(iszero) = [] ; 
    
    % Calculate chi2 distance between P and S
    chi2 = sum(((P-S).^2)./P) ; 
    
    % Run permutations
    Nperms = 4999 ; 
    chi2s = nan(Nperms,1) ; 
    PS = [P,S] ;
    for i = 1:Nperms
        
        % randomly reorder 'empirical' and 'random' labels
        for j = 1:length(PS)
            ind = randperm(2) ; 
            PSs(j,:) = PS(j,ind) ; 
        end
        Ps = PSs(:,1) ; Ss = PSs(:,2) ; 
        
        % Calculate chi2 distance between Ps and Ss
        chi2s(i) = nansum(((Ps-Ss).^2)./Ps) ; 
        
    end
    
    chi2s(i+1) = chi2 ; 
    p = sum(chi2s>=chi2)/length(chi2s) ; 
    
    % Update object
    obj.stats.syntax.chi2_random = chi2 ; 
    obj.stats.syntax.p_random = p ; 
    
    % Append process
    obj = microstate.functions.process_append(obj,'Calculated statistic: syntax matrix randomness') ; 
 
end


    