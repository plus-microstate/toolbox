function obj = stats_markov(obj)
% Calculate Markov matrix and Markov statistics
    % Check inputs
    if isempty(obj.label)
        error('To calculation Markov property, property label is required')
    end
    
    % Number of states and number of samples
    % Number of states
    label = obj.label ; 
    if ~isempty(obj.maps)
        Ns = size(obj.maps,2) ; 
    else
        Ns = max(label) ; 
    end
    Nt = length(label) ; 
    
    
    % Save syntax matrix
    T = zeros(Ns) ; 
    v1 = obj.label(2:end) ; v0 = obj.label(1:end-1) ;
    for i = 1:Ns
        for j = 1:Ns
            T(i,j) = sum((v0==i)&(v1==j)) ; 
        end
    end
    T = T./sum(T,2) ; % normalize rows to unity
    obj.stats.markov.matrix = T ; 
    
    % Deal with missing states
    for i = 1:Ns
        Px(i) = sum(label == i)/Nt ; 
    end
    idx0 = find(Px == 0) ; 
    orgnum = 1:Ns ; orgnum(idx0) = []; 
    
    Ns = Ns-length(idx0) ; 
    Px(idx0) = [] ; 
    
    oldlabel = label ; 
    for i = 1:Ns
        label(oldlabel==orgnum(i)) = i ; 
    end
    
    % 0th order Markov
    fij = zeros(Ns) ; % initialize transitions label(t)->label(t+1)
    fifj = zeros(Ns) ; % initialize number of label(t) * number of label(t+1)
    v0 = label(1:end-1) ; % label(t)
    v1 = label(2:end) ; % label(t+1)
    % loop over states i,j
    for i = 1:Ns 
        for j = 1:Ns
            fij(i,j) = sum((v1==i)&(v0==j)) ; 
            fifj(i,j) = sum(v0==i)*sum(v1==j) ; 
        end
    end
    fij(fij==0) = eps ; 
    G0 = 2*sum(sum(fij.*log(Nt*fij./fifj))) ; % This is the statistic for zero order markov
    p0 = 1-chi2cdf(G0,(Ns-1)^2) ; % null is distributed as chi2 function with (ns-1)^2 d.o.f., where ns is number of microstates


    % 1st order Markov
    fijk = zeros(Ns,Ns,Ns) ; % initialize transitions label(t-1)->label(t)->label(t+1)
    fj = zeros(Ns,Ns,Ns)  ; % initialize number of label(t)
    fij = zeros(Ns,Ns,Ns) ; % initialize number of label(t-1) * number of label(t)
    fjk = zeros(Ns,Ns,Ns)  ; % initialize number of label(t) * number of label(t+1)
    v0 = label(2:end-1) ; % label(t)
    v1 = label(3:end) ; % label(t+1)
    vm1 = label(1:end-2) ; % label(t-1)
    for i = 1:Ns
        for j = 1:Ns
            for k = 1:Ns
                fijk(i,j,k) = sum((v1==i)&(v0==j)&(vm1==k)) ; 
                fj(i,j,k) = sum(v0==j) ; 
                fij(i,j,k) = sum((v1==i)&(v0==j)) ; 
                fjk(i,j,k) = sum((v0==j)&(vm1==k)) ; 
            end
        end
    end
    fijk(fijk==0) = eps ; 
    fij(fij==0) = eps ; 
    fjk(fjk==0) = eps ; 
    G1 = 2*nansum(nansum(nansum(fijk.*log((fijk.*fj)./(fij.*fjk))))) ; 
    p1 = 1-chi2cdf(G1,Ns*(Ns-1)^2) ; % null is distributed as chi2 function with ns*(ns-1)^2 d.o.f.

    
    % Put into object
    obj.stats.markov.G0 = G0 ; 
    obj.stats.markov.p0 = p0 ; 
    obj.stats.markov.G1 = G1 ; 
    obj.stats.markov.p1 = p1 ; 
    
    % Append process
    obj = microstate.functions.process_append(obj,'Calculated statistics: Markov properties') ; 
    
 
end


    