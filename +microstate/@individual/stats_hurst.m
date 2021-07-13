function obj = stats_hurst(obj)
% Calculate Hurst exponent of microstate sequences
    % Check inputs
    if isempty(obj.label)
        error('To calculation microstate coverage, property label is required')
    end
    
    %% Convert microstate sequence to random walk
    
    % Assign states values of +1, -1, or 0
    Ns = max(obj.label) ; % number of states
    if ~mod(Ns,2) % even number of states
        class0 = [] ; % classes worth 0 in the walk
        classplus = randperm(Ns,Ns/2) ; % classes worth +1 in the walk
        classminus = setdiff(1:Ns,classplus) ; % classes worth -1 in the walk
    else % odd number of states
        classes = 1:Ns ; % make a list of the classes
        class0 = randi(length(classes)) ; % select a class to be worth 0 in the walk
        classes(class0) = [] ; % remove this class
        classplus = classes(randperm(Ns-1,(Ns-1)/2)) ; % classes worth +1 in the walk 
        classminus = setdiff(classes,classplus) ; % classes worth -1 in the walk
    end

    % Make a randomwalk
    rw = obj.label ; 
    % replace all classes in classminus with -1
    for i = 1:length(classminus)
        rw(obj.label==classminus(i)) = -1 ; 
    end
    % replace all classes in classplus with +1
    for i = 1:length(classplus)
        rw(obj.label==classplus(i)) = 1 ; 
    end
    % replace all classes in class0 with 0
    for i = 1:length(class0)
        rw(obj.label==class0(i)) = 0 ; 
    end
    % convert to randomwalk
    rw = cumsum(rw) ; 
    
    %% Calculate Hurst exponent using detrended fluctuation analysis
    
    % We will split the random walk into windows of length n samples. But
    % first we need to pick a range of n. Choose powers of 2 for
    % computational efficiency.  
    log2nmin = 6 ; % choose 64=2^6 as a minimum reliable estimate. 
    log2nmax = floor(log2(length(rw)))-2 ; % maximum n which would fit in the length of the randomwalk is floor(log2(length(rw))). Ensure we have >4 segments by subtracting 2. 
    n = 2.^(log2nmin:log2nmax) ; % get our window lengths
    
    if isempty(n)
        obj.stats.hurst = nan ; 
        return
    end
    
    % Loop over window lengths and calculate fluctuation F(n)
    rw = rw(:) ; % ensure it is a column
    for i=1:length(n) 
        
        % Construct piecewise linear sequence
        numwindows=floor(length(rw)/n(i)) ; % number of windows of length n(i) which fit into randomwalk
        Y = nan(length(rw),1) ; % will be filled in with piecewise straight line fits
        for j=1:numwindows % loop over windows
            idx = (n(i)*(j-1)+1):(n(i)*j) ; idx = idx' ; % indices of window
            rwi = rw(idx) ; % random walk in window
            coeff = [ones(n(i),1),idx]\rwi ; % fit linear model
            Y(idx) = coeff(1)+coeff(2)*idx ; % update with prediction from linear model
        end
        
        % calculate fluctuation 
        F(i)=sqrt(nanmean((rw-Y).^2));
        
    end
    
    % Calculate Hurst exponent as gradient of graph of log(n) vs log(F)
    coeff = [ones(length(n),1),log2(n(:))]\log2(F(:)) ; % fit linear model
    obj.stats.hurst = coeff(2) ; % update output
    
    % Append process
    obj = microstate.functions.process_append(obj,'Calculated statistic: Hurst exponent') ; 
 
end

    