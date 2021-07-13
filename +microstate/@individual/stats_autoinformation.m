function obj = stats_autoinformation(obj,maxlag)
% Calculate autoinformation of microstate sequence
    % Check inputs
    if isempty(obj.label)
        error('To calculation microstate autoinformation, property label is required')
    end
    if nargin < 2
        maxlag = 500 ; % 500 samples
    end
    if maxlag ~= round(maxlag) || maxlag <= 0
        error('Maximum lag must be a positive integer, corresponding to a number of samples')
    end
    
    % Calculate autoinformation
    aif = zeros(1,maxlag) ; 
    for i = 1:maxlag
        X = obj.label(1:(end-i)) ; 
        Y = obj.label((1+i):end) ; 
        
        aif(i) = MI(X,Y)/MI(Y,Y) ; 
    end
    
    % Put into object
    obj.stats.autoinformation = aif ; % update output
    
    % Append to process
    options = struct ; 
    options.maxlag = maxlag ; 
    obj = microstate.functions.process_append(obj,'Calculated statistic: autoinformation function',options) ; 
    
    
    % Mutual information function
    function mi = MI(X,Y)

        % Calculate H(Y)
        HY = 0 ;  
        for y = 1:max(max(X,Y))
            Py = sum(Y==y)/length(Y) ;
            if Py ~= 0
                HY = HY - Py*log2(Py) ; 
            end
        end

        % Calculate H(Y|X=x) and Px for all values of x
        HYX = 0 ; 
        for x = 1:max(max(X,Y))
            Px(x) = sum(X==x)/length(X) ; 
            idx = find(X==x) ; 
            Yx = Y(idx) ; % Y|X=x
            HYXx = 0 ; 
            for y = 1:max(max(X,Y))
                Pyx = sum(Yx == y)/length(Yx) ;
                if Pyx ~= 0 
                    HYXx = HYXx+Pyx*log2(Pyx) ;
                end
            end
            HYX = HYX-Px(x)*HYXx ; 
        end

        mi = HY-HYX ; 
    end
 
end


    