function obj = stats_mutualinformation(obj,miseq,atpeaks)
% Calculate mutual information between microstate sequence and another sequence (2nd input)
    if nargin < 3
        atpeaks = false ; 
    end

    % Check inputs
    if isempty(obj.label)
        error('To calculate microstate mutual information, property label is required')
    end
    if nargin < 2
        error('To calculate microstate mutual information, the 2nd input is required. This is the sequence to calculate MI against')
    end

    % get GFP peaks
    if atpeaks
        [~,ind] = obj.cluster_get_gfppeaks ; 
        origseq = obj.label(ind) ; 
        miseq = miseq(ind) ; 
    else
        origseq = obj.label ; 
    end

    % Put into object
    obj.stats.mutualinformation = MI(origseq,miseq)/MI(miseq,miseq) ; % update output
    
    % Append to process
    options = struct ; 
    options.miseq = miseq ; 
    obj = microstate.functions.process_append(obj,'Calculated statistic: mutual information against input signal given in options',options) ; 
    
    
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


    