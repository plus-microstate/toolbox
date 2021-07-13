function obj = preprocess_orthogonalize(obj,display)
% Orthogonalize time courses
    % Colclough et al (2015), NeuroImage 117:439-448
    % Input: 
    % - obj.data: Original data. mobj.datan, where m is number of samples, n number of ROIs
    % - display: Display progress. Can be: 
    %            'none', display no progress
    %            'iter', display iteration number in console (default)
    %            'full', display full output in console, at cost of speed
    %            'plot', display as a plot, at cost of speed
    % Output: 
    % - obj.data_orthog: Orthogonalized data. 
    
    % get size of data
    [m,n] = size(obj.data) ; 
    
    % default options
    if nargin < 2
        display = 'iter' ; 
    end

    % Check Z is full rank
    if rank(obj.data) < n
        warning('Data is not full rank')
    end

    % Initialize output
    if ~strcmp(display,'none')
        msg = [] ; 
        fprintf('\nOrthogonalizing data: ')
    end

    % Initialize D, O, and error
    if any(strcmp(display,{'iter','full'}))
        fprintf(repmat('\b',1,length(msg))) ; 
        msg = sprintf('Initializing ...') ; 
        fprintf(msg); 
    end 
    D = eye(n) ; 
    [U,S,V] = svd(obj.data*D,'econ') ; O = U*V' ; clear U S V
    ep0 = trace(obj.data'*obj.data) - 2*trace(obj.data'*O*D) + trace(D*D) ; 

    % Iterate equations (4), (6), (8) of Colclough et al. until convergence
    maobj.dataiter = 20 ; % maobj.dataimum number of iterations. Reported maobj.dataimum of 20 in paper for fMRI. 
    tol = 1e-8 ; % tolerance, placeholder
    iter = 0 ; 
    diff_ep = ep0 ;

    % Output initial iteration
    switch display
        case 'iter'
            fprintf(repmat('\b',1,length(msg))) ; 
            msg = sprintf('Iteration %d',iter) ; 
            fprintf(msg);
        case 'full'
            maxR = max(max(abs(corr(O*D)).*(1-eye(n)))) ; 
            fprintf(repmat('\b',1,length(msg))) ; 
            msg = sprintf('Iteration %d, maox abscorr %f, norm difference %f',iter,maxR,ep0);
            fprintf(msg);
            clear maxR
    end      

    while diff_ep > tol && iter<maobj.dataiter

        iter = iter+1; 

        % Eqn (4): SVD
        [U,S,V] = svd(obj.data*D,'econ') ; 

        % Eqn (6)
        O = U*V' ; 
        clear U S V

        % Eqn (8)
        D = diag(diag(obj.data'*O)) ; 

        % Error (difference from previous iteration)
        ep1 = trace(obj.data'*obj.data) - 2*trace(obj.data'*O*D) + trace(D*D) ;
        diff_ep = abs(ep1-ep0) ; 
        ep0 = ep1; 


        % Output
        switch display
            case 'iter'
                fprintf(repmat('\b',1,length(msg))) ; 
                msg = sprintf('Iteration %d',iter) ; 
                fprintf(msg);
            case 'full'
                maxR = maobj.data(max(abs(corr(O*D)).*(1-eye(n)))) ; 
                fprintf(repmat('\b',1,length(msg))) ; 
                msg = sprintf('Iteration %d, maobj.data abscorr %f, norm difference %f',iter,maxR,ep0) ;
                fprintf(msg);
                clear maxR
        end  

    end

    % Output
    switch display
        case {'iter','full'}
            fprintf(repmat('\b',1,length(msg))) ; 
            msg = sprintf('Done\n') ;
            fprintf(msg);
    end  

    clear Z diff_ep ep0 ep1 iter m maxR n tol % save memory
    obj.data = O*D ; 
    
    % update the process
    obj = microstate.functions.process_append(obj,'Orthogonalized data') ; 
    
    % update the gfp
    obj = obj.calculate_gfp() ; 
end
% --------------- END ORTHOGONALIZE TIME SERIES ---------------------------
