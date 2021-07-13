function [obj,additionalout] = cluster_global(obj,k,varargin)
% Perform group level microstate clustering

    options = microstate.functions.make_options(varargin) ; 

    % default options
    defaults = {'clustermethod','kmeans'; % kmeans, aahc, pca, ica, hmm
                'kmeans_maxiter',100 ; 
                'kmeans_replicates',20 ;  
                'hmm',struct ; 
                'cohortstat','data' ; % data or maps
                } ;  
    options = microstate.functions.add_options(options,defaults) ; clear defaults
    
    if ~isfield(options,'findpeaks')
        fullsample = false(length(obj.individual),1) ; 
        for i = 1:length(obj.individual)
            if length(obj.individual(i).sample) == length(obj.individual(i).time)
                fullsample(i) = true ; 
            end
        end
        if strcmp(options.clustermethod,'kmeans') && all(fullsample)
            options.findpeaks = true ; 
        else
            options.findpeaks = false ; 
        end
    end
    
    % validate options
    validateattributes(options.kmeans_maxiter,{'numeric'},{'scalar','integer'},'microstate.clustertopogs','options.maxiter')
    validateattributes(options.kmeans_replicates,{'numeric'},{'scalar','integer'},'microstate.clustertopogs','options.replicates')
    validateattributes(options.findpeaks,{'logical'},{'scalar'},'microstate.clustertopogs','options.findpeaks')
    
    % initialize additionalout
    additionalout = cell(0) ; 
    
    
    %% Get data to cluster
    
    switch options.cohortstat
        case 'data'
            X = [] ; 
            for i = 1:length(obj.individual)
                if ~options.findpeaks
                    X = [X ; obj.individual(i).data] ; 
                else
                    obj.individual(i) = obj.individual(i).calculate_gfp ;
                    X = [X ; obj.individual(i).cluster_get_gfppeaks] ;
                end
            end
            
        case 'maps'
            X = [] ; 
            for i = 1:length(obj.individual)
                X = [X ; obj.individual(i).maps] ; 
            end
    end
    
    % make a new artificial individual to perform clustering
    newcoh = microstate.individual(X,obj.individual(1).modality,1:size(X,1)) ; 
    newopt = options ; newopt.findpeaks = false ; 
    
     % cluster
    [newcoh,additionalout] = newcoh.cluster_estimatemaps(k,newopt) ; 
    
    %% Update cohort object
    
    obj.globalmaps = newcoh.maps ; 
    
    str = sprintf('Globally Clustered to estimate %d maps',k) ; 
    obj = microstate.functions.process_append(obj,str,options) ;

   

end
% --------------- END CLUSTER TOPOGRAPHIES --------------------------------

