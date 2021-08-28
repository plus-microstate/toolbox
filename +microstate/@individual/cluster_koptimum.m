function [objopt,kopt,kvec,maps,gev,label] = cluster_koptimum(obj,varargin)
% Calculate optimum number of microstate classes 
% 
% Useage 
% Inputs: 
% - xpeak: EEG maps, identified via global field potential peaks using the
%      function gfppeaks. Alternatively, this can be full EEG time series
%      with options.findMaps = true. NxP, where N is number of points (maps
%      or time points) and P is dimensionality of data (number of
%      channels). 
% Outputs:
% - k: optimum number of microstates
% - copt: kxP matrix of cluster centroids (template maps). 
% - ropt: Variance explained by these centroids. 
% 
    options = microstate.functions.make_options(varargin) ; 
    
    % default options
    defaults = {'clustermethod','kmeans'; 
                'criterion','kneedle' ; % 'KrzanowskiLai','CalinskiHarabasz','DaviesBouldin','gap','silhoutte','kneedle','CrossValidationIndex','hmmFreeEnergy'
                'findpeaks',true ; 
                'kmin',2 ; 
                'kmax',20 ; 
                } ;  
    % ensure that if either hmm or hmmFreeEnergy are input, the other is a
    % default
    if ~isfield(options,'clustermethod') && isfield(options,'criterion')
        if strcmp(options.criterion,'hmmFreeEnergy')
            defaults{1,2} = 'hmm' ; 
        end
    end
    if ~isfield(options,'criterion') && isfield(options,'clustermethod')
        if strcmp(options.clustermethod,'hmm')
            defaults{2,2} = 'hmmFreeEnergy' ; 
        end
    end
    options = microstate.functions.add_options(options,defaults) ; clear defaults
    options.ClusterLabel = true ; % force true, as required for later functioning

    if strcmp(options.criterion,'hmmFreeEnergy') && ~strcmp(options.clustermethod,'hmm')
        error('Criterion hmmFreeEnergy only available for method hmm')
    end
            
    
    % make vector of k values to test
    if ~isfield(options,'kvec')
        options.kvec = options.kmin:options.kmax ; 
    else
        options.kmax = max(options.kvec) ; 
        options.kmin = min(options.kvec) ; 
    end
    kvec = options.kvec ; 
    
    % do clustering
    fprintf('Clustering microstates (%d to %d): \n',options.kmin,options.kmax)
    switch options.criterion
        case {'CalinskiHarabasz','DaviesBouldin','gap','silhouette'}
            label = nan(size(obj.data,1),options.kmax) ; % make a vector of cluster indices for each k value
            badk = [] ;
            for k = options.kvec
                out = obj.cluster_estimatemaps(k,options) ; 
                maps{k} = out.maps ; 
                gev(k) = out.stats.gev ; 
                label(:,k) = out.label ; 
                if length(unique(label(:,k))) < k
                    % warning('Only %d components found when setting k=%d',max(ind(:,k)),k)
                    badk = [badk,k] ; 
                end 
            end
            % remove bad k values
            maps(badk) = [] ; gev(badk) = [] ; label(:,badk) = [] ; 
            badk = sort(badk,'descend') ; 
            for ik = 1:length(badk)
                options.kvec(options.kvec == badk(ik)) = [] ; 
                options.kvec(options.kvec > badk(ik)) = options.kvec(options.kvec>badk(ik))-1 ; 
                kvec(kvec == badk(ik)) = [] ; 
            end 
            
            label = label(:,options.kvec) ;
            
            % get the distance function
            switch obj.modality
                case 'eeg'
                    distfun = @(x,c) 1-abs( ((x-mean(x,2))./vecnorm((x-mean(x,2)),2,2)) * ((c-mean(c,2))./vecnorm((c-mean(c,2)),2,2))' )' ; 
                case 'meg'
                    distfun = @(x,c) 1-abs( (x./vecnorm(x,2,2)) * (c./vecnorm(c,2,2))' )' ; 
                case 'source'
                    distfun = @(x,c) 1-( (x.^2./vecnorm(x.^2,2,2)) * (c.^2./vecnorm(c.^2,2,2))' )' ;
                case 'ampenv'
                    distfun = @(x,c) 1-( (x./vecnorm(x,2,2)) * (c./vecnorm(c,2,2))' )' ; 
            end

            eva = evalclusters(obj.data,label,options.criterion,'Distance',distfun) ; 
            
            kopt = eva.OptimalK ; 
            objopt = obj ; 
            objopt.maps = maps{kopt} ; 
            objopt.stats.gev = gev(kopt) ; 
            objopt.label = label(:,kopt)' ; 
            maps = maps(options.kvec) ; 
            gev = gev(options.kvec) ; 
            label = label(:,options.kvec) ; 
            
        case 'KrzanowskiLai' % NEED TO UPDATE
        
            w = nan(1,options.kmax+1) ; % w(k) is the error given by k clusters
            if options.kmin == 1 % since KL criterion depends on w(k-1), if kmin = 1 then this well be ignored as it is invalid.
                warning('Krzanowski-Lai criterion not valid for 1 cluster, using kmin = 2')
                options.kmin = 2 ; 
            end

            options.ErrorVal = true ; % output error in microstate_clustertopogs

            % Calculate w for kmin-1 and kmax+1, since KL criterion depends on
            % w(k-1) and w(k+1)
            p = size(obj.data,2) ; 
            for k = (options.kmin-1):(options.kmax+1)
                fprintf('k=%d\n',k)
                out = obj.cluster_estimatemaps(k,options) ;
                
                maps{k} = out.maps ; 
                gev(k) = out.stats.gev ; 
                label(:,k) = out.label ; 
                if length(unique(label(:,k))) < k
                    % warning('Only %d components found when setting k=%d',max(ind(:,k)),k)
                    badk = [badk,k] ; 
                end 
                
                D = zeros(k,1) ; n = zeros(k,1) ; 
                for r = 1:k
                    idx = find(out.label == r) ; 

                    % calculate D(r)
                    Dall = microstate.functions.DISS(obj.data(idx,:),obj.data(idx,:),obj.modality) ;
                    Dall = triu(Dall,1) ; 
                    D(r) = sum(Dall(:)) ; 

                    n(r) = length(idx) ; 

                end
                
                W(k) = sum((1./(2*n)).*D) ;
                M(k) = W(k)*(k^(2/p)) ; 
            end
            
            d = M(1:end-1) - M(2:end) ; % d(1-kmax)
            KL = (d(1:end-1)-d(2:end))./(M(1:end-2)) ; 
            
            [~,kopt] = max(KL) ; kopt = kvec(kopt) ; 
                
            gev = gev(options.kvec) ; maps = maps(options.kvec) ; label = label(:,options.kvec) ; 
            
            objopt = obj ; 
            objopt.maps = maps{kvec == kopt} ; 
            objopt.stats.gev = gev(kvec==kopt) ; 
            objopt.label = label(:,kvec==kopt)' ; 
            
        
        case 'kneedle'
            gev = nan(options.kmax,1) ; 
            label = nan(size(obj.data,1),options.kmax) ; % make a vector of cluster indices for each k value
            badk = [] ; 
            for k = options.kvec
                fprintf('k=%d\n',k)
                out = obj.cluster_estimatemaps(k,options) ; 
                maps{k} = out.maps ; 
                gev(k) = out.stats.gev ; 
                label(:,k) = out.label ; 
                if length(unique(label(:,k))) < k
                    % warning('Only %d components found when setting k=%d',max(ind(:,k)),k)
                    badk = [badk,k] ; 
                end 
            end 
            gev = gev(options.kvec) ; maps = maps(options.kvec) ; label = label(:,options.kvec) ; 
            kopt = kneedle(gev,options.kvec,options.kvec) ;  
           
            objopt = obj ; 
            objopt.maps = maps{kvec == kopt} ; 
            objopt.stats.gev = gev(kvec==kopt) ; 
            objopt.label = label(:,kvec==kopt)' ; 
            
        case 'CrossValidationIndex' % NEED TO UPDATE
            gev = nan(options.kmax,1) ; % make a vector of cluster indices for each k value
            for k = options.kvec
                fprintf('k=%d\n',k)
                out = obj.cluster_estimatemaps(k,options) ;
                maps{k} = out.maps ; 
                gev(k) = out.stats.gev ; 
                label(:,k) = out.label ; 
            end 
            gev = gev(options.kvec) ; maps = maps(options.kvec) ; label = label(:,options.kvec) ; 
            p = size(obj.data,2) ; 
            cvi = gev.*((p-1)./(p-1-options.kvec)').^2 ; 
            [~,k] = max(cvi) ;  
            kopt = options.kvec(k) ; 
            
            objopt = obj ; 
            objopt.maps = maps{kvec == kopt} ; 
            objopt.stats.gev = gev(kvec==kopt) ; 
            objopt.label = label(:,kvec==kopt)' ;
            
        case 'hmmFreeEnergy' % need to update
            
            label = nan(size(obj.data,1),options.kmax) ; % make a vector of cluster indices for each k value
            fe = nan(options.kmax,1) ; 
            badk = [] ;
            options.OutputFreeEnergy = true ; 
            for k = options.kvec
                fprintf('k=%d\n',k)
                [out,fe1] = obj.cluster_estimatemaps(k,options) ; 
                fe(k) = fe1{1} ; 
                maps{k} = out.maps ; 
                gev(k) = out.stats.gev ; 
                label(:,k) = out.label ; 
                if length(unique(label(:,k))) < k
                    % warning('Only %d components found when setting k=%d',max(ind(:,k)),k)
                    badk = [badk,k] ; 
                end 
            end
            % remove bad k values
            maps(badk) = [] ; gev(badk) = [] ; label(:,badk) = [] ; fe(badk) = [] ; 
            badk = sort(badk,'descend') ; 
            for ik = 1:length(badk)
                options.kvec(options.kvec == badk(ik)) = [] ; 
                options.kvec(options.kvec > badk(ik)) = options.kvec(options.kvec>badk(ik))-1 ; 
                kvec(kvec == badk(ik)) = [] ; 
            end 
            
            label = label(:,options.kvec) ; fe = fe(options.kvec) ; 
            [~,kopt] = min(fe) ; kopt = options.kvec(kopt) ; 
            maps = maps(options.kvec) ; 
            
            objopt = obj ; 
            objopt.maps = maps{kvec == kopt} ; 
            objopt.stats.gev = gev(kvec==kopt) ; 
            objopt.label = label(:,kvec==kopt)' ; 
            
    end
    
    obj = microstate.functions.process_append(obj,'Calculated optimum number of maps',options) ;
           

end
% --------------- END K OPTIMUM -------------------------------------------




% --------------- KNEEDLE -------------------------------------------------
function Ncomp = kneedle(px,k,kq)
    
    if nargin < 2
        k = (1:length(px))' ; 
    end
    if size(k,2) >1
        k = k' ; 
    end
    if nargin < 3
        kq = (min(k):max(k))' ; 
    end
    if size(kq,2) > 1
        kq = kq' ; 
    end
    % validateattributes(px,{'numeric'},{'column','nondecreasing','real'})
    validateattributes(k,{'numeric'},{'column','nondecreasing','real'})
    % determine number of components using kneedle algorithm
    px = fit(k,px,'smoothingspline') ; px = px(kq) ; % smoothed with spline
    pt = (kq-min(k))/(max(k)-min(k)) ; px = (px-min(px))/(max(px)-min(px)) ; % unit square
    dx = px-pt ; % difference between x and y
    [pks,Ncomp] = findpeaks(dx,kq) ;
    if length(pks)>1
        warning('multiple knees found, taking one with maximum difference')
        [pks,ni] = max(pks) ; 
        Ncomp = Ncomp(ni) ; 
    end
    
    
end
% ------------------ END KNEEDLE ------------------------------------------