function [obj,additionalout] = cluster_estimatemaps(obj,k,varargin)
% Cluster to calculate a specified number of microstate maps 


    options = microstate.functions.make_options(varargin) ; 

    % default options
    defaults = {'clustermethod','kmeans'; % kmeans, aahc, pca, ica, hmm
                'kmeans_maxiter',100 ; 
                'kmeans_replicates',20 ; 
                'findpeaks',true ; 
                'nsample',1 ;
                'hmm',struct ; 
                'keep_polarity',false ; 
                } ;  
    options = microstate.functions.add_options(options,defaults) ; clear defaults
    if options.keep_polarity
        options.findpeaks = false ; 
    end
    
    % validate options
    validateattributes(options.kmeans_maxiter,{'numeric'},{'scalar','integer'},'microstate.clustertopogs','options.maxiter')
    validateattributes(options.kmeans_replicates,{'numeric'},{'scalar','integer'},'microstate.clustertopogs','options.replicates')
    validateattributes(options.findpeaks,{'logical'},{'scalar'},'microstate.clustertopogs','options.findpeaks')
    
    % initialize additionalout
    additionalout = cell(0) ; 
    
    
    %% kmeans clustering
    switch options.clustermethod
        case 'kmeans'
            
            % make sure we have a gfp
            obj = obj.calculate_gfp ; 

            % find peaks
            if options.findpeaks
                x = obj.cluster_get_gfppeaks(options.nsample) ; 
            else
                x = obj.data ; 
            end
            options.nsample = 1 ; % avoid accidentally sampling twice in future calls

            % Sizes of data
            [n,p] = size(x) ; 

            % get the cluster statistic, distance function, and centroid
            % function
            switch obj.modality
                case 'eeg'
                    x = x-mean(x,2) ; % re-reference to average
                    if options.keep_polarity
                        distfun = @(x,c) 1-( (x./vecnorm(x,2,2)) * ((c-mean(c,2))./vecnorm((c-mean(c,2)),2,2))' )' ;
                    else
                        distfun = @(x,c) 1-abs( (x./vecnorm(x,2,2)) * ((c-mean(c,2))./vecnorm((c-mean(c,2)),2,2))' )' ; 
                    end
                    centfun = @(v1) v1 ; 
                case 'meg'
                    % do nothing
                    if options.keep_polarity
                        distfun = @(x,c) 1-( (x./vecnorm(x,2,2)) * (c./vecnorm(c,2,2))' )' ; 
                    else
                        distfun = @(x,c) 1-abs( (x./vecnorm(x,2,2)) * (c./vecnorm(c,2,2))' )' ; 
                    end
                    centfun = @(v1) v1 ; 
                case 'source'
                    x = abs(x) ; % .^2 ; % set to magnitude
                    distfun = @(x,c) 1-( (x./vecnorm(x,2,2)) * (c./vecnorm(c,2,2))' )' ;
                    centfun = @(v1) abs(v1) ; 
                    if options.keep_polarity
                        warning('For source/amplitude envelope data, keep_polarity has no effect and polarity is always ignored.') ; 
                    end
                case 'ampenv'
                    % do nothing
                    distfun = @(x,c) 1-( (x./vecnorm(x,2,2)) * (c./vecnorm(c,2,2))' )' ; 
                    centfun = @(v1) abs(v1) ; 
                    if options.keep_polarity
                        warning('For source/amplitude envelope data, keep_polarity has no effect and polarity is always ignored.') ; 
                    end
            end
    
            % Initialize explained variance
            gev = 0 ; 
            maps = nan(k,p) ;
            
            % Make figure to update progress
            progfig = figure('Name','k-means progress','numbertitle','off') ; 
            set(progfig,'Units','Normalized') ; 
            set(progfig,'Position',[0.4,0.25,0.2,0.5]) ; 
            set(progfig,'MenuBar','none','Toolbar','none') ;
            set(progfig,'Color','w') ; 
            progax = axes ; 
            set(progax,'Units','Normalized') ; 
            set(progax,'Position',[0.45,0.1,0.1,0.8]) ;
            title(sprintf('k-means progress (k=%d)',k))
            
            cla
            bar(progax,1,0,'r','LineStyle','none') ; 
            ylim([0,options.kmeans_maxiter*options.kmeans_replicates])
            xlim([0.6 1.4])
            set(progax,'XTick',[],'YTick',options.kmeans_maxiter*(1:options.kmeans_replicates))
            set(progax,'YTickLabels',1:options.kmeans_replicates)
            grid on
            set(progax,'GridAlpha',1)
            ylabel('# replicates complete')
            set(progax,'FontSize',12)
            drawnow
    
            % Loop over replicates
            for rep = 1:options.kmeans_replicates % loop over replications, initializing, running algorithm, and calculating explained variance

                progcount = options.kmeans_maxiter*(rep-1) ; 
                nconv = 0 ; % initialize number converged

                % initialization -  k-means++
                c = zeros(k,p) ; 
                c(1,:) = x(randi(n),:) ; % centroid of first cluster, chosen at random
                for j = 2:k % loop over classes, estimating initial seeds using k-means++ algorithm
                    d = distfun(x,c(1:(j-1),:)) ; 
                    dmin = min(d,[],1) ; % distance of each point from nearest centroid
                    px = dmin.^2/sum(dmin.^2) ; % probability of each point being selected
                    px = cumsum(px) ; % cumulative sum of probabilities for drawing a number
                    xj = find(rand<px) ; xj = xj(1) ; % drawn from distribution px
                    c(j,:) = x(xj,:)./vecnorm(x(xj,:)) ; % j'th centroid
                end
                d = distfun(x,c) ; [d,idx] = min(d,[],1) ;  % distance of each point from nearest centroid & initial class assignments

                
                % run k-means clustering algorithm
                iter = 0 ; % iterations calculated
                while iter < options.kmeans_maxiter % run k-means clustering algorithm until maximum number of iterations or convergence reached

                    err = false ; % make error flag false by default
                    iter = iter+1 ; % update number of iterations
                    d0 = d ; idx0 = idx ; % save previous iteration, for checking convergence
                    for j = 1:k % Loop over each cluster
                        Xj = x(idx0 == j,:) ; % Maps previously assigned to cluster j
                        if isempty(Xj) % Not found any in cluster
                            err = true ;  % error flag
                        else
                            [V,D] = eig(Xj'*Xj) ; 
                            [~,indmax] = max(diag(D)) ; % should just be the biggest value, but put this in in case
                            v1 = V(:,indmax) ; % get first eigenvector
                            if options.keep_polarity
                                C = v1'*Xj' ; 
                                v1 = sign(mean(C))*v1 ; 
                            end
                        end
                        if err
                            break
                        end
                        c(j,:) = centfun(v1) ; % assign as centroid
                    end
                    
                    if err % if error, break while loop
                        break
                    end
                    
                    d = distfun(x,c) ; % get distance from each point to each centroid
                    [d,idx] = min(d,[],1) ; % distance of each point from nearest centroid & class assignments
                    
%                     % Update figure
%                     progcount = progcount+1 ; 
%                     figure(progfig) ; 
%                     cla
%                     bar(progax,1,progcount,'r','LineStyle','none') ; 
%                     ylim([0,options.kmeans_maxiter*options.kmeans_replicates])
%                     xlim([0.6 1.4])
%                     set(progax,'XTick',[],'YTick',options.kmeans_maxiter*(1:options.kmeans_replicates))
%                     set(progax,'YTickLabels',1:options.kmeans_replicates)
%                     grid on
%                     set(progax,'GridAlpha',1)
%                     ylabel('# replicates complete')
%                     set(progax,'FontSize',12)
%                     drawnow
                    
                    if all(idx == idx0) % all have converged
                        break % no need to continue
                    end
                end
                
                progcount = options.kmeans_maxiter*rep ; 
                % figure(progfig) ; 
                cla
                bar(progax,1,progcount,'r','LineStyle','none') ; 
                ylim([0,options.kmeans_maxiter*options.kmeans_replicates])
                xlim([0.6 1.4])
                set(progax,'XTick',[],'YTick',options.kmeans_maxiter*(1:options.kmeans_replicates))
                set(progax,'YTickLabels',1:options.kmeans_replicates)
                grid on
                set(progax,'GridAlpha',1)
                ylabel('# replicates complete')
                set(progax,'FontSize',12)
                title(sprintf('k-means progress (k=%d)',k))
                drawnow
                
                if err % if error, move onto next replicate
                    continue
                end
                
%                 % if source, sqrt maps
%                 if strcmp(obj.modality,'source')
%                     c = sqrt(c) ; 
%                 end

                % calculate gev
                obj.maps = c' ; % put the maps into the object
                obj = obj.cluster_alignmaps(options.findpeaks,'keep_polarity',options.keep_polarity) ; % align maps to data
                obj = obj.stats_gev ; % get gev
                if obj.stats.gev > gev % if explained variance greater than previous maximum, replace
                    gev = obj.stats.gev ; % update maximum gev
                    maps = c ; % update maximum centroids
                end
            end
            
%             figure(progfig)
%             cla
%             clf(progfig), plot([],[])
            close(progfig) ; 
            
            % assign to object
            obj.maps = maps' ; % assign maps
            [obj,additionalout] = obj.cluster_alignmaps(options.findpeaks,'keep_polarity',options.keep_polarity) ; % align to data % RETEST TO CHECK!!!
            obj = obj.stats_gev ; % get gev
            
        case 'pca'
            
            % get maps 
            maps = pca(obj.data,'NumComponents',k)' ; 
            
            % calculate gev
            obj.maps = maps' ; % put the maps into the object
            obj = obj.cluster_alignmaps ; % align maps to data
            obj = obj.stats_gev ; % get gev
            
        case 'ica'
            
            % do ica
            path = microstate.functions.toolbox_path ;
            warning('Adding fastica to path')
            addpath(fullfile(path,'+external','FastICA_25')) ; 
            [~,~,maps,~] = evalc("fastica(obj.data','numOfIC',k)") ;
            
            % calculate gev
            obj.maps = maps ; % put the maps into the object
            obj = obj.cluster_alignmaps ; % align maps to data
            obj = obj.stats_gev ; % get gev
            
        case 'hmm'
            
            switch obj.modality
                case {'ampenv'}
                    % Set the HMM-MAR options based on example script from Baker et
                    % al. (2014) + run_HMMMAR.m (modality M/EEG power)
                    hmmopts = struct;
                    hmmopts.K = k;
                    hmmopts.Fs = 1/mean(diff(obj.time)) ;
                    hmmopts.covtype = 'full';
                    hmmopts.order = 0 ; 
                    hmmopts.zeromean = 0 ; 
                    hmmopts.standardise = 1;
                    hmmopts.verbose = 1;
                    hmmopts.onpower = 0 ; % already on power
                case {'meg','eeg','source'}
                    % Set the HMM-MAR options based on example script from Viduarre et
                    % al. (2018) Nat Comms and run_HMMMAR.m (modality M/EEG +
                    % ndim>10)
                    Fs = 1/mean(diff(obj.time)) ; % sampling rate
                    embeddedlag = floor(60e-3*Fs/2) ; % 60 ms window, rounded to the nearest odd number of samples (then subtract time zero and divide by 2)
                    hmmopts = struct;
                    hmmopts.order = 0;
                    hmmopts.zeromean = 1;
                    hmmopts.covtype = 'full';
                    hmmopts.embeddedlags = -embeddedlag:embeddedlag;
                    hmmopts.pca = size(obj.data,2)*2;
                    hmmopts.K = k;
                    hmmopts.Fs = Fs;
                    hmmopts.verbose = 1;
                    hmmopts.onpower = 0; 
                    hmmopts.standardise = 1;
                    hmmopts.standardise_pc = hmmopts.standardise;
            end
            
            % HMM computation
            path = microstate.functions.toolbox_path ;
            warning('Adding HMM-MAR to path')
            addpath(genpath(fullfile(path,'+external','HMM-MAR-master'))) ; 
            [hmm,~,~,viterbipath,~,~,fehist] = hmmmar(obj.data,size(obj.data,1),hmmopts);
            if any(strcmp(obj.modality,{'meg','eeg','source'}))
                viterbipath = [nan(embeddedlag,1) ; viterbipath ; nan(embeddedlag,1)] ;
            end
            
            % now get microstate GFP peaks
            [xpeaks,pks] = cluster_get_gfppeaks(obj) ; 
            
            % maps can come with centroid function
            switch obj.modality
                case 'eeg'
                    xpeaks = xpeaks-mean(xpeaks,2) ; % re-reference to average 
                    centfun = @(v1) v1 ; 
                case 'meg'
                    centfun = @(v1) v1 ; 
                case {'source','ampenv'}
                    xpeaks = abs(xpeaks) ; 
                    centfun = @(v1) abs(v1) ; 
            end
            
            % get maps
            idx = viterbipath(pks) ; 
            for j = 1:k % Loop over each cluster
                Xj = xpeaks(find(idx==j),:) ; 
                [V,D] = eig(Xj'*Xj) ; 
                [~,indmax] = max(diag(D)) ; % should just be the biggest value, but put this in in case
                v1 = V(:,indmax) ; 
                maps(j,:) = centfun(v1) ; 
            end
            
            % calculate gev
            obj.maps = maps' ; % put the maps into the object
            obj.label = viterbipath ; % align maps to data
            obj = obj.stats_gev ; % get gev
            
            if isfield(options,'OutputFreeEnergy')
                if options.OutputFreeEnergy
                    additionalout = [additionalout , fehist(end)] ; 
                end
            end
            
            if isfield(options,'OutputHMM')
                if options.OutputHMM
                    additionalout = [additionalout, hmm] ; 
                end
            end
            
            
            
    end
    
    str = sprintf('Clustered to estimate %d maps',k) ; 
    obj = microstate.functions.process_append(obj,str,options) ;


end
% --------------- END CLUSTER TOPOGRAPHIES --------------------------------

