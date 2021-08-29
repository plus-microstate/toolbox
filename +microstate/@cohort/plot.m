function plt = plot(obj,param,varargin)
% Visualize microstate data and statistics
    if nargin == 1
        param = 'data' ; 
    end
    
    % ensure plotting functions are on path
    mspath = microstate.functions.toolbox_path ; 
    addpath(fullfile(mspath,'+external','brewermap'))
    
    switch param
            
        case {'mean_duration','hurst','complexity','complexity_Z','markov_G0','markov_G1','gfp_peaksfreq','gev','syntax_chi2_random'}
            
            % number of conditions
            numcondition = max(obj.condition) ;
            
            % concatenate data
            X = [] ; % initialize data
            G = [] ; % initialize groups
            for i = 1:numcondition
                switch param
                    case 'markov_G0'
                        X = [X;log10(obj.stats(i).markov.G0)] ; 
                        G = [G;i*ones(length(obj.stats(i).markov.G0),1)] ; 
                    case 'markov_G1'
                        X = [X;log10(obj.stats(i).markov.G1)] ; 
                        G = [G;i*ones(length(obj.stats(i).markov.G1),1)] ;
                    case 'complexity' 
                        X = [X;obj.stats(i).complexity.complexity] ; 
                        G = [G;i*ones(length(obj.stats(i).complexity.complexity),1)] ;
                    case 'complexity_raw'
                        X = [X;obj.stats(i).complexity.complexity_raw] ; 
                        G = [G;i*ones(length(obj.stats(i).complexity.complexity_raw),1)] ;
                    case 'syntax_chi2_random'
                        X = [X;obj.stats(i).syntax.chi2_random] ; 
                        G = [G;i*ones(length(obj.stats(i).syntax.chi2_random),1)] ;
                    otherwise
                        X = [X;obj.stats(i).(param)] ; 
                        G = [G;i*ones(length(obj.stats(i).(param)),1)] ; 
                end
            end
            
            if strcmp(param,'mean_duration')
                X = 1000*X ; % convert to ms
            end
              
            
            % make colours for the plot
            idx = find(strcmp(varargin,'cmap')) ; 
            if ~isempty(idx)
                cmap = varargin{idx+1} ; 
            else 
                if numcondition<=9
                    cmap = brewermap(numcondition,'Set1') ; 
                else
                    rng('default')
                    cmap = lhsdesign(numcondition,3) ; 
                end
            end
            
            % plot data
            plt = microstate.external.boxplot_LT(X,'Groups',G,'GroupColors',cmap,...
                'XLabel',obj.conditionlabels,'YLabel',strrep(param,'_',' ')) ; 
            
            if any(strcmp(param,{'markov_G0','markov_G1'}))
                
                % try to get a number of states
                if ~isempty(obj.globalmaps)
                    Ns = size(obj.globalmaps,2) ; 
                elseif isfield(obj.stats(1),'duration')
                    Ns = size(obj.stats(1).duration,2) ; 
                elseif isfield(obj.stats(1),'coverage')
                    Ns = size(obj.stats(1).coverage,2) ; 
                elseif isfield(obj.stats(1),'occurrence')
                    Ns = size(obj.stats(1).duration,2) ; 
                elseif ~isempty(obj.individual(1).label)
                    Ns = max(obj.individual(1).label) ; 
                elseif ~isempty(obj.individual(1).markov_matrix)
                    Ns = size(obj.individual(1).markov_matrix,1) ; 
                elseif ~isempty(obj.individual(1).syntax_matrix)
                    Ns = size(obj.individual(1).syntax_matrix,1) ; 
                else 
                    return
                end
                
                switch param
                    case 'markov_G0'
                        dof = (Ns-1)^2 ; 
                    case 'markov_G1'
                        dof = Ns*(Ns-1)^2 ; 
                end
                
                % find G|p=0.05
                fun = @(G,dof,target) 1-chi2cdf(G,dof)-target ; 
                p05= fzero(@(G) fun(G,dof,0.05),dof) ; 
                xl = xlim ; 
                hold on
                plot(xl,log10([p05,p05]),'k--')
                
            end
            
        case {'markov_p0','markov_p1','syntax_p_random'}
            
             % number of conditions
            numcondition = max(obj.condition) ;
            
            % concatenate data
            X = [] ; % initialize data
            G = [] ; % initialize groups
            for i = 1:numcondition
                switch param
                    case 'markov_p0'
                        P = mafdr(obj.stats(i).markov.p0,'bhfdr','true') ;
                        label = sprintf('Fraction not\nzero-order Markov') ;
                        X(i) = sum(P<0.05)/length(P) ;
                    case 'markov_p1'
                        P = mafdr(obj.stats(i).markov.p1,'bhfdr','true') ; 
                        label = sprintf('Fraction not\nfirst-order Markov') ;
                        X(i) = sum(P<0.05)/length(P) ;
                    case 'syntax_p_random'
                        P = mafdr(obj.stats(i).syntax.p_random,'bhfdr','true') ; 
                        label = sprintf('Fraction with\nnon-random syntax') ; 
                        X(i) = sum(P>=0.05)/length(P) ;
                end
                 
            end
            
            % make colours for the plot
            idx = find(strcmp(varargin,'cmap')) ; 
            if ~isempty(idx)
                cmap = varargin{idx+1} ; 
            else
                rng('default')
                cmap = lhsdesign(numcondition,3) ;
            end
            
            % plot bars
            for i = 1:numcondition
                b = bar(i,X(i),'FaceColor',cmap(i,:)) ; 
                hold on 
            end
            
            % ticks etc
            set(gca,'XTick',[1:numcondition],'XTickLabel',obj.conditionlabels)
            ylabel(label)
            xlim([0.5,numcondition+0.5])
            yl = ylim ; 
            if yl(2) == 1
                yl(2) = 1.1 ;
                ylim(yl)
            end

            
             
            
     
                
            
        case {'duration','coverage','occurrence'}
            
            % number of conditions
            numcondition = max(obj.condition) ;
            
            % number of states
            Ns = size(obj.stats(1).(param),2) ; 
            
            % concatenate data
            X = [] ; % initialize data
            G = 0 ; % initialize groups
            Xax = 0 ; skip = 0 ; % initialize Xaxis locations
            for j = 1:Ns
                for i = 1:numcondition
                    X = [X;obj.stats(i).(param)(:,j)] ; 
                    G = [G;(G(end)+1)*ones(length(obj.stats(i).(param)(:,j)),1)] ; 
                    Xax = [Xax ; Xax(end)+1 + skip] ; 
                    skip = 0 ; 
                end
                skip = 1 ; 
            end
            G(1) = [] ; Xax(1) = [] ; 
            
            % make colours for the plot
            idx = find(strcmp(varargin,'cmap')) ; 
            if ~isempty(idx)
                cmap = varargin{idx+1} ; 
            else
                rng('default')
                cmap = lhsdesign(Ns,3) ;
                al = linspace(0,0.75,numcondition) ; al = repmat(al',Ns,1) ; 
                cmap = (1-al).*repelem(cmap,numcondition,1) + al.*[1,1,1] ;
            end
            
            % plot data
            plt = microstate.external.boxplot_LT(X,'Groups',G,'GroupColors',cmap,...
                'XLabel',obj.conditionlabels,'YLabel',param,'GroupXAxis',Xax') ;   
            xlabel('Microstate')
            
            
        case {'std_duration','std_coverage','std_occurrence'}
            
            % remove std from param
            param = param(5:end) ; 
            
            % number of conditions
            numcondition = max(obj.condition) ;
            
            % number of states
            Ns = size(obj.stats(1).(param),2) ; 
            
            % concatenate data
            X = [] ; % initialize data
            G = [] ; % initialize groups
            for i = 1:numcondition
                X = [X;nanstd(obj.stats(i).(param),[],2)] ; 
                G = [G;i*ones(length(nanstd(obj.stats(i).(param),[],2)),1)] ; 
            end
            
            
            
            % make colours for the plot
            idx = find(strcmp(varargin,'cmap')) ; 
            if ~isempty(idx)
                cmap = varargin{idx+1} ; 
            else
                rng('default')
                cmap = lhsdesign(numcondition,3) ;
            end
            
            param = ['std_' param] ; 
            
            % plot data
            plt = microstate.external.boxplot_LT(X,'Groups',G,'GroupColors',cmap,...
                'XLabel',obj.conditionlabels,'YLabel',strrep(param,'_',' ')) ; 
            
        case 'autoinformation'
            
            % number of conditions
            numcondition = max(obj.condition) ;
            
            % make colours for the plot
            idx = find(strcmp(varargin,'cmap')) ; 
            if ~isempty(idx)
                cmap = varargin{idx+1} ; 
            else 
                if numcondition<=9
                    cmap = brewermap(numcondition,'Set1') ; 
                else
                    rng('default')
                    cmap = lhsdesign(numcondition,3) ; 
                end
            end
            
            
            % make time axis
            t = obj.individual(1).time ; t = t-t(1) ;
            
            % concatenate data
            hold on
            for i = 1:numcondition
                mAIF = mean(obj.stats(i).autoinformation) ; 
                vAIF = std(obj.stats(i).autoinformation) ; 
                
                ti = t(1:length(mAIF)) ; ti = ti(:) ; 
                
                fill([ti ; flipud(ti)] , [mAIF'+vAIF' ; flipud(mAIF'-vAIF')] , cmap(i,:),'facealpha',0.4,'linestyle','none')
                plot(ti,mAIF,'color',cmap(i,:))
            end
            
            
        case 'globalmaps'
            
            if isempty(varargin)  % plot as matrix
                
                imagesc(obj.globalmaps) ; 
                colorbar ; 
                xlabel('Microstate #')
                
                switch obj.individual(1).modality
                    case {'eeg','meg'}
                        colormap(gca,brewermap(64,'RdBu')) ;
                        caxis([-max(abs(obj.globalmaps(:))),max(abs(obj.globalmaps(:)))]) ; 
                        ylabel('Sensor')
                    case {'source','ampenv'}
                        colormap(gca,brewermap(64,'OrRd')) ; 
                        caxis([0 max(abs(obj.globalmaps(:)))]) ; 
                        ylabel('ROI')
                end
                
            elseif isnumeric(varargin{1}) % plot as topography
                % Select x and y coordinates and labels of the channels in the data
                x = varargin{1}(:,1);
                y = varargin{1}(:,2);

                % center at origin
                x = x-mean(x) ;
                y = y-mean(y) ; 

                % Make head mask
                Rhead = 1.05*max(sqrt(x.^2+y.^2)) ; 
                theta = linspace(0,2*pi,101) ; 
                mask(:,1) = Rhead*cos(theta) ; mask(:,2) = Rhead*sin(theta) ; 

                % Make grid for interpolation
                Xi = Rhead*linspace(-1,1,50) ; 
                [Xi,Yi] = meshgrid(Xi,Xi) ; 
                GridMask = nan(size(Xi)) ; GridMask(sqrt(Xi.^2+Yi.^2)<Rhead) = 1 ; 

                % Decide layout for plot
                m = 1 ; n = 1 ; 
                k = size(obj.globalmaps,2) ; 
                while m*n < k
                    n=n+1 ; 
                    if m*n < k
                        m=m+1 ; 
                    end
                end

                % Get colormap
                switch obj.individual(1).modality
                    case {'eeg','meg'}
                        cmap = brewermap(64,'RdBu') ;
                    case {'source','ampenv'}
                        cmap = brewermap(64,'OrRd') ; 
                end

                % Loop over maps
                for i = 1:k
                    F = scatteredInterpolant(x,y,obj.globalmaps(:,i)) ; 
                    Zi = F(Xi,Yi) ; 

                    % mask image
                    Zi = Zi.*GridMask; 
                    Zmax = max(Zi(:)) ; 

                    % Plot
                    subplot(m,n,i), hold on
                    surf(Xi,Yi,Zi,'linestyle','none','facecolor','interp')
                    view([0,0,1])
                    plot3(x,y,repmat(Zmax+1,length(x)),'ok','MarkerSize',5/max(m,n))
                    plot3(mask(:,1),mask(:,2),repmat(Zmax+1,length(mask)),'k','LineWidth',2)

                    % nose
                    noseX = [-0.2,0,0.2] ; 
                    noseY = [sqrt(1-0.2^2),0.2+sqrt(1-0.2^2),sqrt(1-0.2^2)] ; 
                    plot3(noseX*Rhead,noseY*Rhead,repmat(Zmax+1,length(noseX)),'k','LineWidth',2) ; 

                    % ears
                    R = 0.2 ; 
                    phi = acos(R/2) ; 
                    phi = linspace(phi,2*pi-phi,101) ; 
                    earX = R*cos(phi) ; 
                    xscale = 0.25*(cos(2*pi*((1:length(earX))-1)/length(earX)))+0.75 ; 
                    earX = earX.*xscale-1 ; 
                    earY = R*sin(phi) ;
                    plot3(earX*Rhead,earY*Rhead,repmat(Zmax+1,length(earX)),'k','LineWidth',2) ; 
                    plot3(-earX*Rhead,earY*Rhead,repmat(Zmax+1,length(earX)),'k','LineWidth',2) ; 
                    axis off
                    hold off
                    axis equal

                    ax = gca ; 
                    colormap(ax,cmap) ;
                    caxis([-max(abs(obj.globalmaps(:,i))),max(abs(obj.globalmaps(:,i)))])
                end

            elseif isstruct(varargin{1})
                
                % LAYOUT FILE
                layout = varargin{1} ; 
                
                % Decide layout for plot
                m = 1 ; n = 1 ; 
                k = size(obj.globalmaps,2) ; 
                while m*n < k
                    n=n+1 ; 
                    if m*n < k
                        m=m+1 ; 
                    end
                end
                
                % Get modality of layout structure
                issens = isfield(layout,'outline') && isfield(layout,'pos') ; 
                issource = isfield(layout,'tissue') ; 

                % Get colormap
                switch obj.individual(1).modality
                    case {'eeg','meg'}
                        cmap = brewermap(64,'RdBu') ;
                        if issource
                            error('The supplied layout file is a source space atlas, while the modality is M/EEG') ; 
                        end
                    case {'source','ampenv'}
                        cmap = brewermap(64,'OrRd') ; 
                        if issens
                            error('The supplied layout file is a sensor space topography, while the modality is source/ampenv') ; 
                        end
                end
                
                % Make grid for interpolation
                if issens
                    Xi = linspace(-0.5,0.5,50) ; 
                    [Xi,Yi] = meshgrid(Xi,Xi) ; 
                    GridMask = inpolygon(Xi,Yi,layout.outline{1}(:,1),layout.outline{1}(:,2)) ;
                elseif issource
                    p = load(fullfile(mspath,'+external','fieldtrip','template','anatomy','inflated')) ;
                    p.facecolor = 'interp' ; 
                    p.linestyle = 'none' ; 
                    
                    p0 = p ; 
                    p0.facecolor = 0.8*ones(1,3) ; 
                end
                     
                % Loop over maps
                for i = 1:k
                    if issens
                        F = scatteredInterpolant(layout.pos(:,1),layout.pos(:,2),obj.globalmaps(:,i)) ; 
                        Zi = F(Xi,Yi) ; 
                        Zi(Zi>max(obj.globalmaps(:,i))) = max(obj.globalmaps(:,i)) ; 
                        Zi(Zi<min(obj.globalmaps(:,i))) = min(obj.globalmaps(:,i)) ; 

                        % mask image
                        Zi(~GridMask) = nan ; 
%                         Zi = Zi.*GridMask; 
                        Zmax = max(Zi(:)) ; 

                        % Plot
                        subplot(m,n,i), hold on
                        surf(Xi,Yi,Zi,'linestyle','none','facecolor','interp')
                        view([0,0,1])
                        plot3(layout.pos(:,1),layout.pos(:,2),repmat(Zmax+1,length(layout.pos)),'ok','MarkerSize',5/max(m,n))
                        
                        % outline
                        for j = 1:length(layout.outline)
                            plot3(layout.outline{j}(:,1),layout.outline{j}(:,2),repmat(Zmax+1,length(layout.outline{j})),'k','LineWidth',2) ; 
                        end
                        
                        axis off
                        hold off
                        axis equal

                        ax = gca ; 
                        colormap(ax,cmap) ;
                        caxis([-max(abs(obj.globalmaps(:))),max(abs(obj.globalmaps(:)))])
                        
                    elseif issource
                        
                        % Plot
                        subplot(m,n,i)
                        patch(p0) ;
                        
                        % See if colour scale is specified
                        idx_colscale = find(strcmpi(varargin,'cscale')) ;
                        if isempty(idx_colscale)
                            cscale = [0,1] ; 
                        else
                            cscale = varargin{idx_colscale+1} ; 
                        end
                        
                        C = nan(length(p.vertices),1) ; 
                        for j = 1:length(layout.tissue) ; 
                            ind = layout.Vertices{j} ; 
                            C(ind) = obj.globalmaps(j,i) ; 
                        end

                        p.facevertexcdata = C; 
                        patch(p)
                        
                        axis off
                        hold off
                        axis equal

                        ax = gca ; 
                        colormap(ax,cmap) ;
                        caxis(cscale*max(abs(obj.globalmaps(:,i))))
                        
                        camlight headlight ; lighting gouraud ; material dull
                        
                    else
                        error('The supplied layout file is not a valid format') ; 
                    end
                end
                
                
            else
                error('If plotting maps, the first input must either be a Nx2 double or a mesh structure')
            end
        
            
           
            
            
        otherwise
            warning('Unable to recognise option %s',param)
    end
    
    
end
