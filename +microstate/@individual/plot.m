function plot(obj,param,varargin)
% Visualise properties of the data and microstates
    if nargin == 1
        param = 'data' ; 
    end
    
    % Get axis
    idx = find(strcmpi(varargin,'ax')) ; 
    if isempty(idx)
        ax = gca ; 
    else
        ax = varargin{idx+1} ; 
    end
    
    % ensure plotting functions are on path
    mspath = microstate.functions.toolbox_path ; 
    addpath(fullfile(mspath,'+external','brewermap'))
    
    % Look for colour or colour map
    idx = find(strcmpi(varargin,'color')) ; 
    if ~isempty(idx)
        col = varargin{idx+1} ; 
    else
        col = []; 
    end
    idx = find(strcmpi(varargin,'colormap')) ; 
    if ~isempty(idx)
        cmap = varargin{idx+1} ; 
    else
        cmap = [] ; 
    end
    
    switch param
        case 'data'
            if isempty(col)
                col = 'k' ; 
            end
            off = 5*std(obj.data(:))*(1:size(obj.data,2)) ; 
            if size(col,1) == 1
                plot(ax,obj.time,obj.data + off,'color',col) ;
            else
                plot(ax,obj.time,obj.data + off) ;
                set(ax,'ColorOrder',col)
            end
            xlim(ax,[obj.time(1) , obj.time(end)])
            ylim(ax,[min(min(obj.data+off)) , max(max(obj.data+off))]) ;
            xlabel(ax,'time')
            switch obj.modality
                case 'eeg'
                    ylabel('EEG')
                case 'meg'
                    ylabel('MEG')
                case 'source'
                    ylabel('Source data')
                case 'ampenv'
                    ylabel('Amplitude envelope')
            end
            set(ax,'YTick',off) ; 
            set(ax,'YTickLabel',[])
            
        case 'datalabel'
            
            off = 5*std(obj.data(:))*(1:size(obj.data,2)) ;  
            xlim([obj.time(1) , obj.time(end)])
            yl = [min(min(obj.data+off)) , max(max(obj.data+off))] ; 
            ylim(yl) ;
            
            k = max(obj.label) ;
            if isempty(col)
                if k<=9
                    col = brewermap(k,'Set1') ; 
                else
                    rng('default')
                    col = lhsdesign(k,3) ; 
                end
            end
            
            dv = [1,find(diff(obj.label)),length(obj.label)] ; 
            hold on
            for i = 1:length(dv)-1
                fill([obj.time(dv(i)),obj.time(dv(i)),obj.time(dv(i+1)),obj.time(dv(i+1))],...
                    [yl(1),yl(2),yl(2),yl(1)],...
                    col(obj.label(dv(i)+1),:),...
                    'LineStyle','None',...
                    'FaceAlpha',0.4)
                drawnow
            end
            
            for i = 1:max(obj.label)
                plot(obj.time,obj.data(:,i) + off(i),'Color',col(i,:)) ;
            end
            xlabel('time')
            ylabel(obj.modality)
            set(gca,'YTick',off) ; 
            set(gca,'YTickLabel',[])

        case 'spectrum'
            
            k = max(obj.label) ; 
            if isempty(col)
                if k<=9
                    col = brewermap(k,'Set1') ; 
                else
                    rng('default')
                    col = lhsdesign(k,3) ; 
                end
            end

            [S,f] = pwelch(obj.data-mean(obj.data),[],[],[],mean(1./diff(obj.time))) ; 
            if size(col,1) == 1
                plot(f,S,'color',col) ; 
            else
                plot(f,S) ; 
                set(gca,'ColorOrder',col)
            end
            xlabel('Frequency [Hz]')
            ylabel('Power') ; 
            
        case 'gfp'
            if isempty(col)
                col = 'b' ; 
            end
            plot(obj.time,obj.gfp,'color',col)
            xlabel('time')
            ylabel('GFP')
            
        case 'label'
            if isempty(col)
                col = 'k' ; 
            end
            plot(obj.time,obj.label,'color',col)
            xlabel('time')
            ylabel('label')
            
        case 'gfplabel'
            
            k = max(obj.label) ; 
            if isempty(col)
                if k<=9
                    col = brewermap(k,'Set1') ; 
                else
                    rng('default')
                    col = lhsdesign(k,3) ; 
                end
            end
            dv = [1,find(diff(obj.label)),length(obj.label)] ; 
            basevalue = floor(0.9*min(obj.gfp)) ; 
            
            idx = find(strcmp(varargin,'tmax')) ; 
            if ~isempty(idx)
                tmax = varargin{idx+1} ; 
            else
                tmax = inf ; 
            end
            
            hold on
            for i = 1:length(dv)-1
                area(obj.time(dv(i):dv(i+1)),obj.gfp(dv(i):dv(i+1)),basevalue,'FaceColor',col(obj.label(dv(i)+1),:))
                drawnow
                if obj.time(dv(i+1))>=tmax
                    xlim([0 tmax])
                    break
                end
            end
            xlabel('time')
            ylabel('GFP')
            
        case 'syntax_graph'
            rng('default')
            k = max(obj.label) ; 
            if isempty(col)
                if k<=9
                    col = brewermap(k,'Set1') ; 
                else
                    rng('default')
                    col = lhsdesign(k,3) ; 
                end
            end
            G = digraph(obj.stats.syntax.matrix) ; 
            plot(G,'LineWidth',5*G.Edges.Weight/max(G.Edges.Weight),...
                'EdgeColor','k','MarkerSize',15,'NodeColor',col)
            
        case 'markov_graph'
            rng('default')
            k = max(obj.label) ; 
            if isempty(col)
                if k<=9
                    col = brewermap(k,'Set1') ; 
                else
                    rng('default')
                    col = lhsdesign(k,3) ; 
                end
            end
            G = digraph(obj.stats.markov.matrix) ; 
            plot(G,'LineWidth',5*G.Edges.Weight/max(G.Edges.Weight),...
                'EdgeColor','k','MarkerSize',15,'NodeColor',col)
            
        case 'syntax_matrix'
            imagesc(obj.stats.syntax.matrix,'alphadata',~eye(size(obj.stats.syntax.matrix)))
            c = colorbar ; ylabel(c,'P(transition)')
            xlabel('state(t+1)')
            ylabel('state(t)')
            if isempty(cmap)
                colormap(gca,brewermap(64,'OrRd')) ; 
            else
                colormap(gca,cmap) ; 
            end
            
        case 'markov_matrix'
            imagesc(obj.stats.markov.matrix)
            c = colorbar ; ylabel(c,'P(transition)')
            xlabel('state(t+1)')
            ylabel('state(t)')
            if isempty(cmap)
                colormap(gca,brewermap(64,'OrRd')) ; 
            else
                colormap(gca,cmap) ; 
            end
            
            
        case 'maps'
            
            if isempty(varargin)  % plot as matrix
                
                imagesc(obj.maps) ; 
                colorbar ; 
                xlabel('Microstate #')
                
                switch obj.modality
                    case {'eeg','meg'}
                         if isempty(cmap)
                            colormap(gca,brewermap(64,'RdBu')) ; 
                        else
                            colormap(gca,cmap) ; 
                        end
                        caxis([-max(abs(obj.maps(:))),max(abs(obj.maps(:)))]) ; 
                        ylabel('Sensor')
                    case {'source','ampenv'}
                         if isempty(cmap)
                            colormap(gca,brewermap(64,'OrRd')) ; 
                        else
                            colormap(gca,cmap) ; 
                        end
                        caxis([0 max(abs(obj.maps(:)))]) ; 
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
                k = size(obj.maps,2) ; 
                while m*n < k
                    n=n+1 ; 
                    if m*n < k
                        m=m+1 ; 
                    end
                end

                % Get colormap
                switch obj.modality
                    case {'eeg','meg'}
                        cmap = brewermap(64,'RdBu') ;
                    case {'source','ampenv'}
                        cmap = brewermap(64,'OrRd') ; 
                end

                % Loop over maps
                for i = 1:k
                    F = scatteredInterpolant(x,y,obj.maps(:,i)) ; 
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
                    caxis([-max(abs(obj.maps(:,1))),max(abs(obj.maps(:,i)))])
                end

            elseif isstruct(varargin{1})
                
                % LAYOUT FILE
                layout = varargin{1} ; 
                
                % Decide layout for plot
                m = 1 ; n = 1 ; 
                k = size(obj.maps,2) ; 
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
                switch obj.modality
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
                        F = scatteredInterpolant(layout.pos(:,1),layout.pos(:,2),obj.maps(:,i)) ; 
                        Zi = F(Xi,Yi) ; 
                        Zi(Zi>max(obj.maps(:,i))) = max(obj.maps(:,i)) ; 
                        Zi(Zi<min(obj.maps(:,i))) = min(obj.maps(:,i)) ; 

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
                        caxis([-max(abs(obj.maps(:))),max(abs(obj.maps(:)))])
                        
                    elseif issource
                        
                        % Plot
                        subplot(m,n,i)
                        patch(p0) ; 
                        
                        C = nan(length(p.vertices),1) ; 
                        for j = 1:max(layout.tissue) ; 
                            ind = find(layout.tissue == j) ; 
                            C(ind) = obj.maps(j,i) ; 
                        end

                        p.facevertexcdata = C; 
                        patch(p)
                        
                        axis off
                        hold off
                        axis equal

                        ax = gca ; 
                        colormap(ax,cmap) ;
                        caxis([0,max(abs(obj.maps(:)))])
                        
                        camlight headlight ; lighting gouraud ; material dull
                        
                    else
                        error('The supplied layout file is not a valid format') ; 
                    end
                end
                
                
            else
                error('If plotting maps, the first input must either be a Nx2 double or a mesh structure')
            end
             
        case 'autoinformation'
            
            aif = obj.stats.autoinformation ; 
            if isempty(col)
                col = 'k' ; 
            end
            plot(obj.time(1:length(aif)),aif,'color',col) ; 
            
        otherwise
            warning('Unable to recognise option %s',param)
    end
    
    
end