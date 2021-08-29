function plt = networks_plot(nets,layout,varargin)

if nargin<2
    layout = [] ;
end

% ensure plotting functions are on path
mspath = microstate.functions.toolbox_path ;

% Decide layout for plot
m = 1 ; n = 1 ; 
k = numel(nets) ; 
while m*n < k
    n=n+1 ; 
    if m*n < k
        m=m+1 ; 
    end
end

% Check for baseline nets
ind_baseline = find(strcmpi(varargin,'Baseline')) ; 
if ~isempty(ind_baseline)
    baseline = varargin{ind_baseline+1} ; 
    % correct for baseline
    for i = 1:length(nets)
        thisnet = nets{i} ; 
        thisnet = (thisnet-baseline)./(1-baseline) ; 
        nets{i} = thisnet ; 
    end
end

% Check for density
ind_density = find(strcmpi(varargin,'Density')) ; 
if isempty(ind_density)
    if isempty(layout)
        density = 1 ; 
    else
        density = 0.01 ;
    end
else
    density = varargin{ind_density+1} ; 
end

% Threshold at density
for i = 1:length(nets)
    thisnet = nets{i} ; 
    q = size(thisnet,1) ; 
    thisnet = thisnet+diag(nan(q,1)) ; 
    vals = sort(thisnet(:),'descend') ; 
    vals(isnan(vals)) = [] ; 
    dens_ind = floor(density*(q^2-q)); 
    thisnet(thisnet<vals(dens_ind)) = 0 ; 
    thisnet = 0.5*(thisnet+thisnet') ;
    
    
    % Reset diags to zero
    for j = 1:q
        thisnet(j,j) = 0 ; 
    end
    
    nets{i} = thisnet ; 
end


if isempty(layout) % plot as matrix
    
    for i = 1:k
        subplot(m,n,i)
        imagesc(abs(nets{i}))
        % caxis([0,1])
        axis image
    end
    
else % plot as graph
    
    % Get modality of layout structure
    issens = isfield(layout,'outline') && isfield(layout,'pos') ; 
    issource = isfield(layout,'tissue') ;
    
    % Make grid for interpolation
    if issens
        Xi = linspace(-0.5,0.5,50) ; 
        [Xi,Yi] = meshgrid(Xi,Xi) ; 
        GridMask = inpolygon(Xi,Yi,layout.outline{1}(:,1),layout.outline{1}(:,2)) ;
    elseif issource
        p = load(fullfile(mspath,'+external','fieldtrip','template','anatomy','inflated')) ;
        p.facecolor = 0.5*ones(1,3) ;
        p.facealpha = 0.2 ; 
        p.linestyle = 'none' ; 
        
        for i = 1:length(layout.tissue)
            verts = layout.Vertices{i} ; 
            verts = p.vertices(verts,:) ; 
            centr(i,:) = mean(verts) ; 
        end
    end
                     
    % Loop over maps
    for i = 1:k
        if issens

            % Plot
            subplot(m,n,i), hold on
            plot(layout.pos(:,1),layout.pos(:,2),'ok','MarkerSize',5/max(m,n))

            % outline
            for j = 1:length(layout.outline)
                plot(layout.outline{j}(:,1),layout.outline{j}(:,2),'k','LineWidth',2) ; 
            end
            
            
            axis off
            hold off
            axis equal

        elseif issource

            % Plot
            subplot(m,n,i), hold on
            patch(p) ; 
            
            % network
            G = graph(nets{i},'upper') ; 
            plot(G,'XData',centr(:,1),'YData',centr(:,2),'ZData',centr(:,3)) ; 

            

            axis off
            hold off
            axis equal

            camlight headlight ; lighting gouraud ; material dull

        else
            error('The supplied layout file is not a valid format') ; 
        end
    end
end             
                
    