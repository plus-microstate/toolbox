function pl = boxplot_LT(x,varargin) ; 
    % Boxplot function. By default will plot a boxplot with points next to it
    % representing each data point. 
    % 
    % BASIC PLOTTING -----------
    % Inputs: 
    % - x: Data to be plotted. This can be a column vector of all points, or a
    %      NxM matrix. In the latter case, if there is no additional input
    %      'Groups', it is assumed that there are N data points in each of M
    %      groups. 
    % - Optional inputs: name value pairs. 
    %      To see all optional inputs, see the HELP section below
    % 
    % Outputs: 
    % - pl: A boxplot structure, contains all handles and all options for the
    %   plot. Can be used for updating the plot (see the UPDATING A PLOT
    %   section below)
    % 
    % Example usage: 
    %    x1 = 2+randn(15,1) ; % generate data for group 1
    %    x2 = 1.5+randn(10,1) ; % generate data for group 2
    %    x = [x1 ; x2] ; % concatenate
    %    g = [zeros(15,1) ; ones(10,1)] ; % make group labels
    %    % plot data in black and red
    %    pl = good_boxplot(x,'BoxFaceColors',[0,0,0 ; 1,0,0],'Groups',g) ; 
    %    % Notice that you cannot see the median line for the first group.
    %    % Below, we demonstrate how to update the plot. 
    % 
    % UPDATING A PLOT ---------
    % Inputs: 
    % - x: A boxplot structure, which is the output of a previous plot. 
    % 
    % Outputs: 
    % - pl: An updated boxplot structure
    % 
    % Example usage: 
    %     % Call the example usage script above to obtain two boxplots that are
    %     % red and black, and a boxplot structure pl. 
    %     pl.MedianLineWidth = 2 ; % make medians very thick
    %     pl.MedianLineColors = [1,1,1] ; % make medians white
    %     pl.update(pl) ; % update the plot

    % Empty - call help 
    if nargin == 0 || strcmp(x,'help')
        help boxplot_LT
        return
    end
    % Options help 
    if strcmp(x,'help_options')
        describe_options() ; 
        return
    end

    [x,pl,Ngroups] = check_inputs(x,varargin) ; 
    pl.update = @(a) update(a) ; 
    pl.handle.axis = gca ; 

    % We need to turn hold on, so find out if currently hold is off so we can
    % turn it back off again at the end
    currhold = ishold ; 
    if ~currhold
        cla
        hold on
    end

    % get group labels
    grouplbls = unique(pl.Groups) ; 
    
    
    %% Plot box and whisker
    
    % Loop over groups
    for i = 1:Ngroups

        idx = find(pl.Groups == grouplbls(i)) ; 
        xm = x(idx) ; 
        y = pl.Box_BoxFun{i}(xm) ; % prctile(x(j),[25;75]) ; 
        m = pl.Box_MedianFun{i}(xm) ; % median(x(j)) ;
        mx = pl.Box_WhiskerMinFun{i}(xm) ; % min(x(j)) ; 
        Mx = pl.Box_WhiskerMaxFun{i}(xm) ; % max(x(j)) ; 

        % Make a slightly rounded box for boxplot
        pl.handle.box(i) = rectangle('Position',[pl.GroupXAxis(i)-pl.Box_Width(i)/2 , y(1) , pl.Box_Width(i) , diff(y)],'curvature',pl.Box_Curvature(i,:),'EdgeColor',pl.Box_BoxLineColor(i,:),'FaceColor',pl.Box_FaceColor(i,:),'LineWidth',pl.Box_BoxLineWidth(i)) ; 
        % Add median
        pl.handle.median(i) = plot(pl.GroupXAxis(i)+[-1,1]*pl.Box_Width(i)/2,[m,m],'Color',pl.Box_MedianLineColor(i,:),'LineWidth',pl.Box_MedianLineWidth(i)) ; 
        % Add whiskers
        pl.handle.whiskersmin(i) = plot([pl.GroupXAxis(i),pl.GroupXAxis(i)],[mx,y(1)],'Color',pl.Box_WhiskerLineColor(i,:),'LineWidth',pl.Box_WhiskerLineWidth(i)) ; % lower
        pl.handle.whiskersmax(i) = plot([pl.GroupXAxis(i),pl.GroupXAxis(i)],[y(2),Mx],'Color',pl.Box_WhiskerLineColor(i,:),'LineWidth',pl.Box_WhiskerLineWidth(i)) ;  % upper
        % Add whisker caps
        pl.handle.whiskermincaps(i) = plot([pl.GroupXAxis(i),pl.GroupXAxis(i)]+pl.Box_WhiskerCapWidth(i)/2*[-1,1],[mx,mx],'Color',pl.Box_WhiskerLineColor(i,:),'LineWidth',pl.Box_WhiskerLineWidth(i)) ; 
        pl.handle.whiskermaxcaps(i) = plot([pl.GroupXAxis(i),pl.GroupXAxis(i)]+pl.Box_WhiskerCapWidth(i)/2*[-1,1],[Mx,Mx],'Color',pl.Box_WhiskerLineColor(i,:),'LineWidth',pl.Box_WhiskerLineWidth(i)) ; 
    end

    
    %% Plot points
    
    % Loop over groups
    for i = 1:Ngroups

        idx = find(pl.Groups == grouplbls(i)) ; 
        xm = x(idx) ;

        % draw plot point locations from within the grid
        y = rand(length(xm),1)-0.5 ; 
        
        % make distribution
        dst = fitdist(xm,'kernel') ; 
        yq = pdf(dst,xm) ; yq = 0.98*yq/max(yq) ; 
        y = y.*yq ; 

        % plot
        if iscell(pl.Points_MarkerFaceColor)
            facecol = pl.Points_MarkerFaceColor{i} ; 
        else
            facecol = pl.Points_MarkerFaceColor(i,:) ; 
        end
        pl.handle.points(:,i) = plot(pl.GroupXAxis(i) + y,xm,pl.Points_Marker{i},'Color',pl.Points_MarkerLineColor(i,:),'MarkerFaceColor',facecol,'MarkerSize',pl.Points_MarkerSize(i)) ; 
    end
    
    
    
    %% Decorate
    
    % set x ticks and x limits to fit the data
    set(pl.handle.axis,'XTick',pl.GroupXAxis,'LineWidth',pl.AxesLineWidth,'FontSize',pl.FontSize)
    % xlim([min(pl.BoxCentre-1.8*pl.BoxWidth/2), max(pl.BoxCentre+1.8*pl.BoxWidth/2)])
    % xlim([min(min(pl.PointsLocations,pl.BoxCentre-pl.BoxWidth/2) - 0.2*pl.BoxWidth) , max(max(pl.PointsLocations,pl.BoxCentre+pl.BoxWidth/2) + 0.2*pl.BoxWidth)])
    xlim([min(pl.GroupXAxis')-0.5 , max(pl.GroupXAxis'+0.5)])
    xtickangle(pl.XLabelAngle)
    if pl.AxisBox
        box on
    end
    
    % Background colour
    set(gca,'Color',pl.AxisColor) ; 
    if pl.AxisGrid
        grid on
    end
    set(gca,'GridColor',pl.GridColor,'GridAlpha',1)
    

    % xlabels
    if ~isempty(pl.XLabel)
        set(gca,'XTickLabels',pl.XLabel)
    end

    % ylabels
    if ~isempty(pl.YLabel)
        ylabel(pl.YLabel)
    end

    % if hold was originally off, turn it back off
    if ~currhold
        hold off
    end

    % suppress the output 
    if nargout < 1 
        clear pl
    end


    %% Check INPUTS
    function [x,pl,Ngroups] = check_inputs(x,inpl)

        % Check options are in the correct format
        if length(inpl) == 1 && iscell(inpl{1}) % in the case where options were input as a cell
            inpl = inpl{1} ; 
        end
        if mod(length(inpl),2) % check its even numbered
            error('Options must be name-value pairs')
        end
        inpl = reshape(inpl',2,length(inpl)/2)' ; % make a column vector
        for i = 1:size(inpl,1)
            if ~ischar(inpl{i,1})
                error('Options must be name-value pairs') 
            end
        end

        % Check x is not w row vector
        if size(x,1) == 1
            x = x' ; 
        end

        % Groups --- 
        if size(x,2) == 1 % case of a column vector
            pl.Data = x ; 
            pl.Groups = ones(size(x,1),1) ; % default - one group
        elseif size(x,2) > 1 % case of a matrix
            groups = repmat([1:size(x,2)],size(x,1),1) ; % assign columns as groups
            x = x(:) ; % make x a column 
            pl.Data = x ; 
            pl.Groups = groups(:) ; % reshape groups accordingly
        end
        if any(strcmp(inpl(:,1),'Groups')) % check if groups were input
            ind = find(strcmp(inpl(:,1),'Groups')) ; % index of input
            pl.Groups = inpl{ind,2} ; % set groups to input groups
            pl.Groups = pl.Groups(:) ; % ensure column, since we ensured x was a column
            inpl(ind,:) = [] ; 
        end
        Ngroups = length(unique(pl.Groups)) ; % useful for later
        pl.GroupXAxis = 1:Ngroups ; % default to each group being at a fixed location
        pl.GroupColors = default_cmap(Ngroups) ; 
        % Groups done --- 

        % Some general options --- 
        pli = struct ; 
        pli.plotLineWidth = 0.5 ; 
        pli.plotLineColor = [0,0,0] ; 
        [pli,inpl] = update_options(pli,inpl) ; % update
        pl = insert_plotpl(pl,pli,Ngroups) ; 
        ax = gca ;
        pl.XLabel = [] ;
        pl.YLabel = ax.YLabel.String ; 
        pl.AxesLineWidth = median(pl.plotLineWidth) ; 
        pl.FontSize = ax.FontSize ; 
        pl.XLabelAngle = 0 ; 
        pl.groupAxis = 'x' ; 
        pl.AxisBox = true ; 
        pl.AxisColor = [0.92,0.92,0.95] ; 
        pl.AxisGrid = true ; 
        pl.GridColor = [1,1,1]; 
        [pl,inpl] = update_options(pl,inpl) ; % update
        % General options done --- 

        % Points options --- 
        pli = struct ; 
        pli.Points_MarkerSize = 4 ; 
        pli.Points_MarkerFaceColor = pl.GroupColors ; 
        pli.Points_MarkerLineColor = pl.plotLineColor ; 
        pli.Points_Marker = 'o' ; 
        [pli,inpl] = update_options(pli,inpl) ; % update here as variables depend on this later
        pl = insert_plotpl(pl,pli,Ngroups) ;
        % Points options done --- 

        % Box and whisker options --- 
        % meaning of boxes and whiskers
        pli = struct ; 
        pli.Box_BoxFun = @(x) prctile(x, [25; 75]) ;
        pli.Box_MedianFun = @(x) nanmedian(x) ;
        pli.Box_WhiskerMinFun = @(x) max(min(x),nanmedian(x)-1.5*iqr(x)) ; 
        pli.Box_WhiskerMaxFun = @(x) min(max(x),nanmedian(x)+1.5*iqr(x)) ; 
        pli.Box_FaceColor = [1,1,1] ; 
        pli.Box_Curvature = [0.2,0.2] ; 
        pli.Box_BoxLineColor = pl.plotLineColor ;
        pli.Box_BoxLineWidth = pl.plotLineWidth ;
        pli.Box_MedianLineColor = pl.GroupColors ; 
        pli.Box_MedianLineWidth = 4*pl.plotLineWidth ; 
        pli.Box_WhiskerLineColor = pl.plotLineColor ; 
        pli.Box_WhiskerLineWidth = pl.plotLineWidth ; 
        pli.Box_WhiskerCapWidth = 0 ; % 0.3*diff(pli.Box_GridBox,[],2) ; 
        pli.Box_Width = 0.8 ; 
        [pli,inpl] = update_options(pli,inpl) ; % update here as variables depend on this later
        pl = insert_plotpl(pl,pli,Ngroups) ; % insert plot options into options structure
        % Box and whisker plot options done ---

        % Check no leftover inputs
        for m = 1:size(inpl,1)
            warning(sprintf('Unknown input: %s',inpl{m,1})) ; 
        end

        % Nested function to update the options. Given a pl structure and the
        % input pl structure, this function overwrites the default pls with the
        % inputs. 
        function [pl,inpl] = update_options(pl,inpl)
    %         if isempty(inpl)
    %             return
    %         end

            % add new options
            str = fieldnames(pl) ; 
            indrm = [] ; 
            for m = 1:size(inpl,1)
                if any(strcmp(inpl(m,1),str))
                    ind = find(strcmp(inpl(m,1),str)) ; 
                    if iscell(inpl(m,2))
                        pl.(str{ind}) = inpl{m,2} ;
                    else
                        pl.(str{ind}) = inpl(m,2) ;
                    end
                    indrm = [indrm;m] ; 
                end
            end
            inpl(indrm,:) = []; 

            str = fieldnames(pl) ; 
            for m = 1:length(str)
                if ischar(pl.(str{m})) || isa(pl.(str{m}),'function_handle')
                    tmp = pl.(str{m}) ; 
                    pl.(str{m}) = cell(1) ; 
                    pl.(str{m}){1} = tmp ; 
                end
            end

        end
        % END OF UPDATE_OPTIONS ----------

        % Nested function to generate the default colormap. Based on the 'bone'
        % colormap, but discarding the first third of the colormap as it is too
        % dark to see default black lines. 
        function cmap = default_cmap(Ngroups)
    %          cmap = colormap('bone') ;
    %          cmap = cmap(24:end,:) ; 
    %          cmap = cmap(round(linspace(1,41,Ngroups)),:) ;
            cmap = colormap(gca,'jet') ; 
            if length(cmap) > 64
                cmap = cmap(round(linspace(1,length(cmap),64)),:) ; 
            end
            cmap = cmap(8:56,:) ; 
            mincmap = 0.4 ; 
            cmap = (1-mincmap)*cmap+mincmap ; 
            cmap = cmap(round(linspace(1,size(cmap,1),Ngroups)),:) ;
        end
        % END OF DEFAULT_CMAP -----------

        % Nested function to insert fields of a structure plpl into the
        % structure pl. Importantly, all fields of plpl must have an element
        % per group, so if there are too few elements we repeat elements until
        % there is one per group. Then we place the fields into the pl
        % structure. 
        function pl = insert_plotpl(pl,plpl,Ngroups)
            % first, check all plot options are repeated for each plot
            str = fieldnames(plpl) ; 
            for m = 1:length(str)
                if size(plpl.(str{m}),1) < Ngroups
                    plpl.(str{m}) = repmat(plpl.(str{m}),Ngroups,1) ;
                    plpl.(str{m}) = plpl.(str{m})(1:Ngroups,:) ; 
                end
                pl.(str{m}) = plpl.(str{m}) ; 
            end
        end
        % END OF INSERT_PLOTpl ------------

        % Nested function to generate default locations for th points to go. 
        function [locs,inpl] = default_PointsLocations(pl,inpl,Ngroups)
            locs = pl.BoxCentre-1.4*pl.BoxWidth/2 ; % default locations

            % Now - if we have only two groups, by default plots the points to
            % the left of the first group and the right of the 2nd. This can be
            % ignored by using UseLRPoints = false as an input. 
            if Ngroups == 2
                lrlocs = true ;
                if any(strcmp(inpl(:,1),'UseLRPoints'))
                    ind = find(strcmp(inpl(:,1),'UseLRPoints')) ; 
                    lrlocs = logical(inpl(ind,2)) ; 
                    inpl(ind,:) = [] ; 
                end
                if lrlocs
                    locs(2) = pl.BoxCentre(2)+(pl.BoxCentre(2)-locs(2)) ; % pl.BoxCentre(2)+pl.BoxWidth(2)/2+(1-pl.BoxWidth(2))/4 ;
                end
            end
        end
        % END OF DEFAULT_POINTSLOCATIONS -----------

    end


    % Additional subfunction to update an already plotted structure. See
    % function description for useage. 
    function pl = update(pl)
        % remove old boxplot
        delete(pl.handle.box)
        delete(pl.handle.median)
        delete(pl.handle.whiskersmin)
        delete(pl.handle.whiskersmax)
        if isfield(pl.handle,'points')
            delete(pl.handle.points) 
        end

        % remove fields from pl structure which are not valid inputs
        pl = rmfield(pl,'handle') ; 
        pl = rmfield(pl,'update') ; 
        x = pl.Data ; % get data
        pl = rmfield(pl,'Data') ; 

        % Concatenate into input cell
        str = fieldnames(pl) ; 
        newplin = cell(1,2*length(str)) ; 
        for i = 1:length(str)
            newplin{1,2*i-1} = str{i} ; 
            newplin{1,2*i} = pl.(str{i}) ; 
        end

        % Replot updated boxplot
        pl = plot_data(x,newplin) ; 
    end

    function describe_options()
        fprintf('%s\n',...
    ["Optional Inputs (name-value pairs):" , ...
    "- 'Groups': A vector or matrix of the same size as x, denoting which" , ...
    "            group each data point belongs to. Each group is plotted as a" , ...
    "            separate boxplot. Default: If x is a column vector, all" , ...
    "            points are in a single group. If x is a matrix NxM matrix, it" , ...
    "            is assumed each column is a group." , ...
    "- 'Xlabel': Labels for each group, to go on the x-axis. Default: []" , ...
    "- 'Ylabel': The label for the y-axis. Default: []" , ...
    "- 'BoxFun','MedianFun','WhiskerMinFun','WhiskerMaxFun': Functions to set" , ...
    "            the meaning of the box and whiskers. Default: BoxFun calls" , ...
    "            the 1st and 3rd quartile (i.e. box is interquartile range)," , ...
    "            MedianFun calls the median, and WhiskerMinFun/WhiskerMaxFun" , ...
    "            calls the minimum and maximum respectively." , ...
    "- 'LineWidth': Width of the lines in the boxplot. This can be overwritten" , ...
    "            for the box, median, and whiskers independently using inputs" , ...
    "            'BoxLineWidth', 'MedianLineWidth', and 'WhiskerLineWidth'." , ...
    "            Can be a scalar (all gropus the same line width) or Ngroupsx1" , ...
    "            to adjust each group individually. Default: 0.5." , ...
    "- 'LineColors': Colours of lines in the boxplot. This can be overwritten" , ...
    "            for the box, median, and whiskers independently using inputs" , ...
    "            'BoxLineColors', 'MedianLineColors', and 'WhiskerLineColors'." , ...
    "            Can be 1x3 (all groups the same colour) or Ngroupsx3 to" , ...
    "            adjust each group individually. Default: [0,0,0] (black). " , ...
    "- 'BoxCentre': Location (on the x-axis) of the centres of the boxes." , ...
    "            Should be Ngroupsx1. Default: (1:Ngroups)'. " , ...
    "- 'BoxWidth': Width of the box. Should be a scalar (all groups equal) or" , ...
    "            Ngroupsx1 (adjust each group individually). Default: 0.55. " , ...
    "- 'BoxCurvature': Adjusts the curvature of the box. [0,0] is rectangle," , ...
    "            [1,1] is circular. Should be either 1x2 (all groups equal) or" , ...
    "            Ngroupsx2 (adjust each group individually). Default:" , ...
    "            [0.2,0.2]. " , ...
    "- 'BoxFaceColors': Colors of the faces ..." , ...

        ])
    end


    
end
