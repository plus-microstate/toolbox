function erp_TANOVA_plot(info_TANOVA,plottype)

if nargin<2
    plottype = 'clusterperm' ; 
end

hold on
switch plottype
    case 'clusterperm'
        area(info_TANOVA.time,info_TANOVA.DISS) ; % plot the DISS vs time
        plot([info_TANOVA.time(1),info_TANOVA.time(end)],...
            [info_TANOVA.threshold, info_TANOVA.threshold],'--r') ; % plot the threshold for cluster permutation
        yl = ylim ; % get y-limits
        fill([info_TANOVA.maxcluster_time(1),info_TANOVA.maxcluster_time(1),info_TANOVA.maxcluster_time(2),info_TANOVA.maxcluster_time(2)],...
            [yl(1),yl(2),yl(2),yl(1)], 'k' ,'facealpha',0.2,'linestyle','none') ; % plot the significant cluster
        ylabel('DISS')
        xlabel('Time')
        
    case 'millisecond'

        % Plot millisecond-by-millisecond TANOVA p-values
        area(info_TANOVA.time,1-info_TANOVA.p_sample) ; % plot 1-pval
        ylim([0.95,1])
        xlabel('Time [s]')
        
end
drawnow
hold off