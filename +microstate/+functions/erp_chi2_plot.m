function erp_chi2_plot(info_chi2,plottype)

if nargin<2
    plottype = 'chi2vstime' ; 
end

latency = find(~isnan(info_chi2.chi2)) ; 
info_chi2.time = info_chi2.time(latency) ; 
info_chi2.chi2 = info_chi2.chi2(latency) ; 
info_chi2.uncorrected_pvals = info_chi2.uncorrected_pvals(latency) ; 
info_chi2.pearsonResiduals = info_chi2.pearsonResiduals(:,latency) ; 
hold on
switch plottype
    case 'chi2vstime'
        
        area(info_chi2.time,info_chi2.chi2) ;
        ylabel('\chi^2')
        xlabel('Time')
        
    case 'clusterperm'
        area(info_chi2.time,info_chi2.chi2) ; % plot the DISS vs time
        plot([info_chi2.time(1),info_chi2.time(end)],...
            [info_chi2.threshold, info_chi2.threshold],'--r') ; % plot the threshold for cluster permutation
        yl = ylim ; % get y-limits
        fill([info_chi2.maxcluster_time(1),info_chi2.maxcluster_time(1),info_chi2.maxcluster_time(2),info_chi2.maxcluster_time(2)],...
            [yl(1),yl(2),yl(2),yl(1)], 'k' ,'facealpha',0.2,'linestyle','none') ; % plot the significant cluster
        ylabel('\chi^2')
        xlabel('Time')
        
    case 'millisecond'

        % Plot millisecond-by-millisecond TANOVA p-values
        area(info_chi2.time,1-info_chi2.uncorrected_pvals) ; % plot 1-pval
        ylim([0.95,1])
        xlabel('Time [s]')
        
    case 'residual'
        
        % Plot pearson residuals
        if isfield(info_chi2,'maxcluster_time')
            idx = (info_chi2.time >= info_chi2.maxcluster_time(1)) & (info_chi2.time <= info_chi2.maxcluster_time(2)) ;
            info_chi2.chi2(~idx) = nan ; 
        end
            
        [~,peak] = max(info_chi2.chi2) ; 

        resid = info_chi2.pearsonResiduals(:,peak) ; 
        bar(resid)
        ylabel('Pearson Residual')
        
        % Make Bonferroni corrected threshold
        alphaBF = 0.05/length(resid) ; 
        thresh = fminsearch(@(x) abs(2*normcdf(-abs(x))-alphaBF),3) ; 
        
        hold on
        plot([0,10],thresh*[1,1],'r--')
        plot([0,10],-thresh*[1,1],'r--')
        
end
drawnow
hold off