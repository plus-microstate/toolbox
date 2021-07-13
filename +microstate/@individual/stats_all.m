function obj = stats_all(obj)
% Calculate all available statistics 
    obj.stats = [] ; 
    obj = obj.stats_duration ; 
    obj = obj.stats_coverage ; 
    obj = obj.stats_occurrence ; 
    obj = obj.stats_hurst ; 
    obj = obj.stats_complexity ; 
    obj = obj.stats_autoinformation ; 
    obj = obj.stats_markov ; 
    obj = obj.stats_syntax ; 
    obj = obj.stats_gfp_peaksfreq ; 
    obj = obj.stats_gev ;

   
end


    