function [ordered_maps,map_similarity,order] = align_maps2template(maps,templates,modality)

simfun = microstate.functions.map_similarity_funhandle(modality) ; 

for iter = 1:size(templates,2)
    R = simfun(maps',templates') ; 
    [j,i] = find(R == max(R(:))) ; 
    order(j) = i ; 
    ordered_maps(:,j) = maps(:,i) ; 
    map_similarity(j) = max(R(:)) ; 
    
    maps(:,i) = nan ; 
    templates(:,j) = nan ; 
end

idx = find(vecnorm(ordered_maps) == 0) ; 
ordered_maps(:,idx) = nan ; 

remmaps = ~any(isnan(maps)) ; 
ordered_maps = [ordered_maps , maps(:,remmaps)] ; 
order = [order,find(remmaps)] ; 