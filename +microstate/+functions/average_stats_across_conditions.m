function [stats_avg,conditions] = average_stats_across_conditions(stats,conditions)

if nargin < 2
    conditions = {1:length(stats)} ; 
end
if ~iscell(conditions)
    conditions = {conditions(:)} ; 
end

stats_avg = struct ; 

for i = 1:length(conditions)
    
    if i == 1
        stats_avg = avgstruct(stats(conditions{i})) ; 
    else
        stats_avg(i) = avgstruct(stats(conditions{i})) ; 
    end
    
end

function savg = avgstruct(s) ; 

    fn = fieldnames(s) ; 
    for j = 1:length(fn)
        if isnumeric(s(1).(fn{j}))
            dim = length(size(s(1).(fn{j}))) ; 
            X = cat(dim+1,[],s.(fn{j})) ; 
            savg.(fn{j}) = nanmean(X,dim+1) ; 
        elseif isstruct(s(1).(fn{j}))
            savg.(fn{j}) = avgstruct([s.(fn{j})]) ; 
        end
    end
    
end

end