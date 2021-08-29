function [avg_net,net_avg_all] = average_nets(nets) ; 

nets_avg = cellfun(@(C) mean(C,3),nets,'UniformOutput',false) ; 
for j = 1:size(nets,2)
    nets_j = nets_avg(:,j) ; 
    nets_j = cell2mat(permute(nets_j,[3,2,1])) ; 
    avg_net{j} = nanmean(nets_j,3) ; 
end

% get avg network
net_avg_all = nanmean(cell2mat(permute(nets(:),[3,2,1])),3) ; 
