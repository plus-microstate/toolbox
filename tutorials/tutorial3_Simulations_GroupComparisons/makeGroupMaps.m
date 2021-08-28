function GroupMaps = makeGroupMaps ; 
% This function makes group-level ground truth microstate maps in
% example3_script
    Z = zeros(5,1) ; F = @(x) sign(randn)*(0.5+0.5*rand)*normpdf([2;1;0;1;2]) ;
    GroupMaps = [F() ,  Z  ,  Z  , F() ; Z  , F() ,  Z  , Z ; Z  ,  Z  , F() , Z ; 
        Z  ,  Z  ,  Z  , F()] ; 
    GroupMaps = GroupMaps./vecnorm(GroupMaps) ; 
end
