function [p,Dxy] = paired_tanova(X,Y,modality,Nperm)

% Check this is paired data
if size(X,2) ~= size(Y,2)
    error('paired_tanova can only be used for paired data')
end


if nargin < 4
    Nperm = 5000 ; 
end

switch modality
    case 'eeg'
        X = X-mean(X) ; % re-reference to average
        Y = Y-mean(Y) ; 
        dissfun = @(vX,vY) std((vX./vecnorm(vX)) - (vY./vecnorm(vY)))  ;
    case {'meg','ampenv'}
        % do nothing
        dissfun = @(vX,vY) vecnorm((vX./vecnorm(vX)) - (vY./vecnorm(vY)))  ; 
    case 'source'
        X = abs(X) ; % .^2 ; % set to magnitude
        Y = abs(Y) ; 
        dissfun = @(vX,vY) vecnorm((vX./vecnorm(vX)) - (vY./vecnorm(vY)))  ;
end


% Remove columns with NAN
idx = isnan(sum(X)) | isnan(sum(Y)) ; 
X(:,idx) = [] ; 
Y(:,idx) = [] ; 

if isempty(X) || isempty(Y)
    warning('No data remaining after removal of nans')
    p = nan ; Dxy = nan ; 
    return
end

% Get centroid maps
vX = microstate.functions.map_centroids(X,modality) ; 
vY = microstate.functions.map_centroids(Y,modality) ; 

% Get distance
Dxy = dissfun(vX,vY) ; 

% return if no permutation testing
if Nperm == 0
    p = nan ;
    return
end

% Loop over permutation
Dxys = zeros(1,Nperm) ; Dxys(1) = Dxy ;
Npart = size(X,2) ; 
for i = 1:(Nperm-1)
    
    doswitch = rand(Npart,1)<0.5 ; 
    
    Xs = X ; Xs(:,doswitch) = Y(:,doswitch) ; 
    Ys = Y ; Ys(:,doswitch) = X(:,doswitch) ; 
    
    vXs = microstate.functions.map_centroids(Xs,modality) ;
    vYs = microstate.functions.map_centroids(Ys,modality) ;
    Dxys(i+1) = dissfun(vXs,vYs) ; 
    
end

% Get p-value
p = 1-sum(Dxy>Dxys)/Nperm ; 

% % Get p-value
% [~,idx] = sort(Dxys,'descend') ;
% idx = find(idx == 1) ; 
% p = idx/Nperm ; 


end