function obj = preprocess_spatialfilter(obj,pos,kneigh)
% CARTOOL style spatial filter, see
% https://doi.org/10.3389/fneur.2019.00325 for reference

if nargin < 3
    kneigh = 6 ; 
end

if strcmp(obj.modality,'source') || strcmp(obj.modality,'ampenv')
    error('Spatial filtering should only be applied to sensor signals')
    return
end

% Are electrodes in 2d or 3d? 
dim = length(size(pos)) ; % should be either 2d or 3d

% Get distance between each pair of electrodes
for i = 1:dim
    dx(:,:,i) = pos(:,i)-pos(:,i)' ; 
end
dx = vecnorm(dx,2,3) ; 

% Check valid number of neighbours
[T,N] = size(obj.data) ; 
if kneigh > 0.5*N; 
    kneigh = floor(0.5*N) ; 
    warning('Using %d neighbours',kneigh) ; 
end
if kneigh < 4
    error('Too few electrodes for spatial filtering, at least 8 required')
end

% Loop over electrodes
X = zeros(size(obj.data)) ; 
for i = 1:N
    
    dxi = dx(i,:) ; 
    dxi(dxi==0)=1 ; 
    [dxi,ord] = sort(dxi) ; 
    dxi = dxi(1:kneigh+1) ; 
    ord = ord(1:kneigh+1) ; 
    
    Xi = obj.data(:,ord) ; 
    [Xi,idx] = sort(Xi,2) ; 
    
    dxi = repmat(dxi,T,1) ; 
    for t = 1:T
        dxi(t,:) = dxi(t,idx(t,:)) ;
    end
    
    Xi(:,1) = [] ; dxi(:,1) = [] ; 
    Xi(:,end) = [] ; dxi(:,end) = [] ; 
    
    X(:,i) = sum(Xi./dxi,2) ./ sum(1./dxi,2) ; 
end

obj.data = X ; 
options = struct ; 
options.kneigh = kneigh ; 
obj = microstate.functions.process_append(obj,'Applied spatial filter',options) ;