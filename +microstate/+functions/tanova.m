function [p,Dxy] = tanova(X,Y,modality,Nperm)

if nargin < 4
    Nperm = 5000 ; 
end

switch modality
    case 'eeg'
        X = X-mean(X) ; % re-reference to average
        Y = Y-mean(Y) ; 
        distfun = @(X,Y) 1-abs( (X./vecnorm(X))' * (Y./vecnorm(Y)) )' ; 
        centfun = @(v1) v1 ; 
    case 'meg'
        % do nothing
        distfun = @(X,Y) 1-abs( (X./vecnorm(X))' * (Y./vecnorm(Y)) )' ; 
        centfun = @(v1) v1 ; 
    case 'source'
        X = abs(X) ; % .^2 ; % set to magnitude
        Y = abs(Y) ; 
        distfun = @(X,Y) 1-( (X./vecnorm(X))' * (Y./vecnorm(Y)) )' ;
        centfun = @(v1) abs(v1) ; 
    case 'ampenv'
        % do nothing
        distfun = @(X,Y) 1-( (X./vecnorm(X))' * (Y./vecnorm(Y)) )' ; 
        centfun = @(v1) abs(v1) ; 
end

% Remove columns with NAN
X(:,isnan(sum(X))) = [] ; X = X' ; 
Y(:,isnan(sum(Y))) = [] ; Y = Y' ; 

if isempty(X) || isempty(Y)
    warning('No data remaining after removal of nans')
    p = nan ; Dxy = nan ; 
    return
end

% Get centroid maps
[V,D] = eig(X'*X) ; 
[~,indmax] = max(diag(D)) ; % should just be the biggest value, but put this in in case
vX = centfun(V(:,indmax)) ; % get first eigenvector

[V,D] = eig(Y'*Y) ; 
[~,indmax] = max(diag(D)) ; % should just be the biggest value, but put this in in case
vY = centfun(V(:,indmax)) ; % get first eigenvector

% Get distance
Dxy = distfun(vX,vY) ; 

% return if no permutation testing
if Nperm == 0
    p = nan ;
    return
end

% Loop over permutation
XY = [X;Y] ; lbl = [ones(size(X,1),1) ; 2*ones(size(Y,1),1)] ; 
Dxys = zeros(1,Nperm) ; Dxys(1) = Dxy ; 
for i = 1:(Nperm-1)
    
    ind = randperm(size(XY,1)) ; 
    XYs = XY(ind,:) ; 
    
    Xs = XYs(1:size(X,1),:) ; 
    Ys = XYs(size(X,1)+(1:size(Y,1)),:) ; 
    
    [~,Dxys(i+1)] = microstate.functions.tanova(Xs',Ys',modality,0) ; 
end

% Get p-value
[~,idx] = sort(Dxys,'descend') ;
idx = find(idx == 1) ; 
p = idx/Nperm ; 


end