function C = map_centroids(X,modality)

switch modality
    case 'eeg'
        X = X-mean(X) ; % re-reference to average
        centfun = @(v1) v1 ; 
    case 'meg'
        % do nothing
        centfun = @(v1) v1 ; 
    case 'source'
        X = abs(X) ; 
        centfun = @(v1) abs(v1) ; 
    case 'ampenv'
        % do nothing
        centfun = @(v1) abs(v1) ; 
end

X(:,isnan(sum(X))) = [] ; 
X = X' ; 

% Get centroid maps
[V,D] = eig(X'*X) ; 
[~,indmax] = max(diag(D)) ; % should just be the biggest value, but put this in in case
C = centfun(V(:,indmax)) ; % get first eigenvector