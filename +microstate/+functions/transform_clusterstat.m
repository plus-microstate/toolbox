function x_transform = transform_clusterstat(x,modality)

switch modality
    case 'eeg'
        x = x-mean(x,2) ; % re-reference to average
    case 'meg'
        % do nothing
    case 'source'
        x = abs(x) ; % .^2 ; % set to magnitude
    case 'ampenv'
        % do nothing
end