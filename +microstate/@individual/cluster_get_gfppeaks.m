function [xpeak,peaks,peak2sample] = cluster_get_gfppeaks(obj,nsample)
% Extract GFP peaks of the data. 
% 
% Useage: 
% Inputs: 
% - x: M/EEG data. TxP, where T is the number of time points and P is the
%      number of channels. 
% Outputs: 
% - xpk: Maps at peaks. NxP, where N is the number of peak maps. 
% - ind: Designation of maps to time series.
% 
% Global field potential is defined as the standard deviation of the EEG
% data at each time point. Peaks in GFP correspond to regions of highest
% signal to noise ratio, so are used to extract microstate maps [1]. 
% 
% References: 
% [1] Koenig et al (1999). Eur Arch Psychiatry Clin Neurosci 249. 

    
    % Calculate GFP
    if isempty(obj.gfp)
        obj = obj.calculate_gfp ; 
    end
    if isempty(obj.gfp)
        error('GFP cannot be calculated')
    end
    
    % Check nsample input
    if nargin < 2
        nsample = 1 ; 
    end
    
    % Check it is valid
    if ~isnumeric(nsample) || length(nsample) ~= 1
        error('Input nsampl must be either 0, a float 0-1, or an integer >= 1')
    end
    
    % Find the location of the GFP peaks
    [~,peaks] = findpeaks(obj.gfp) ; % peaks of gfp
    
    
    % Get the data at these peaks
    if isempty(obj.data)
        error('Data must be supplied in the microstate.individual structure')
    end
    xpeak = obj.data(peaks,:) ; 
    
    % Get peaks2ind
    if nsample ~= 1 && nargout == 3 % only valid case is if nsample = 1
        error('Output peaks2ind only valid if nsample = 1')
    elseif nsample == 1 && nargout == 3
        dx = floor(sum([peaks(2:end) ; peaks(1:end-1)],1)/2) ; % points at which microstates switch
        dx = [0 , dx , length(obj.gfp)] ; % add first and final points
        peak2sample = zeros(length(obj.gfp),1) ; % initialise vector of microstate classes
        for i = 1:length(dx)-1 % assign each class to the vector
            peak2sample(dx(i)+1:dx(i+1)) = i ; 
        end
    end
    
    % Deal with sampling ---
    % nsample=0: Do not get peaks
    if nsample == 0
        xpeak = [] ; 
        return
        
    % nsample=1: Use all peaks
    elseif nsample == 1
        return % everything is already calculated
        
    % nsample between 0 and 1: Get 100*nsample% of gfp peaks
    elseif (nsample > 0)&&(nsample < 1) 
        numpeaks = round(nsample*size(xpeaks,1)) ; % get number of peaks
        peaks = randperm(size(xpeaks,1),numpeaks) ; % select some peaks
        xpeak = xpeak(peaks,:) ; % extract those peaks

    elseif nsample > 1 && (round(nsample) == nsample) % integer greater than 1

        % check nsample < number of peaks
        if nsample > size(xpeak,1)
            warning('Fewer GFP peaks than number of samples selected in nsample, keeping all peaks')
        end
        
        numpeaks = nsample ; % get number of peaks
        peaks = randperm(size(xpeaks,1),numpeaks) ; % select some peaks
        xpeak = xpeak(peaks,:) ; % extract those peaks


    else % error in nsample

        error('Input nsample must be either 0, a float 0-1, or an integer >= 1')

    end % finished with nsample input
    

end