function obj = cluster_globalmaps2individual(obj,makeIndividualMaps)
% Align global maps to each individual 

if nargin<2
    makeIndividualMaps = false ; 
end
validateattributes(makeIndividualMaps,{'logical'},{'scalar'},'cluster_globalmaps2individual','makeIndividualMaps') ; 

if isempty(obj.globalmaps)
    error('Expected cohort object to have global maps. Use cluster_global or cluster_globalkoptimum to obtain global maps')
end

for i = 1:length(obj.individual)
    obj.individual(i).maps = obj.globalmaps ; % assign the global maps to 
                                          % the individual
    obj.individual(i) = obj.individual(i).cluster_alignmaps ; % get labels
    if makeIndividualMaps
        obj.individual(i) = obj.individual(i).cluster_label2maps ; % make individual maps
    end
end

obj = microstate.functions.process_append(obj,'Aligned global maps to individual') ;
   
end

