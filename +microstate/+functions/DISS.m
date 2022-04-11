function DISS = DISS(x,y,modality,keep_polarity)

if nargin < 4
    keep_polarity = false ; 
end

simfun = microstate.functions.map_similarity_funhandle(modality,keep_polarity) ; 
DISS = real(sqrt(2*(1-simfun(x,y)))) ; 