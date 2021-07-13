function DISS = DISS(x,y,modality)

simfun = microstate.functions.map_similarity_funhandle(modality) ; 
DISS = real(sqrt(2*(1-simfun(x,y)))) ; 