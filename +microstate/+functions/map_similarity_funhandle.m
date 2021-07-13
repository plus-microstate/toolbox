function simfun = map_similarity_funhandle(modality,keeppolarity)

if nargin < 2
    keeppolarity = false ; 
end

if keeppolarity 
    switch modality
        case 'eeg'
            simfun = @(x,c) ( (x./vecnorm(x,2,2)) * ((c-mean(c,2))./vecnorm((c-mean(c,2)),2,2))' )' ; 
        case 'meg'
            simfun = @(x,c) ( (x./vecnorm(x,2,2)) * (c./vecnorm(c,2,2))' )'  ; 
        case 'source'
            simfun = @(x,c) ( (x./vecnorm(x,2,2)) * (c./vecnorm(c,2,2))' )' ;
        case 'ampenv'
            simfun = @(x,c) ( (x./vecnorm(x,2,2)) * (c./vecnorm(c,2,2))' )' ; 
    end
    
else

    switch modality
        case 'eeg'
            simfun = @(x,c) abs((x./vecnorm(x,2,2)) * ((c-mean(c,2))./vecnorm((c-mean(c,2)),2,2))')' ; 
        case 'meg'
            simfun = @(x,c) abs((x./vecnorm(x,2,2)) * (c./vecnorm(c,2,2))')'  ; 
        case 'source'
            simfun = @(x,c) ( (x./vecnorm(x,2,2)) * (c./vecnorm(c,2,2))' )' ;
        case 'ampenv'
            simfun = @(x,c) ( (x./vecnorm(x,2,2)) * (c./vecnorm(c,2,2))' )' ; 
    end
    
end
