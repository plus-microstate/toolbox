function bad_samples = artifacts_json2badsamples(filename)

artfctdef_jsonstr = fileread(filename) ; 
artfctdef = jsondecode(artfctdef_jsonstr) ; 

bad_samples = [] ;
for mthd = {'clip','jump','zscore'}
    [m,artfctdef] = size_artfctdef(artfctdef,mthd{1}) ;
    for i = 1:m
        bad_samples = [bad_samples , (artfctdef.(mthd{1}).artifact(i,1)-5):(artfctdef.(mthd{1}).artifact(i,2)+5)] ; 
    end
end
bad_samples = unique(bad_samples) ;
bad_samples(bad_samples<1) = [] ; 
bad_samples(bad_samples>122880)=[] ; 

function [m,artfctdef] = size_artfctdef(artfctdef,str)
    [m,n] = size(artfctdef.(str).artifact) ; 
    if m>1 && n==1
        artfctdef.(str).artifact = artfctdef.(str).artifact' ; 
        [m,n] = size(artfctdef.(str).artifact) ; 
    end
end
    
end