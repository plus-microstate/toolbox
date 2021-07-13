% Append to process
function obj = process_append(obj,process,info)
    if nargin < 3
        info = struct ; 
    end
    process = convertCharsToStrings(process) ; 
    obj.process = [obj.process;table(process,{info},'VariableNames',{'Process','Info'})] ; 
end