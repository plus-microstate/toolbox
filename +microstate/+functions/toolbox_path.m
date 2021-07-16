function [path,version] = toolbox_path() ; 

path = fileparts(fileparts(which('microstate.functions.toolbox_path'))) ; 
version = 'v1.1' ; 
