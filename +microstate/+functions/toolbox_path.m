function [path,version,versionmatfiles] = toolbox_path() ; 

path = fileparts(fileparts(which('microstate.functions.toolbox_path'))) ; 
version = 'v1.4' ; 

% Version updates
% - Major update to tutorials into order matching manuscript, including
%   addition of resting-state group MEG tutorial. This also allows for
%   download of resting and auditory stimulated MEG from the Open Science
%   Framework
% - Updated handling of downloads of tutorial data. Before all downloads of
%   the main toolbox and tutorial data were handled by the install and
%   toolbox_path functions. These functions now only handle download of the
%   main toolox. Tutorial data is handled for download within each tutorial
%   directory, meaning each tutorial dataset can be downloaded separately.
% - Added import_edf function, including the external function edfread.m 


