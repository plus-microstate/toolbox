function [path,version,versionmatfiles] = toolbox_path() ; 

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


path = fileparts(fileparts(which('microstate.functions.toolbox_path'))) ; 
version = 'v1.4' ; 

% List the .mat files in this version of the toolbox
% column 1: Location in toolbox
% column 2: URL for download
% column 3: Force download (e.g. if a new version of the file is included
% with this update)
versionmatfiles = {...
    fullfile('+microstate','+external','HCP230','atlas','hcp230.mat'),[],false ; ...
    fullfile('+microstate','+external','LORETA','LORETA_dipole_positions.mat'),[],false ; ...
    fullfile('+microstate','+external','LORETA','sLORETA_dipole_positions.mat'),[],false ; ...
    fullfile('+microstate','+external','fieldtrip','template','anatomy','inflated.mat'),[],false ; ...
    fullfile('+microstate','+external','fieldtrip','template','atlas','aal.mat'),[],false; ...
    fullfile('+microstate','+external','fieldtrip','template','atlas','afni_brick0.mat'),[],false ;...
    fullfile('+microstate','+external','fieldtrip','template','atlas','afni_brick1.mat'),[],false ;...
    fullfile('+microstate','+external','fieldtrip','template','atlas','brainnetome.mat'),[],false ;...
    fullfile('+microstate','+external','fieldtrip','template','layout','biosemi256.mat'),[],false ; ...
    fullfile('+microstate','+external','fieldtrip','template','layout','CTF151.mat'),[],false ;...
    fullfile('+microstate','+external','fieldtrip','template','layout','CTF275.mat'),[],false ;...
    fullfile('+microstate','+external','fieldtrip','template','layout','EEG1005.mat'),[],false ;...
    fullfile('+microstate','+external','fieldtrip','template','layout','neuromag306.mat'),[],false} ; 

for mat = 1:size(versionmatfiles,1)
    
    url = ['https://github.com/plus-microstate/toolbox/raw/master/' versionmatfiles{mat}] ; 
    url = replace(url,filesep,'/') ; 
    url = replace(url,'+','%2B') ; 
    
    versionmatfiles{mat,2} = url ; 
    versionmatfiles{mat,1} = fullfile(fileparts(path),versionmatfiles{mat,1}) ; 
end