function [path,version,versionmatfiles] = toolbox_path() ; 

% Version updates
% - Update allowing option to keep polarity when clustering in sensor
%   space, useful for ERPs/ERFs.
% - Updated methods for choosing filter order for smooth simulated random
%   walk sequence. Default remains unchanged, but you can now use a GA or
%   staircase method as well, which may be faster for long mean durations.
% - Microstate networks now run MVPA by default and give p-value and
%   confusion matrices as first and 2nd output. 
% - Spatial filter added for sensor data. 
% - Tutorials are updated to reflect these changes. 
% - Fix to resampling function to avoid edge artifacts from non-zero mean
%   time series. 
% - Filter allows option to demean, and new preprocessing function for
%   demeaning added. 


path = fileparts(fileparts(which('microstate.functions.toolbox_path'))) ; 
version = 'v1.6' ; 

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