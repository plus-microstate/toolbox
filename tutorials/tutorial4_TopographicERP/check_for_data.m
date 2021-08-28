% Get location of the toolbox and specify data directory
toolbox_path = microstate.functions.toolbox_path ; % location where the toolbox is installed
datadir = fullfile(fileparts(toolbox_path),'tutorials','tutorial4_TopographicERP') ; % path to folder containing the data
addpath(datadir) % add this to the path

% Check if the data is in the data directory
if ~exist(fullfile(datadir,'tutorial4_data.mat'),'file')
    answer = questdlg('Data for tutorial 4 has not been downloaded. Download now?',...
        'Download data?',...
        'Yes','No','Yes') ;

    switch answer
        case 'Yes'

            % Install the data
            datafile = fullfile(datadir,'tutorial4_data.mat') ;
            url = ['https://github.com/plus-microstate/toolbox/raw/master/' datafile] ; 
            url = replace(url,filesep,'/') ; 
            url = replace(url,'+','%2B') ; 
            webopts = weboptions('Timeout',120);
            websave(datafile,url,webopts) ; 
            clear datafile url webopts answer

        case 'No'

            % Finish running the data
            errordlg('Cannot run tutorial without downloading data','Cancel tutorial')
            clear answer
            return

    end

end