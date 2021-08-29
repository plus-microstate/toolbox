% Get location of the toolbox and specify data directory
toolbox_path = microstate.functions.toolbox_path ; % location where the toolbox is installed
datadir = fullfile(fileparts(toolbox_path),'tutorials','tutorial2_RestingSource_Group') ; % path to folder containing the data
addpath(datadir) % add this to the path

% Check if the data is in the data directory
if ~exist(fullfile(datadir,'MEG-rest'),'dir')
    answer = questdlg('Data for tutorial 2 has not been downloaded. Download now?',...
        'Download data?',...
        'Yes','No','Yes') ;

    switch answer
        case 'Yes'

            % Install the data
            disp('Downloading data, this may be time consuming...')
            datafile = fullfile(datadir,'MEG-rest.zip') ; 
            url = 'https://files.de-1.osf.io/v1/resources/db9u4/providers/osfstorage/6128c83738595800dc2ebf8f/?zip=' ; 
            webopts = weboptions('Timeout',120);
            websave(datafile,url,webopts) ; 
            unzip(datafile,fullfile(datadir,'MEG-rest')) ;  
            delete(datafile)
            
            datafile = fullfile(datadir,'cluster_globalkoptimum_output.mat') ;
            url = 'https://github.com/plus-microstate/toolbox/raw/master/tutorials/tutorial2_RestingSource_Group/cluster_globalkoptimum_output.mat' ; 
            websave(datafile,url,webopts) ; 
            
            datafile = fullfile(datadir,'microstate_segmented_wpli.mat') ;
            url = 'https://github.com/plus-microstate/toolbox/raw/master/tutorials/tutorial2_RestingSource_Group/microstate_segmented_wpli.mat' ; 
            websave(datafile,url,webopts) ; 
            clear datafile url webopts answer

        case 'No'

            % Finish running the data
            errordlg('Cannot run tutorial without downloading data','Cancel tutorial')
            clear answer
            return

    end

end

datadir = fullfile(datadir,'MEG-rest') ; 
    