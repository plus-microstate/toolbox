% Get location of the toolbox and specify data directory
toolbox_path = microstate.functions.toolbox_path ; % location where the toolbox is installed
datadir = fullfile(fileparts(toolbox_path),'tutorials','tutorial5_Chi2_Stimulus_Response') ; % path to folder containing the data
addpath(datadir) % add this to the path

% Check if the data is in the data directory
d = dir(fullfile(datadir,'tutorial5_data.mat')) ;
if d.bytes<150
    answer = questdlg('Data for tutorial 5 has not been downloaded. Download now?',...
        'Download data?',...
        'Yes','No','Yes') ;

    switch answer
        case 'Yes'

            % Install the data
            datafile = fullfile(datadir,'tutorial5_data.mat') ;
            url = 'https://github.com/plus-microstate/toolbox/raw/master/tutorials/tutorial5_Chi2_Stimulus_Response/tutorial5_data.mat' ; ;  
            webopts = weboptions('Timeout',120);
            websave(datafile,url,webopts) ; 
            
            % Install the layout file
            datafile = fullfile(datadir,'layout.mat') ;
            url = 'https://github.com/plus-microstate/toolbox/raw/master/tutorials/tutorial5_Chi2_Stimulus_Response/layout.mat' ; ;  
            websave(datafile,url,webopts) ; 
            clear datafile url webopts answer

        case 'No'

            % Finish running the data
            errordlg('Cannot run tutorial without downloading data','Cancel tutorial')
            clear answer
            return

    end

end
    