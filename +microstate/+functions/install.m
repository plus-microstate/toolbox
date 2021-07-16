function install(option)
% Function to install the latest version of the toolbox. 
% 
% Useage: 
% microstate.functions.install or microstate.functions.install('full') to
% run an update and full install. 
% 
% microstate.functions.install('-data') to not include
% example data for tutorials in download to reduce file size. 
% 
% microstate.functions.install('data') to only install example data.
% 
% microstate.functions.install(filepath) to download the most recent
% version of a particular file, where filepath is the path to the file. 

if nargin<1
     option = 'full' ; 
end
option = lower(option) ; 
webopts = weboptions('Timeout',120);

switch option
    %% Full install
    case {'full','-data'}
        % Get path to toolbox and current version
        [path,currentversion] = microstate.functions.toolbox_path ; 
        
        % Get latest version
        tmpfile = sprintf('tmp%09d.txt',randi(10^9-1)) ; 
        tmpfile = fullfile(fileparts(path),tmpfile) ; 
        websave(tmpfile,'https://raw.githubusercontent.com/plus-microstate/toolbox/master/%2Bmicrostate/%2Bfunctions/toolbox_path.m',webopts) ; 
        fid = fopen(tmpfile) ; 
        str = fscanf(fid,'%s') ; 
        fclose(fid) ; 
        delete(tmpfile) ; 
        newestversion = regexp(str,"v\d.\d",'match') ; 
        newestversion = newestversion{1} ; 
        
        if strcmpi(newestversion,currentversion)
            % do not need to update any files
            fprintf('Current version (%s) up to date...\n',newestversion)
            
        else % update the files which need updating
            
            fprintf('Updating +microstate from version %s to %s...\n',currentversion,newestversion)
            
            tmpfile = sprintf('update-file-tmp%09d-DO-NOT-DELETE.zip',randi(10^9-1)) ; 
            tmpfile = fullfile(fileparts(path),tmpfile) ; 
            websave(tmpfile,'https://github.com/plus-microstate/toolbox/archive/refs/heads/master.zip',webopts) ; 
            
            % Get all files in the zip folder using java
            zipFile = java.io.File(tmpfile) ; 
            zipFile = org.apache.tools.zip.ZipFile(zipFile) ; 
            filearray = zipFile.getEntries ; 
            
            % Make stream copier
            streamCopier = com.mathworks.mlwidgets.io.InterruptibleStreamCopier.getInterruptibleStreamCopier ; 
 
            while filearray.hasMoreElements % Loop over files
                
                file = filearray.nextElement ; 
                filename = file.getName.toCharArray' ; 
                
                % Skip directories
                if strcmp(filename(end),filesep)
                    continue
                end
                
                % Skip example data
                if ~isempty(regexp(filename,'example\d_data.mat'))
                    continue
                end
                
                fprintf('Updating file %s\n',filename) ; 
                
                % Download
                replacefilename = strsplit(filename,filesep) ; replacefilename = strjoin(replacefilename(2:end),filesep) ; 
                replacefilename = fullfile(fileparts(path),replacefilename) ; 
                unzipFile = java.io.File(replacefilename) ; 
                fileOutputStream = java.io.FileOutputStream(unzipFile) ; 
                fileInputStream = zipFile.getInputStream(file) ; 
              
                streamCopier.copyStream(fileInputStream,fileOutputStream);  
                fileOutputStream.close ; 
                
                
            end
            
            zipFile.close ;
            delete(tmpfile)
        end
        
        % Deal with .mat files
        if strcmp(option,'full')
            disp('Downloading files from GitHub-lfs including tutorial data files...')
        elseif strcmp(option,'-data')
            disp('Downloading files from GitHub-lfs excluding tutorial data files...')
        end
        [~,~,versionmatfiles] = microstate.functions.toolbox_path ; 
        
        for mat = 1:size(versionmatfiles,1)
            
            dodownload = ~versionmatfiles{mat,2} || (versionmatfiles{mat,2} && strcmp(option,'full')) ; 
            if dodownload
                
                fprintf('Downloading file %s\n',versionmatfiles{mat,3}) ; 
                
                websave(versionmatfiles{mat,1},versionmatfiles{mat,3},webopts) ; 
            end
            
        end
            
    case 'data'
        
        disp('Downloading tutorial data files from GitHub-lfs...')
        [~,~,versionmatfiles] = microstate.functions.toolbox_path ; 
        
        for mat = 1:size(versionmatfiles,1)
            
            dodownload = versionmatfiles{mat,2} ; 
            if dodownload
                websave(versionmatfiles{mat,1},versionmatfiles{mat,3},webopts) ; 
            end
            
        end
        
end    
            
        
        
