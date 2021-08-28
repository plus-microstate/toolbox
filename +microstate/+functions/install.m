function install()
% Function to install the latest version of the toolbox. 

% Make options for webpage - timeout after 2 minutes
webopts = weboptions('Timeout',120);

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
        if contains(filename,'.mat') && ~isempty(regexp(filename,'example\d_','ONCE'))
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
        
        
        
end    
