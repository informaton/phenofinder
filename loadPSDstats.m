function [data, column_names, psdFilenames, artifact_data] = loadPSDstats(pathname)
% [data, column_names, files, artifact_data] = loadPSDstats(pathname)
% helper function that loads in PSD statistics.
% Pathname is the full directory path which specifies
% where the PSD files (with extension .psd.txt) are located.  This function
% will be called from the main PhenoFinder function/program.
%
% column_names are the column headers found in the first non-# delimited
% row of the input file
% data is a cell array of cells. The outter cell contains the same number
% of elements as files found in the pathname directory.  
% The next inner cell array contains as many columns as found in each file
% (same number as the number of elements in column_names).
% The innermost values at this point correspond to the rows found in each
% column and represent the time-ordered power density values.
%
% files is a structure array corresponding to the filenames of where the data in
% each cell of data (output variable) was pulled from.
%
% Hyatt Moore IV
% October 23, 2010
%> @note Example PSD header content:
%> - #Power Spectral Density values from FFTs with the following parameters: (Batch ID: 2016Jan28_06_16_04)
%> - #	CHANNEL:	C3
%> - #	window length (seconds):	5.0
%> - #	FFT length (samples):	500
%> - #	FFT interval (taken every _ seconds):	5.0
%> - #	Initial Sample Rate(Hz):	256
%> - #	Final Sample Rate(Hz):	100
%> - 0.0 	0.2 	0.4 	0.6 	0.8 	1.0 	1.2 	1.4 	1.6 	1.8 	2.0 	2.2 	2.4 	2.6	...	49.8	50.0	Slow	Delta	Theta	Alpha	Sigma	Beta	Gamma	Mean0_30	Sum0_30	A	A_type	S	E



if(nargin<1)
    pathname = pwd;
end;
[psdFilenames, psdFullFilenames] = getFilenamesi(pathname,'psd.txt');%dir(fullfile(pathname,'*.psd.txt'));
numPSDFiles = numel(psdFilenames);
delimiter = '#';

[~, artifactFullFilenames] = getFilenamesi(pathname,'stats.txt'); %dir(fullfile(pathname,'*.stats.txt'));
numArtifactFiles = numel(artifactFullFilenames);
if(numArtifactFiles>0 && numArtifactFiles~=numPSDFiles)
    errordlg('Inconsistent number of artifact and psd files found!  There may be errors!','oops');
end

data = cell(numPSDFiles,1);
column_names = [];

% We will allocate space for each PSD file, though may not actually have
% artifact data for it.  This prevents requirement of artifact data.
artifact_data = cell(numPSDFiles,1);


% We are expecting a one:one match between artifact output and non artifact
% output.
for k=1:numPSDFiles
    fprintf(1,'Parsing %s\n',psdFilenames{k});
    fid = fopen(psdFullFilenames{k},'r');
    
    hdr_data = fgetl(fid);
    while(isempty(hdr_data)||strcmp(hdr_data(1),delimiter))
        hdr_data = fgetl(fid);
    end;
    
%     if(k==1)
        %should be at the row with column headers/labels row
        all_column_names = regexp(hdr_data,'\s+','split');
        maxFrequencyToDisplay = 30;
        maxFrequencyToDisplayIndex = find(str2double(all_column_names)<=maxFrequencyToDisplay,1,'last');
        maxFrequencyIndex = find((max(str2double(all_column_names))==str2double(all_column_names)));
        numFrequencyFieldsAfterMaxDisplayFrequency = maxFrequencyIndex-maxFrequencyToDisplayIndex;  % there are ignored.
        
        column_names = all_column_names([1:maxFrequencyToDisplayIndex,maxFrequencyIndex+1:end]);
        numMetaDataFields = numel(all_column_names)-maxFrequencyIndex;
        
        
        % Create the scan string expression, which is float for all
        % frequencies of interest, skipped floats for the non interested
        % frequencies, some floats for the meta data and the artifact flag,
        % an ignored string (%*s) for the artifact_type column, and then
        % floats for the remaing 'S'tage and 'E'poch columns
        scanStr = [repmat('%f\t',1,maxFrequencyToDisplayIndex),repmat('%*f\t',1,numFrequencyFieldsAfterMaxDisplayFrequency),repmat('%f\t',1,numMetaDataFields-3),'%*s\t%f\t%f'];
%     end;

    tic
    data{k} = cell2mat(textscan(fid,scanStr));
    toc
    fclose(fid);
end

for k=1:numArtifactFiles

    fid = fopen(artifactFullFilenames{k},'r');
    artifact_data{k} = cell2mat(textscan(fid,'%f %f %f %f %f %f %f','headerlines',2));
    fclose(fid);
    toc
end;