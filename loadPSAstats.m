function [data, column_names, files, artifact_data] = loadPSAstats(pathname)
% [data, column_names, files, artifact_data] = loadPSAstats(pathname)
% helper function that loads in PSA statistics.
% Pathname is the full directory path which specifies
% where the PSD files (with extension .psa.txt) are located.  This function
% will be called from the main PSD_analyzer function/program.
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

if(nargin<1)
    pathname = pwd;
end;
files = dir(fullfile(pathname,'*.psa.txt'));
delimiter = '#';
% Example PSA file header contents

%#Power Spectral Density values from FFTs with the following parameters:
%#	CHANNEL:	O1-M2 (5)
%#	window length (seconds):	2.0
%#	FFT length (samples):	200
%#	FFT interval (taken every _ seconds):	2.0
%#	Initial Sample Rate:	100 Hz
%#	Final Sample Rate:	100 Hz
%0.0 	0.5 	1.0 	1.5 	2.0 	2.5 	3.0 	3.5 	4.0 	4.5 	5.0 	5.5 	6.0 	6.5 	7.0 	7.5 	8.0 	8.5 	9.0 	9.5	10.0	10.5	11.0	11.5	12.0	12.5	13.0	13.5	14.0	14.5	15.0	15.5	16.0	16.5	17.0	17.5	18.0	18.5	19.0	19.5	20.0	20.5	21.0	21.5	22.0	22.5	23.0	23.5	24.0	24.5	25.0	25.5	26.0	26.5	27.0	27.5	28.0	28.5	29.0	29.5	30.0	30.5	31.0	31.5	32.0	32.5	33.0	33.5	34.0	34.5	35.0	35.5	36.0	36.5	37.0	37.5	38.0	38.5	39.0	39.5	40.0	40.5	41.0	41.5	42.0	42.5	43.0	43.5	44.0	44.5	45.0	45.5	46.0	46.5	47.0	47.5	48.0	48.5	49.0	49.5	50.0	Slow	Delta	Theta	Alpha	Sigma	Beta	Gamma	Mean0_30	Sum0_30	A	O	S	E
artifact_files = dir(fullfile(pathname,'*.stats.txt'));
if(numel(artifact_files)~=numel(files))
    errordlg('Inconsistent number of artifact and psa files found!  There may be errors!','oops');
    num_files = min(numel(artifact_files),numel(files));
else
    num_files = numel(files);
end;


data = cell(num_files,1);
column_names = [];
artifact_data = cell(num_files,1);

for k=1:num_files
    files(k).name
    fid = fopen(files(k).name,'r');
    
    
    hdr_data = fgetl(fid);
    while(isempty(hdr_data)||strcmp(hdr_data(1),delimiter))
        hdr_data = fgetl(fid);
    end;
    
%     if(k==1)
        %should be at the row with column headers/labels row
        all_column_names = regexp(hdr_data,'\s+','split');
        column_names = all_column_names([1:61,end-8:end]);
        
        scanStr = [repmat('%f\t',1,61),repmat('%*f\t',1,numel(all_column_names)-69),repmat('%f\t',1,9)];
        scanStr = scanStr(1:end-2);% = 'n';
%     end;

    tic
    data{k} = cell2mat(textscan(fid,scanStr));
    fclose(fid);
    
    fid = fopen(artifact_files(k).name,'r');
    artifact_data{k} = cell2mat(textscan(fid,'%f %f %f %f %f %f %f','headerlines',2));
    fclose(fid);
    toc
end;