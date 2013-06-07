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