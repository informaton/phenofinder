function [demographic_data] = loadDemographics(pathname,filename)
% [demographic_data] = loadDemographics(pathname,filename)
% helper function that loads in demographic data for phenofinder.
% Pathname is the full directory path which specifies
%where the PSD files (with extension .psa.txt) are located.  This function
%will be called from the main PSD_analyzer function/program.
% filename (optional) is a string specifying the filename to load.  If it
% is not included the default filename 'WSC_DemographSlim is used.
%
%
% demographic_data is a cell with the following elements:
%  column_names are the column headers found in the first non-# delimited
% row of the input file
%  data is a cell array of cells. The outter cell contains the same number
% of elements as files found in the pathname directory.  
%  The next inner cell array contains as many columns as found in each file
% (same number as the number of elements in column_names).
%  The innermost values at this point correspond to the rows found in each
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
if(nargin<2)
    filename='WSC_DemographSlim.txt';
end;
commentStyle = '#';

fid = fopen(fullfile(pathname,filename),'r');
    
    
% demographic_data = textscan(fid,'%s %*s %c %f','commentstyle',commentStyle,'headerlines',1);
demographic_data = textscan(fid,'%s %c %f','commentstyle',commentStyle,'headerlines',0);
fclose(fid);