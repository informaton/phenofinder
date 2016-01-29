function varargout = PhenoFinder(varargin)
% PHENOFINDER M-file for PhenoFinder.fig
%      PHENOFINDER, by itself, creates a new PHENOFINDER or raises the existing
%      singleton*.
%
%      H = PHENOFINDER returns the handle to a new PHENOFINDER or the handle to
%      the existing singleton*.
%
%      PHENOFINDER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PHENOFINDER.M with the given input arguments.
%
%      PHENOFINDER('Property','Value',...) creates a new PHENOFINDER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PhenoFinder_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PhenoFinder_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PhenoFinder

% Last Modified by GUIDE v2.5 13-May-2013 10:13:10


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PhenoFinder_OpeningFcn, ...
                   'gui_OutputFcn',  @PhenoFinder_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT


% --- Executes just before PhenoFinder is made visible.
function PhenoFinder_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PhenoFinder (see VARARGIN)

global UNDO_CELL;
showBusy(handles);
cur_path = pwd;
set([handles.axes_pie,handles.text_piechart],'visible','off');
% set(handles.pb_pathname,'string',cur_path);

if(exist('psd_data.mat','file'))
    handles.user = load('psd_data.mat');
else
    [handles.user.data, handles.user.column_names, handles.user.files,...
        handles.user.artifact_data] = loadPSAstats(cur_path);
    handles.user.dataPathname = cur_path;
    data = handles.user.data;
    column_names = handles.user.column_names;
    files = handles.user.files;
    artifact_data = handles.user.artifact_data;
    dataPathname = cur_path;
    save('psd_data.mat','data','column_names', 'files','artifact_data','dataPathname');
end;

demographics_cell = loadDemographics();

if(isempty(demographics_cell))
    warndlg(sprintf('No input data found.\nPlease visit ''http://www.stanford.edu/~hyatt4'' \nto obtain a tutorial data set.'));
% else
    
    num_files = numel(handles.user.files);
    
    age_demographics = zeros(num_files,1);
    gender_demographics = char(num_files,1);
    
    
    for k = 1:numel(handles.user.files)
        index = find(strcmp(handles.user.files(k).name(1:7),demographics_cell{1}),1);
        if(isempty(index))
            disp([handles.user.files(k).name,' was not found, dummy demographic being used.']);
            age_demographics(k) = 55;
            gender_demographics(k) = 'M';
        else
            age_demographics(k) = demographics_cell{3}(index);
            gender_demographics(k) = demographics_cell{2}(index);
        end;
        
    end;
    handles.user.demographics.male_file_indices = 'M'==gender_demographics;
    handles.user.demographics.female_file_indices = 'F'==gender_demographics;
    handles.user.demographics.age_50less_file_indices = age_demographics<50;
    handles.user.demographics.age_50to65_file_indices = age_demographics>=50&age_demographics<65;
    handles.user.demographics.age_65plus_file_indices = age_demographics>=65;
    
    handles.user.demographics.age = age_demographics;
    handles.user.demographics.gender = gender_demographics;
    
    set(handles.uipanel_settings_comparisons,'title','Comparisons','foregroundcolor',[0.5 0.5 0.5]);
    
    handles.user.max_undos = 5;
    UNDO_CELL = cell(handles.user.max_undos,1);
    
    handles.user.num_ticks = 31;
    handles.user.epochs2hold = 1;
    handles.user.min_duration = true; %min/max state of epochs necessary setting
    set(handles.edit_epochs2hold,'string',handles.user.epochs2hold);
    set(handles.slider_epochs2hold,'min',1,'max',240,'sliderstep',[.01 .1],'value',1);
    
    set(handles.text_comparisons,'string','');
    
    stagesStr = {'All Night','Awake','All nonREM','Stage 1','Stage 2','Stage 3-4','Stage 5 (REM)','Difference','Relative Difference','Ratio'};
    spectrumStr = {'Frequency Spectrum','All Frequency Bands','Delta Band','Theta Band','Alpha Band','Sigma Band'};
    
    set(handles.pop_spectrum,'string',spectrumStr,'value',1);
    set(handles.pop_stage,'string',stagesStr,'value',1);
    
    set(handles.pop_ratio_top,'string',stagesStr(1:end-2),'value',1,'enable','off');
    set(handles.pop_ratio_bottom,'string',stagesStr(1:end-2),'value',2,'enable','off');
    
    handles.user.selected_handle = [];
    handles.user.settings.mean = 1;
    handles.user.settings.median = 0;
    handles.user.settings.lines_to_hold = 0;
    handles.user.num_files = numel(handles.user.data);
    
    line_width = [2;
        1.75;
        1.5;
        1.25;
        1.;
        0.75];
    
    handles.user.min_freq = 0;
    handles.user.max_freq = 30;
    handles.user.cur_freq = 0;
    handles.user.xlim = [handles.user.min_freq-0.5 handles.user.max_freq];
    
    %put the linehandles and powerlines first so that Matlab's legend function
    %can take care of this and put it in the correct order - and I can just
    %give a cell to include the first few lines as necessary..
    
    numLines = numel(line_width);
    
    handles.user.linehandles  = zeros(numLines,1);
    
    for k = 1:numLines
        handles.user.linehandles(k) = line('parent',handles.axes_main,'visible','on'...
            ,'linewidth',line_width(k),'xdata',[],'ydata',[],'displayname','Study Average - All Night',...
            'buttondownfcn',@pri_line_buttonDownFcn);
        if(k>(handles.user.settings.lines_to_hold+1))
            tmp = get(handles.user.linehandles(k),'annotation');
            tmp.LegendInformation.IconDisplayStyle = 'off';
        end;
    end;
    
    
    handles.user.bandpower_linehandle = line('parent',handles.axes_main,'visible','off'...
        ,'linewidth',line_width(2),'xdata',[],'ydata',[],'color',[1 0.5 0.2],'displayname','Band Power');
    
    tmp = get(handles.user.bandpower_linehandle,'annotation');
    tmp.LegendInformation.IconDisplayStyle = 'off';
    
    % set(handles.user.linehandles(1),'ButtonDownFcn',@pri_line_buttonDownFcn)
    
    dummy_line_for_patient_data_reference_in_legend = ...
        line('parent',handles.axes_main,'visible','off',...
        'linewidth',0.5,'linestyle',':','color',[0.5 0.5 0.5],'xdata',[],'ydata',[],...
        'ButtonDownFcn',@aux_line_buttonDownFcn,'userdata',0,'displayname',...
        'Individual patient data','xdata',[0 100],'ydata',[-100 -100]);
    
    
    %put these second so that my legend function does not have to include each
    %of these lines...
    handles.user.aux_linehandles  = zeros(handles.user.num_files,1);
    
    for k = 1:handles.user.num_files
        handles.user.aux_linehandles(k) = line('parent',handles.axes_main,'visible','off',...
            'linewidth',0.5,'linestyle',':','color',[0.5 0.5 0.5],'xdata',[],'ydata',[],...
            'ButtonDownFcn',@aux_line_buttonDownFcn,'userdata',k,'displayname',...
            handles.user.files(k).name);
        tmp = get(handles.user.aux_linehandles(k),'annotation');
        tmp.LegendInformation.IconDisplayStyle = 'off';
        
    end
    
    
    
    % hg = hggroup('parent',handles.axes_main,'displayname','Patient Data');
    % set(get(get(hg,'Annotation'),'LegendInformation'),...
    %     'IconDisplayStyle','on');
    % set(handles.user.aux_linehandles,'parent',hg);
    %
    % handles.user.hgroup = hg;
    
    % handles.user.distribution.show = false;
    
    handles.user.distribution.vertical_linehandle = line('parent',handles.axes_main,'visible','off',...
        'linewidth',2,'color',[0 0 0],'linestyle','--','displayname','Distribution Marker','xdata',[1 1]);
    
    tmp = get(handles.user.distribution.vertical_linehandle,'annotation');
    tmp.LegendInformation.IconDisplayStyle = 'off';
    
    handles.user.sem_patchhandle = patch('parent',handles.axes_main,'facecolor','g','edgecolor','k',...
        'edgealpha',0.35,'facealpha',0.35,'displayname','Standard Error of the Mean',...
        'visible','off','hittest','off');
    tmp = get(handles.user.sem_patchhandle,'annotation');
    tmp.LegendInformation.IconDisplayStyle = 'off';
    
    handles.output = handles.user.data;
    
    handles = processData(handles);
    
    set(handles.axes_main,'xlim',handles.user.xlim,'xminortick','off','fontsize',10,'nextplot','add');
    setTickMarks(handles);
    
    set(handles.axes_main_inset,'xticklabelmode','manual','xticklabel',[],...
        'xtickmode','manual','xtick',[],'xminortick','off','fontsize',10,...
        'nextplot','replace','visible','off','box','on','color',[0 1 0.7]);
    
    set(handles.cb_autoscale,'value',1);
    
    cb_autoscale_Callback(handles.cb_autoscale, eventdata, handles); %take care of the y-label,etc.
    
    % pop_layers_Callback(hObject, eventdata, handles); %take care of the number of lines to hold onto
    
    if(handles.user.settings.mean)
        set(get(handles.axes_main,'ylabel'),'string','Mean Power (\muV^2/Hz)','fontsize',13);
    elseif(handles.user.settings.median)
        set(get(handles.axes_main,'ylabel'),'string','Median Power (\muV^2/Hz)','fontsize',13);
    end;
    % set(get(handles.axes_main,'xlabel'),'string','Hz','fontsize',13);
    
    set(handles.fig,'visible','on','toolbar','figure');
    
    %just keep the initially specified number of layers/line visible
    set(handles.user.linehandles(:),'visible','off');
    
    
    set(handles.user.linehandles(1:handles.user.settings.lines_to_hold+1),'visible','on');
    
    save_for_undo = false;
    handles.user.distribution.show = false;
    
    handles = plotCurrentSettingsData(handles,save_for_undo);
    handles.user.distribution.show = true;
end

guidata(hObject,handles);

% uitoolbar;



% UIWAIT makes PhenoFinder wait for user response (see UIRESUME)
% uiwait(handles.fig);

function drawDistribution(handles)
%freq_index is the index to use for the distribution and to know where to
%draw the vertical line

if(handles.user.distribution.show)
    x = handles.user.cur_freq;
    x = min(max(x,handles.user.x(1)),handles.user.x(end));
    freq_index = (x-handles.user.x(1))*2+1;
    %     freq_index = min(max(0,(x-handles.user.x(1))*2+1),size(handles.user.power,1));  %make sure we don't get too low
    % data = handles.user.power(:,repmat(freq_index,1,256));
    try
        data = handles.user.power(handles.user.selected_file_indices,freq_index);
    catch me
        showME(me);
    end
        
%     if(handles.user.distribution.show)
%         tmp = get(handles.user.distribution.vertical_linehandle,'annotation');
%         tmp.LegendInformation.IconDisplayStyle = 'on';
%         
%         set(handles.axes_main_inset,'visible','on');
%     end;
%     
    hist(handles.axes_main_inset,data);
    hist_handle = get(handles.axes_main_inset,'children');
    num_faces = numel(get(hist_handle,'facevertexcdata'));
    set(hist_handle,'facevertexcdata',repmat([0 0 0],num_faces,1));
    inset_title_h = get(handles.axes_main_inset,'title');
    studies_count = sum(handles.user.selected_file_indices);
    set(inset_title_h,'string',['Distribution of Power for ',num2str(studies_count),' studies at ', num2str(x,'%0.1f'),'  Hz']);
    
    ylim = get(handles.axes_main,'ylim');
    
    if(handles.user.distribution.show)
        set(handles.user.distribution.vertical_linehandle,'visible','on',...
            'xdata',[x x],'ydata',ylim);
    end;
end;

function hideDistribution(handles)
%hide the inset histogram and visibility of the vertical line placed at the
%frequency of interest

tmp = get(handles.user.distribution.vertical_linehandle,'annotation');
tmp.LegendInformation.IconDisplayStyle = 'off';

delete(get(handles.axes_main_inset,'children'));
set(handles.axes_main_inset,'visible','off');
set(handles.user.distribution.vertical_linehandle,'visible','off','ydata',[0.0001, 0.0002]);

function handles = plotCurrentSettingsData(handles,save_for_undo)
if(nargin==1)
    save_for_undo = true;
end
showBusy(handles);
stage_selection = get(handles.pop_stage,'string');
stage_selectionStr = stage_selection{get(handles.pop_stage,'value')};

if(strcmp(stage_selectionStr,'Ratio')||strcmp(stage_selectionStr,'Difference')||strcmp(stage_selectionStr,'Relative Difference'))
    
    ratio_top_selection = get(handles.pop_ratio_top,'string');
    ratio_top_selectionStr = ratio_top_selection{get(handles.pop_ratio_top,'value')};
    
    handles_top = process_stage_selection(handles,ratio_top_selectionStr);
    
    ratio_bottom_selection = get(handles.pop_ratio_bottom,'string');
    ratio_bottom_selectionStr = ratio_bottom_selection{get(handles.pop_ratio_bottom,'value')};
    
    handles_bottom = process_stage_selection(handles,ratio_bottom_selectionStr);
    
    if(strcmp(stage_selectionStr,'Difference'))
        handles.user.power = handles_top.user.power-handles_bottom.user.power;
    elseif(strcmp(stage_selectionStr,'Relative Difference'))
        handles.user.power = (handles_top.user.power-handles_bottom.user.power)./(handles_top.user.power+handles_bottom.user.power);
    elseif(strcmp(stage_selectionStr,'Ratio'))
        handles.user.power = handles_top.user.power./handles_bottom.user.power;
    end;
    
    if(isempty(handles_top.user.bandPower)||isempty(handles_bottom.user.bandPower))
        handles.user.bandPower = [];
    else
        if(strcmp(stage_selectionStr,'Difference'))
            handles.user.bandPower = handles_top.user.bandPower-handles_bottom.user.bandPower;        
        elseif(strcmp(stage_selectionStr,'Relative Difference'))
            handles.user.bandPower = (handles_top.user.bandPower-handles_bottom.user.bandPower)./(handles_top.user.bandPower+handles_bottom.user.bandPower);
        elseif(strcmp(stage_selectionStr,'Ratio'))
            handles.user.bandPower = handles_top.user.bandPower./handles_bottom.user.bandPower;
        end;
    end;
    
    handles.user.sem_data = handles_top.user.sem_data+handles_bottom.user.sem_data;
    handles.user.epoch_count = min(handles_top.user.epoch_count,handles_bottom.user.epoch_count);    
    
    handles.user.cur_color = [0.42, 0.25, 0.39]; %the color for comparisons...
    
    if(handles_top.user.x==handles_bottom.user.x)
        handles.user.x = handles_top.user.x;
    else
        warndlg('x value mismatch','oops');
        handles.user.x = handles_top.user.x;
    end;
    
else
    handles = process_stage_selection(handles,stage_selectionStr);
    
end;
handles = plotData(handles,save_for_undo);
% handles = plotData(handles);


function handles = plotData(handles,varargin)
%varargin{1} is a boolean flag that determines if the a state should be
%saved for future undos or not...
showBusy(handles);
if(isempty(varargin))
    save_undo_state = true;
else
    save_undo_state = varargin{1};
end;

if(handles.user.distribution.show)
    hideDistribution(handles);
end;

% if(~isempty(handles.user.selected_handle))
%     set(handles.user.distribution.vertical_linehandle,'visible','off');
% end;
% if(isempty(handles.user.selected_handle))
%     
% %     pri_line_buttonDownFcn(handles.user.linehandles(1),[]);
% end;

epochs2minutesScalar= 2/60; %2 second epochs/ 60 seconds per minute
if(handles.user.min_duration)
    indices_with_correct_duration = (handles.user.epoch_count*epochs2minutesScalar)>=handles.user.epochs2hold;
else
    indices_with_correct_duration = (handles.user.epoch_count*epochs2minutesScalar)<=handles.user.epochs2hold;
end;


% files_indices_with_correct_demographics

selected_file_indices = indices_with_correct_duration;


if(~get(handles.cb_male,'value'))
    selected_file_indices = selected_file_indices&~handles.user.demographics.male_file_indices;
end;
if(~get(handles.cb_female,'value'))
    selected_file_indices = selected_file_indices&~handles.user.demographics.female_file_indices;
end;
if(~get(handles.cb_age_50,'value'))
    selected_file_indices = selected_file_indices&~handles.user.demographics.age_50less_file_indices;
end;
if(~get(handles.cb_age_50_65,'value'))
    selected_file_indices = selected_file_indices&~handles.user.demographics.age_50to65_file_indices;
end;
if(~get(handles.cb_age_65_plus,'value'))
    selected_file_indices = selected_file_indices&~handles.user.demographics.age_65plus_file_indices;
end;

handles.user.selected_file_indices = selected_file_indices;

if(handles.user.settings.mean)
    y_power = mean(handles.user.power(selected_file_indices,:));
elseif(handles.user.settings.median)
    y_power = median(handles.user.power(selected_file_indices,:));
end;

set(handles.user.linehandles(1),'xdata',handles.user.x,'ydata',y_power,'visible','on','color',handles.user.cur_color);

if(~isempty(handles.user.bandPower))
    if(handles.user.settings.mean)
        y_bandPower = mean(handles.user.bandPower(selected_file_indices,:));
    elseif(handles.user.settings.median)
        y_bandPower = mean(handles.user.bandPower(selected_file_indices,:));
    end;
    set(handles.user.bandpower_linehandle,'xdata',handles.user.x,'ydata',y_bandPower,'visible','on');
    
    %     x = [x,nan,x];
    %     y = [y,nan,mean(bandPower(not_nan_indices,:))];
    %
else
    set(handles.user.bandpower_linehandle,'visible','off');
end;

set(handles.user.aux_linehandles(~selected_file_indices),'visible','off');
    
for k = 1:handles.user.num_files
    if(selected_file_indices(k))
        set(handles.user.aux_linehandles(k),'xdata',handles.user.x,'ydata',handles.user.power(k,:)','visible','on');
    end;
end;

xlim = [handles.user.x(1) handles.user.x(end)];
set(handles.axes_main,'xlim',xlim);

% if(~isempty(handles.user.selected_handle))
%     ylim = get(handles.axes_main,'ylim');
%     set(handles.user.distribution.vertical_linehandle,'visible','on','ydata',ylim);
% end;

if(handles.user.distribution.show)
    drawDistribution(handles);
end

updateLegend(handles);

if(save_undo_state)
    saveUndoState(handles);
end;

showReady(handles);
guidata(handles.fig,handles);



function handles = process_stage_selection(handles,stage_selectionStr)

if(get(handles.cb_art_eeg,'value'))
    stage_indices = handles.user.stage_indices;
else
    stage_indices = handles.user.stage_indices_without_artifact;
end


switch(stage_selectionStr)
    case 'All Night'
        color = handles.user.colororder(1,:);
        for k = 1:handles.user.num_files;
            stage_indices{k} = sum(stage_indices{k}(:,2:end),2); %skip the first stage index (which is awake stage)
        end;
    case 'Awake'
        color = handles.user.colororder(2,:);
        for k = 1:handles.user.num_files;
            stage_indices{k} = stage_indices{k}(:,1);
        end;

    case 'All nonREM'
        color = handles.user.colororder(7,:);
        for k = 1:handles.user.num_files;
            stage_indices{k} = (sum(stage_indices{k}(:,2:4),2))>0;
        end;
    case 'Stage 1'
        color = handles.user.colororder(3,:);

        for k = 1:handles.user.num_files;
            stage_indices{k} = stage_indices{k}(:,2);
        end;
    case 'Stage 2'
        color = handles.user.colororder(4,:);
        for k = 1:handles.user.num_files;
            stage_indices{k} = stage_indices{k}(:,3);
        end;
    case 'Stage 3-4'
        color = handles.user.colororder(5,:);
        for k = 1:handles.user.num_files;
            stage_indices{k} = stage_indices{k}(:,4);
        end;
    case 'Stage 5 (REM)'
        color = handles.user.colororder(6,:);
        for k = 1:handles.user.num_files;
            stage_indices{k} = stage_indices{k}(:,5);

        end;        
    otherwise
        disp(stage_selectionStr)
        disp 'line 389'
        warndlg('This selection was not handled correctly in its callback function','woops');
end

spectrum_selection = get(handles.pop_spectrum,'string');
spectrum_indices = 1:handles.user.num_freq_bins;

bandPower_column = [];
switch(spectrum_selection{get(handles.pop_spectrum,'value')})
    case 'Frequency Spectrum'
        spectrum_indices = 1:handles.user.num_freq_bins;
        x = linspace(0,30,numel(spectrum_indices));
        
    case 'All Frequency Bands'
        set(handles.user.bandpower_linehandle,'displayname','Band Power {All Bands}');

        spectrum_indices = [handles.user.freq_bands.delta.indices,...
            handles.user.freq_bands.theta.indices,...
            handles.user.freq_bands.alpha.indices,...
            handles.user.freq_bands.sigma.indices];
       
%        color = [handles.user.freq_bands.delta.color;
%            handles.user.freq_bands.theta.color;
%            handles.user.freq_bands.sigma.color;
%            handles.user.freq_bands.alpha.color];
       
       x = linspace(0.5,16,numel(spectrum_indices));
       bandPower_column = [handles.user.freq_bands.delta.column;
                            handles.user.freq_bands.theta.column;
                            handles.user.freq_bands.alpha.column;
                            handles.user.freq_bands.sigma.column];
       
    case 'Delta Band'
        set(handles.user.bandpower_linehandle,'displayname','Band Power {Delta}');

        spectrum_indices = handles.user.freq_bands.delta.indices;     
%         color = handles.user.freq_bands.delta.color;
        x = linspace(0.5,4,numel(spectrum_indices));
        bandPower_column = handles.user.freq_bands.delta.column;
        
    case 'Theta Band'
        set(handles.user.bandpower_linehandle,'displayname','Band Power {Theta}');
        spectrum_indices = handles.user.freq_bands.theta.indices;
%         color = handles.user.freq_bands.theta.color;
        x = linspace(4.5,8,numel(spectrum_indices));
        bandPower_column = handles.user.freq_bands.theta.column;

    case 'Alpha Band'
        set(handles.user.bandpower_linehandle,'displayname','Band Power {Alpha}');
        spectrum_indices = handles.user.freq_bands.alpha.indices;
%         color = handles.user.freq_bands.alpha.color;
        x = linspace(8.5,12,numel(spectrum_indices));
        bandPower_column = handles.user.freq_bands.alpha.column;
        
    case 'Sigma Band'
        set(handles.user.bandpower_linehandle,'displayname','Band Power {Sigma}');
        spectrum_indices = handles.user.freq_bands.sigma.indices;
%         color = handles.user.freq_bands.sigma.color;
        x = linspace(12.5,16,numel(spectrum_indices));
        bandPower_column = handles.user.freq_bands.sigma.column;
    otherwise
        warndlg('This case was left unspecified - defaulting to all frequencies - a 1,000 apologies ...','woops');
        x = 0:30;
end;

power = zeros(handles.user.num_files,numel(spectrum_indices));
sem_data = zeros(handles.user.num_files,numel(spectrum_indices));

data = handles.user.data;

if(~isempty(bandPower_column))
    bandPower = zeros(handles.user.num_files,numel(spectrum_indices));
else
    bandPower = [];
end

if(~isempty(bandPower))
    numRows = size(bandPower_column,1);
    num_indices_per_band = 8;
end;

epoch_count = zeros(numel(stage_indices),1);

for k = 1:handles.user.num_files
    epoch_count(k) = sum(stage_indices{k}(:));
    if(epoch_count(k)>0)
        
        n = sum(stage_indices{k});
        sem_data(k,:) = std(data{k}(stage_indices{k}>0,spectrum_indices))/sqrt(n);
        
        if(handles.user.settings.mean)
            power(k,:)=mean(data{k}(stage_indices{k}>0,spectrum_indices));
        elseif(handles.user.settings.median)
            power(k,:)=median(data{k}(stage_indices{k}>0,spectrum_indices));
        end;
        
        if(~isempty(bandPower))
            for row_index = 1:numRows
                start = 1+(row_index-1)*num_indices_per_band;
                stop = row_index*num_indices_per_band;
                if(handles.user.settings.mean)
                    bandPower(k,start:stop) = mean(handles.user.data{k}(stage_indices{k}>0,bandPower_column(row_index))); %should be a straight line of values.
                elseif(handles.user.settings.median)
                    bandPower(k,start:stop) = median(handles.user.data{k}(stage_indices{k}>0,bandPower_column(row_index))); %should be a straight line of values.                 
                end;
            end;
        end;
    end;
end

handles.user.epoch_count = epoch_count;
handles.user.power = power;
handles.user.bandPower = bandPower;
handles.user.sem_data = sem_data;
handles.user.x = x;
handles.user.cur_color = color;

function shiftLines(handles)
%call this each time you update plot data... shifts the data back between
%the different layers for transparency and holding...
for k=numel(handles.user.linehandles):-1:2
    set(handles.user.linehandles(k),'xdata',get(handles.user.linehandles(k-1),'xdata'),...
        'ydata',get(handles.user.linehandles(k-1),'ydata'),...
        'color',get(handles.user.linehandles(k-1),'color'),...
        'displayname',get(handles.user.linehandles(k-1),'displayname'));
end;

function saveUndoState(handles)
global UNDO_CELL;
for k=1:numel(UNDO_CELL)-1
   UNDO_CELL{k} = UNDO_CELL{k+1}; 
end
UNDO_CELL{end} = handles;


function handles = undo(handles)
global UNDO_CELL;
if(isequal(UNDO_CELL{end},handles))
    tmpHandles = UNDO_CELL{end-1};
else
    tmpHandles = UNDO_CELL{end};
end;

if(~isempty(tmpHandles))
    handles = tmpHandles;
    for k = numel(UNDO_CELL):-1:2
        UNDO_CELL{k}=UNDO_CELL{k-1};
    end;
    UNDO_CELL{1} = [];
end;

function handles = processData(handles)
%this function will calculate the means needed for plotting from data
%parsed from SEV output psd files.
%header information was stripped, it used to be:
%index 1==0.0, 0.5,	1.0, ... index 61==30.0 then 62: Delta	63: Theta	64: Alpha	65: Sigma	66: SigmaAlpha	67: A 68:S 69:E
% data = handles.user.data;

% bins_per_epoch = 15;
num_freq_bins = 61; 

% handles.user.colormap = colormap('hot');
% handles.user.colormap(end,:)=[1 .8 0.5]; %make sure that our last result is not white, which is too hard to see...
% colormap(handles.user.colormap);

handles.user.num_stages_to_look_at = 5; %0,1,2,3-4,and 5

handles.user.colororder = zeros(handles.user.num_stages_to_look_at+2,3);
handles.user.colororder(1,:) = [0 0 1]; %all night
handles.user.colororder(2,:) = [0 1 1]; %awake
handles.user.colororder(3,:) = [0, 0.5, 0];
handles.user.colororder(4,:) = [1 0 0];
handles.user.colororder(5,:) = [1.0000    0.7917         0];
handles.user.colororder(6,:) = [0, 0.75, 0.75];
handles.user.colororder(7,:) = [0.7396    0.3979    0.1250]; %non-rem

handles.user.colormap = [handles.user.colororder(3:6,:)];

colormap(handles.axes_pie,handles.user.colormap);

% handles.user.colororder(3:handles.user.num_stages_to_look_at+1,:) = handles.user.colormap(linspace(1,size(handles.user.colormap,1),handles.user.num_stages_to_look_at-1),:); %-1 b/c do not want to include the awake stage
% mean(handles.user.colororder(3:handles.user.num_stages_to_look_at+1,:))
% handles.user.colororder(end,:)=mean(handles.user.colororder(3:handles.user.num_stages_to_look_at+1,:)); %non-rem


num_files = handles.user.num_files;
x = handles.user.data;

stage_column = strmatch('S',handles.user.column_names,'exact');
if(isempty(stage_column))
    warndlg('Could not find the stage column name (S) - proceeding to crash now...','oops!');
end;
artifact_column = strmatch('A',handles.user.column_names,'exact');
if(isempty(stage_column))
    warndlg('Could not find the artifact column name (A) - proceeding to crash now...','oops!');
end;


handles.user.freq_bands.dc.column = strmatch('0.0',handles.user.column_names,'exact');

handles.user.freq_bands.delta.indices = strmatch('0.5',handles.user.column_names,'exact'):strmatch('4.0',handles.user.column_names,'exact');
handles.user.freq_bands.theta.indices = strmatch('4.5',handles.user.column_names,'exact'):strmatch('8.0',handles.user.column_names,'exact');
handles.user.freq_bands.alpha.indices = strmatch('8.5',handles.user.column_names,'exact'):strmatch('12.0',handles.user.column_names,'exact');
handles.user.freq_bands.sigma.indices = strmatch('12.5',handles.user.column_names,'exact'):strmatch('16.0',handles.user.column_names,'exact');

handles.user.freq_bands.delta.column = strmatch('Delta',handles.user.column_names,'exact');
handles.user.freq_bands.theta.column = strmatch('Theta',handles.user.column_names,'exact');
handles.user.freq_bands.alpha.column = strmatch('Alpha',handles.user.column_names,'exact');
handles.user.freq_bands.sigma.column = strmatch('Sigma',handles.user.column_names,'exact');

% handles.user.freq_bands.delta.color = [255 165 0]; %orange (ROYGBIV)
% handles.user.freq_bands.theta.color = [255 215 0]; %gold
% handles.user.freq_bands.alpha.color = [127 255 0]; %green (chartreuse
% handles.user.freq_bands.sigma.color = [135 206 250];       %blue (light sky)
% handles.user.freq_bands.rest.color = [238 130 238]; %violet

stage_indices = cell(1,num_files);
stage_indices_without_artifact = cell(1,num_files);


total_stage_count = zeros(1,handles.user.num_stages_to_look_at-1); %-1 because I do not want stage 0 to be included right now (awake stage)
total_stage_count_without_artifact = zeros(1,handles.user.num_stages_to_look_at-1);

for k = 1:num_files;
    stage_indices{k} = zeros(size(x{k},1),handles.user.num_stages_to_look_at);

    stage_indices{k}(:,1) = x{k}(:,stage_column)==0;
    stage_indices{k}(:,2) = x{k}(:,stage_column)==1;
    stage_indices{k}(:,3) = x{k}(:,stage_column)==2;
    stage_indices{k}(:,4) = (x{k}(:,stage_column)==3)|((x{k}(:,stage_column)==4));
    stage_indices{k}(:,5) = x{k}(:,stage_column)==5;

    total_stage_count = total_stage_count + sum(stage_indices{k}(:,2:end));
%     for s = 1:handles.user.num_stages_to_look_at
%         stage_indices{k}(:,s) = x{k}(:,stage_column)==s;        
%     end;
    
    stage_indices_without_artifact{k} = stage_indices{k}(~x{k}(:,artifact_column),:);
    total_stage_count_without_artifact = total_stage_count_without_artifact + sum(stage_indices_without_artifact{k}(:,2:end),1);

    %     power.allnight(k,:)=mean(x{k}(:,1:num_freq_bins));
    

end;

handles.user.stage_count = total_stage_count;
handles.user.stage_count_without_artifact = total_stage_count_without_artifact;
handles.user.stage_indices = stage_indices;
handles.user.stage_indices_without_artifact = stage_indices_without_artifact;
handles.user.num_freq_bins = num_freq_bins;



% --- Outputs from this function are returned to the command line.
function varargout = PhenoFinder_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.user.data;

% delete(hObject);

% handles.output;


% --- Executes on selection change in pop_stage.
function pop_stage_Callback(hObject, eventdata, handles)
% hObject    handle to pop_stage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_stage contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_stage
showBusy(handles);
shiftLines(handles); %allow for saving/holding of previous plots...

stage_selection = get(handles.pop_stage,'string');
stage_selectionStr = stage_selection{get(handles.pop_stage,'value')};
if(strcmp(stage_selectionStr,'Ratio')||strcmp(stage_selectionStr,'Difference')||strcmp(stage_selectionStr,'Relative Difference'))
    if(strcmp(stage_selectionStr,'Difference'))
        set(handles.uipanel_settings_comparisons,'title','Difference','foregroundcolor',[0 0 0]);
        set(handles.text_comparisons,'string','MINUS');
    elseif(strcmp(stage_selectionStr,'Relative Difference'))
        set(handles.uipanel_settings_comparisons,'title','Difference','foregroundcolor',[0 0 0]);
        set(handles.text_comparisons,'string','MINUS');
    elseif(strcmp(stage_selectionStr,'Ratio'))
        set(handles.uipanel_settings_comparisons,'title','Ratio','foregroundcolor',[0 0 0]);
        set(handles.text_comparisons,'string','OVER');
    end;
    set(handles.pop_ratio_top,'enable','on');
    set(handles.pop_ratio_bottom,'enable','on');
else
    set(handles.pop_ratio_top,'enable','off');
    set(handles.pop_ratio_bottom,'enable','off');    
    set(handles.text_comparisons,'string','');
    set(handles.uipanel_settings_comparisons,'title','Comparisons','foregroundcolor',[0.5 0.5 0.5]);
end;
set(handles.user.linehandles(1),'displayname',['Study Average - ',stage_selectionStr]);

plotCurrentSettingsData(handles);


% --- Executes during object creation, after setting all properties.
function pop_stage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_stage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function setTickMarks(handles)
%used to obtain a good tick mark resolution depending on how many freqs are
%being looked at...

column_names = cell(handles.user.num_ticks,1);
cur_x_lim = get(handles.axes_main,'xlim');
large_x_axes = (diff(cur_x_lim)-diff(handles.user.xlim))>=-5;
for k = 0:numel(column_names)-1
    if(mod(k,2)&&large_x_axes)
        column_names{k+1} = '';
    else        
        column_names{k+1} = num2str(k,'%0.1f');
    end;
end
set(handles.axes_main,'xticklabel',...
    column_names,'xtick',0:numel(column_names)-1);


% --- Executes on selection change in pop_spectrum.
function pop_spectrum_Callback(hObject, eventdata, handles)
 % event  data  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_spectrum contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_spectrum
showBusy(handles);
tmp = get(handles.user.bandpower_linehandle,'annotation');

contents = cellstr(get(hObject,'String'));
if(strcmp(contents{get(hObject,'Value')},'Frequency Spectrum'))
    tmp.LegendInformation.IconDisplayStyle = 'off';
else
    tmp.LegendInformation.IconDisplayStyle = 'on';
end;

handles = clearSelection(handles);

plotCurrentSettingsData(handles);
showBusy(handles);
setTickMarks(handles);
cb_autoscale_Callback(handles.cb_autoscale, eventdata, handles); %take care of the y-label,etc.
showReady(handles);

% --- Executes during object creation, after setting all properties.
function pop_spectrum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_spectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_art_eeg.
function cb_art_eeg_Callback(hObject, eventdata, handles)
% hObject    handle to cb_art_eeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_art_eeg
showBusy(handles);
plotCurrentSettingsData(handles); %need to select which indices to use - so call this function which process the stage selection, which process the indices to use...

% if(get(hObject,'value'))
%    disp 'Checked' 
% else
%     disp 'Unchecked'
% end



function pb_pathname_Callback(hObject, eventdata, handles)
% hObject    handle to pb_pathname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pb_pathname as text
%        str2double(get(hObject,'String')) returns contents of pb_pathname as a double


% --- Executes during object creation, after setting all properties.
function pb_pathname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pb_pathname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pb_pathname.
function pb_pathname_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pb_pathname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pathname = uigetdir('','Select Directory with Data Files');
if(pathname && ~isequal(get(handles.pb_pathname,'string'),pathname))
    set(handles.pb_pathname,'string',pathname);
    [handles.user.data, handles.user.column_names, handles.user.files] = loadPSAstats(pathname);
    guidata(hObject,handles);
end;


% --- Executes on button press in cb_autoscale.
function cb_autoscale_Callback(hObject, ~, handles)
% hObject    handle to cb_autoscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_autoscale
if(get(handles.cb_autoscale,'value'))
    set(handles.axes_main,'yticklabelmode','auto','ytickmode','auto','ylimmode','auto');    
else
    set(handles.axes_main,'yticklabelmode','manual','ytickmode','manual','ylimmode','manual');
end;

function aux_line_buttonDownFcn(hObject,eventdata)
%what to do when the user clicks on a line...

handles = guidata(hObject);
showBusy(handles);
if(~strcmp(get(handles.fig,'selectiontype'),'alt')) %'alt' refers to a ctrl-left mouse click or a right mouse click - which will activate the uicontextmenu

    if(~isempty(handles.user.selected_handle))
        set(handles.user.selected_handle,'linewidth',0.5);
        tmp = get(handles.user.selected_handle,'annotation');
        tmp.LegendInformation.IconDisplayStyle = 'off';
    else
        tmp = get(handles.user.distribution.vertical_linehandle,'annotation');
        tmp.LegendInformation.IconDisplayStyle = 'on';
        tmp = get(handles.user.sem_patchhandle,'annotation');
        tmp.LegendInformation.IconDisplayStyle = 'on';
    end;
    
    handles.user.selected_handle = hObject;      
    tmp = get(handles.user.selected_handle,'annotation');
    tmp.LegendInformation.IconDisplayStyle = 'on';

    set(handles.user.aux_linehandles,'color',[0.7 0.7 0.7]); %make it much lighter
    
    
    k = get(hObject,'userdata');
    file = handles.user.files(k);
    set(handles.tb_meta_filename,'string',['Name: ',file.name],'visible','on');
    set(handles.tb_meta_filedate,'string',['Modified: ',file.date],'visible','on');

    gender = handles.user.demographics.gender(k);
    age = num2str(handles.user.demographics.age(k)+randn(1),'%0.1f'); %randomize age with gaussian distribution
    if(gender=='M')
        gender='Male - ';
    else
        gender = 'Female - ';
    end;
    set(handles.user.selected_handle,'linewidth',1.5,'color','k','displayname',['Study ',file.name(1:7),' (',gender,age,' years)']);


    chan_name = {'C3-EEG'};
    chan_number = 3;
    win_len = 2;
    fft_len = 200;
    fft_int = 1;
    
%     [chan_name, chan_number] = textread(file.name,'#\tCHANNEL:%s\t(%d)',1,'headerlines',1);
%     win_len = textread(file.name,'#\twindow length (seconds):\t%d',1,'headerlines',2);
%     fft_len = textread(file.name,'#\tFFT length (samples):\t%d',1,'headerlines',3);
%     fft_int = textread(file.name,'#\tFFT interval (taken every _ seconds):\t%d',1,'headerlines',4);
    
    x = get(hObject,'xdata');
    y = get(hObject,'ydata');
    half_sem_data = handles.user.sem_data(k)/2;
    
    
    autoscale = get(handles.cb_autoscale,'value');
    if(~autoscale)
        ylim = get(handles.axes_main,'ylim');
    end;
    
    ydata = [y(:)+half_sem_data(:); flipud(y(:)-half_sem_data(:))];
    logscale= get(handles.cb_y_log_scale,'value');
    if(logscale)
        ydata = max(ydata(:),[y(:);y(:)]);
    end;
    set(handles.user.sem_patchhandle,'visible','on','xdata',[x(:); flipud(x(:))],'ydata',ydata);

    if(~autoscale)
        set(handles.axes_main,'ylim',ylim);  %setting the ylim overwrites the autosetting
    end;

    currentpoint = get(handles.axes_main,'currentpoint');
    %     freq_index = round(currentpoint(1)*2)+1; %*2 because we have 0.5 hz resolution; +1 because matlab is 1 based and not zero based for their indexing
    x = round(currentpoint(1)*2)/2;
    handles.user.cur_freq = x;
    
    drawDistribution(handles);

    table_data = handles.user.artifact_data{k}([1:6,end],2:end);
    [r c] = size(table_data);
    duration = table_data(:,1);
    delta = datenum([0,0,0,0,0,1]);

    table_data = num2cell(table_data);    
    tmp = datestr(duration*delta,'HH:MM:SS');
    tmp_size = size(tmp);
    tmp = [repmat(' ',tmp_size),tmp];
    table_data(:,1) = mat2cell(tmp,ones(r,1),size(tmp,2));
%     num_epochs = duration/30; %30 second epochs...
    set(handles.uitable_meta,'data',table_data); 

    pos = get(handles.uitable_meta,'position');
    extent = get(handles.uitable_meta,'extent');
    set(handles.uitable_meta,'position',[pos(1:2), extent(3:4)]);
    
    set(handles.tb_channel_name,'string',['Name: ',chan_name{:}],'visible','on');
    set(handles.tb_channel_num,'string',['Number: ',num2str(chan_number)],'visible','on');
    set(handles.tb_channel_winlen,'string',['Window Length (sec): ',num2str(win_len)],'visible','on');
    set(handles.tb_channel_fftlen,'string',['FFT Length (samples): ',num2str(fft_len)],'visible','on');
    set(handles.tb_channel_fftint,'string',['FFT Interval (sec): ',num2str(fft_int)],'visible','on');
    
    if(get(handles.cb_art_eeg,'value'))
        stage_indices = handles.user.stage_indices{k}(:,2:end); %skip the awake stage...
    else
        stage_indices = handles.user.stage_indices_without_artifact{k}(:,2:end);
    end
    
    %clean-up any previous plots
    child = get(handles.axes_pie,'children');
    if(~isempty(child))
        delete(child);
    end;
    explode_vec = zeros(1,size(stage_indices,2));
  
    stage_selection = get(handles.pop_stage,'string');
  
    switch(stage_selection{get(handles.pop_stage,'value')})
        case 'All Night'
        case 'Awake'
        case 'Ratio'
        case 'Difference'
        case 'Relative Difference'
        case 'All nonREM'
            explode_vec(1:3)=1;
        case 'Stage 1'
            explode_vec(1)=1;
        case 'Stage 2'
            explode_vec(2)=2;
        case 'Stage 3-4'
            explode_vec(3)=3;
        case 'Stage 5 (REM)'
            explode_vec(4)=4;
        otherwise
            disp(stage_selection{get(handles.pop_stage,'value')});
            disp 'line 1090'
            warndlg('This selection was not handled correctly in its callback function','woops');
    end
    labelStr = {'Stage I','Stage II','Stage III-IV','Stage V (REM)'};
    pie(handles.axes_pie,sum(stage_indices,1)+.1,explode_vec);
    legend(handles.axes_pie,'show',labelStr,'Location','SouthEastOutside');
    set(handles.text_piechart,'visible','on');
    %     set(handles.axes_pie,'visible','on');

    saveUndoState(handles);
    updateLegend(handles);
    guidata(hObject,handles);
    
%     metaData = regext('# CHANNEL:\s+(?<channel>[^\n\r]+)[^#]*#window legnth (seconds):\s+(?<winlen>[0-9]*\.?[0-9]*)');
%     #Power Spectral Density values from FFTs with the following parameters:
% #	CHANNEL:	C3-M2 (3)
% #	window length (seconds):	2.0
% #	FFT length (samples):	200
% #	FFT interval (taken every _ seconds):	2.0

end
showReady(handles);

function handles = clearSelection(handles)

if(~isempty(handles.user.selected_handle))
    set(handles.user.selected_handle,'linewidth',0.5);
    tmp = get(handles.user.selected_handle,'annotation');
    tmp.LegendInformation.IconDisplayStyle = 'off';
    tmp = get(handles.user.sem_patchhandle,'annotation');
    tmp.LegendInformation.IconDisplayStyle = 'off';
end;

set(handles.user.aux_linehandles,'color',[0.5 0.5 0.5]); %reset to normal color

set(handles.user.sem_patchhandle,'visible','off');

handles.user.selected_handle = [];

%remove the inset axes and vertical line associated with the histogram
%distribution
% hideDistribution(handles);

function pri_line_buttonDownFcn(hObject,~)
%what to do when the user clicks on a line...
handles = guidata(hObject);
showBusy(handles);
if(~strcmp(get(handles.fig,'selectiontype'),'alt')) %'alt' refers to a ctrl-left mouse click or a right mouse click - which will activate the uicontextmenu
    
    handles = clearSelection(handles);

    currentpoint = get(handles.axes_main,'currentpoint');
%     freq_index = round(currentpoint(1)*2)+1; %*2 because we have 0.5 hz resolution; +1 because matlab is 1 based and not zero based for their indexing
    x = round(currentpoint(1)*2)/2;
    handles.user.cur_freq = x;
    
    drawDistribution(handles);

    set(handles.tb_meta_filename,'string','Name: ','visible','on');
    set(handles.tb_meta_filedate,'string','Modified: ','visible','on');
    
    
    set(handles.uitable_meta,'data',[]);
    
    pos = get(handles.uitable_meta,'position');
    extent = get(handles.uitable_meta,'extent');
    set(handles.uitable_meta,'position',[pos(1:2), extent(3:4)]);
    
    set(handles.tb_channel_name,'string','Name: ','visible','on');
    set(handles.tb_channel_num,'string','Number: ','visible','on');
    set(handles.tb_channel_winlen,'string','Window Length (sec): ','visible','on');
    set(handles.tb_channel_fftlen,'string','FFT Length (samples): ','visible','on');
    set(handles.tb_channel_fftint,'string','FFT Interval (sec): ','visible','on');
    
    if(get(handles.cb_art_eeg,'value'))
        stage_count = handles.user.stage_count;
    else
        stage_count = handles.user.stage_count_without_artifact;
    end
    
    child = get(handles.axes_pie,'children');
    if(~isempty(child))
        delete(child);
    end;
    explode_vec = zeros(1,size(stage_count,2));
    
    stage_selection = get(handles.pop_stage,'string');
    
    switch(stage_selection{get(handles.pop_stage,'value')})
        case 'All Night'
        case 'Awake'
        case 'Ratio'
        case 'Difference'
        case 'Relative Difference'
        case 'All nonREM'
            explode_vec(1:3)=1;
        case 'Stage 1'
            explode_vec(1)=1;
        case 'Stage 2'
            explode_vec(2)=2;
        case 'Stage 3-4'
            explode_vec(3)=3;
        case 'Stage 5 (REM)'
            explode_vec(4)=4;
        otherwise
            disp 'line 866'
            stage_selection{get(handles.pop_stage,'value')}
            warndlg('This selection was not handled correctly in its callback function','woops');
    end
    
    labelStr = {'Stage I','Stage II','Stage III-IV','Stage V (REM)'};
    pie(handles.axes_pie,stage_count+.1,explode_vec);
    legend(handles.axes_pie,'show',labelStr,'Location','SouthEastOutside');
    
    saveUndoState(handles);
    updateLegend(handles);

    guidata(hObject,handles);
    
end
showReady(handles);

function updateLegend(handles)
legH = legend(handles.axes_main);
if(~isempty(legH))
    delete(legH);
end;
legend(handles.axes_main,'off');
legh= legend(handles.axes_main,'show','location','north');
leg_pos = get(legh,'position');
axes_pos = get(handles.axes_main,'position');
leg_pos(2) = sum(axes_pos([2,4]))*0.95-leg_pos(4);
leg_pos(1) = sum(axes_pos([1,3]))*0.35;
set(legh,'position',leg_pos,'interpreter','none','hittest','off');
% handles.axes_main_legend_h = legh;
% guidata(gcf,handles);

% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over axes background.
function axes_main_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pri_line_buttonDownFcn(hObject,eventdata);


% --------------------------------------------------------------------
function menu_file_settings_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.user.settings = settings();

if(handles.user.settings.mean)
    set(get(handles.axes_main,'ylabel'),'string','Mean Power (\muV^2/Hz)','fontsize',13);
elseif(handles.user.settings.median)
    set(get(handles.axes_main,'ylabel'),'string','Median Power (\muV^2/Hz)','fontsize',13);
end;

set(handles.user.linehandles(:),'visible','off');
set(handles.user.linehandles(1:handles.user.settings.lines_to_hold),'visible','on');

hAnnotation = get(handles.user.linehandles(:),'Annotation');
for k=1:handles.user.settings.lines_to_hold
    hAnnotation{k}.LegendInformation.IconDisplayStyle = 'on';
end;
for k=handles.user.settings.lines_to_hold+1:numel(hAnnotation)
    hAnnotation{k}.LegendInformation.IconDisplayStyle = 'off';
end;

plotCurrentSettingsData(handles);


% --------------------------------------------------------------------
function menu_file_load_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pathname = uigetdir(handles.user.dataPathname,'Select Directory with Data Files');
if(pathname)
    handles.user.dataPathname = pathname;
    try
        % Changed loadPSAstats to loadPSDstats to keep up with SEV filename
        % conventions.
        [handles.user.data, handles.user.column_names, handles.user.files] = loadPSDstats(handles.user.dataPathname);
    catch me
        showME(me);
    end
    guidata(hObject,handles);
end;


% --------------------------------------------------------------------
function menu_file_about_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('Sleep data obtained from Stanford Sleep Lab with permission.','About');

% --- Executes on button press in cb_y_log_scale.
function cb_y_log_scale_Callback(hObject, eventdata, handles)
% hObject    handle to cb_y_log_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_y_log_scales
showBusy(handles);
if(get(handles.cb_y_log_scale,'value'))
    set(handles.axes_main,'yscale','log');
else
    set(handles.axes_main,'yscale','linear');
end;
showReady(handles);


% --- Executes on selection change in pop_ratio_bottom.
function pop_ratio_bottom_Callback(hObject, eventdata, handles)
% hObject    handle to pop_ratio_bottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_ratio_bottom contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_ratio_bottom
showBusy(handles);
plotCurrentSettingsData(handles);

% --- Executes during object creation, after setting all properties.
function pop_ratio_bottom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_ratio_bottom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_ratio_top.
function pop_ratio_top_Callback(hObject, eventdata, handles)
% hObject    handle to pop_ratio_top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_ratio_top contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_ratio_top
showBusy(handles);
plotCurrentSettingsData(handles);

% --- Executes during object creation, after setting all properties.
function pop_ratio_top_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_ratio_top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_art_ocular.
function cb_art_ocular_Callback(hObject, eventdata, handles)
% hObject    handle to cb_art_ocular (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_art_ocular
handles = plotData(handles);
guidata(hObject,handles);

% --- Executes on button press in cb_age_50.
function cb_age_50_Callback(hObject, eventdata, handles)
% hObject    handle to cb_age_50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_age_50
handles = plotData(handles);
guidata(hObject,handles);


% --- Executes on button press in cb_age_50_65.
function cb_age_50_65_Callback(hObject, eventdata, handles)
% hObject    handle to cb_age_50_65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = plotData(handles);
guidata(hObject,handles);

% Hint: get(hObject,'Value') returns toggle state of cb_age_50_65


% --- Executes on button press in cb_age_65_plus.
function cb_age_65_plus_Callback(hObject, eventdata, handles)
% hObject    handle to cb_age_65_plus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_age_65_plus
showBusy(handles);
handles = plotData(handles);
showReady(handles);
guidata(hObject,handles);



% --- Executes on button press in cb_male.
function cb_male_Callback(hObject, eventdata, handles)
% hObject    handle to cb_male (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_male
showBusy(handles);
handles = plotData(handles);
showReady(handles);
guidata(hObject,handles);


% --- Executes on button press in cb_female.
function cb_female_Callback(hObject, eventdata, handles)
% hObject    handle to cb_female (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_female
showBusy(handles);
handles = plotData(handles);
showReady(handles);
guidata(hObject,handles);


% --- Executes on slider movement.
function slider_epochs2hold_Callback(hObject, eventdata, handles)
% hObject    handle to slider_epochs2hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
showBusy(handles);
cur_value = floor(get(hObject,'value'));
handles.user.epochs2hold = cur_value;
set(handles.edit_epochs2hold,'string',cur_value);
plotData(handles);

% --- Executes during object creation, after setting all properties.
function slider_epochs2hold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_epochs2hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_epochs2hold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_epochs2hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_epochs2hold as text
%        str2double(get(hObject,'String')) returns contents of edit_epochs2hold as a double
showBusy(handles);
epochs2hold = str2double(get(hObject,'String'));
if(epochs2hold>=get(handles.slider_epochs2hold,'min') &&epochs2hold<=get(handles.slider_epochs2hold,'max'))
    epochs2hold = ceil(epochs2hold);
    handles.user.epochs2hold = epochs2hold;
    set(handles.slider_epochs2hold,'value',epochs2hold);
    handles = plotData(handles);

else %this case occurs for negative numbers or strings that get converted to NaN's.
    set(hObject,'string',handles.user.epochs2hold);
end;

showReady(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_epochs2hold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_epochs2hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on fig and none of its controls.
function fig_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to fig (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

key=eventdata.Key;

if(strcmp(key,'leftarrow') && handles.user.distribution.show)
    handles = guidata(hObject);
    showBusy(handles);
    x = get(handles.user.distribution.vertical_linehandle,'xdata');

    %         x = round(currentpoint(1)*2)/2-1;
    x = x(1)-0.5;
    x = max(x,handles.user.min_freq); %don't go off to the left into negative frequences
    handles.user.cur_freq = x;

    drawDistribution(handles);
    saveUndoState(handles);

    drawnow();
    showReady(handles);
    guidata(hObject,handles);
elseif(strcmp(key,'rightarrow') && handles.user.distribution.show)
    handles = guidata(hObject);
    showBusy(handles);
    x = get(handles.user.distribution.vertical_linehandle,'xdata');
    x=x(1)+0.5;
    x = min(x,handles.user.max_freq); %don't go beyond the number of frequenceis we have.
    %         x = round(currentpoint(1)*2)/2+1;
    %         freq_index = min(handles.user.num_ticks*2-1,x(1)*2+2); %(there is a +1 for rightarrow and a +1 for one-base indexing that sums here)
    handles.user.cur_freq = x;
    
    drawDistribution(handles);
    saveUndoState(handles);

    drawnow();
    showReady(handles);
    guidata(hObject,handles);
end;

%nuke it
if(strcmp(eventdata.Key,'c') && strcmp(eventdata.Modifier,'control'))
    delete(hObject);
end;

%take screen capture of figure
if(strcmp(eventdata.Key,'p') && strcmp(eventdata.Modifier,'control'))
    screencap(gcf);
    disp('screen capture');   
end;

%take screen capture of axes
if(strcmp(eventdata.Key,'a') && strcmp(eventdata.Modifier,'control'))
    screencap(handles.uipanel_axes);
    disp('screen capture');
end;

%undo last step...
if(strcmp(eventdata.Key,'z') && strcmp(eventdata.Modifier,'control'))
    handles = guidata(hObject);
    showBusy(handles);

    handles = undo(handles);
    save_undo_state = false;
    plotCurrentSettingsData(handles,save_undo_state);
    drawnow();
    showReady(handles);
    guidata(hObject,handles);
end;

function showBusy(handles)
set(handles.text_working,'visible','on');
drawnow();

function showReady(handles)
set(handles.text_working,'visible','off');
drawnow();

% --- Executes during object creation, after setting all properties.
function uipanel_settings_duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_settings_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

uicontextmenu_handle = uicontextmenu();%,get(parentAxes,'parent'));
uimenu(uicontextmenu_handle,'Label','Minimum Duration','separator','off','callback',@min_duration_callback)
uimenu(uicontextmenu_handle,'Label','Maximum Duration','separator','off','callback',@max_duration_callback);
set(hObject,'uicontextmenu',uicontextmenu_handle);

function min_duration_callback(hObject,eventdata)
handles = guidata(hObject);
showBusy(handles);
handles.user.min_duration = true;
handles.user.max_duration  = false;
set(handles.uipanel_settings_duration,'title','Minimum Duration');
handles = plotCurrentSettingsData(handles);
guidata(hObject,handles);

function max_duration_callback(hObject,eventdata)
handles = guidata(hObject);
showBusy(handles);
handles.user.min_duration = false;
handles.user.max_duration = true;
set(handles.uipanel_settings_duration,'title','Maximum Duration');
handles = plotCurrentSettingsData(handles);
guidata(hObject,handles);


% --------------------------------------------------------------------
function Settings_Duration_contextmenu_Callback(hObject, eventdata, handles)
% hObject    handle to Settings_Duration_contextmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function duration_minimum_Callback(hObject, eventdata, handles)
% hObject    handle to duration_minimum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function duration_maximum_Callback(hObject, eventdata, handles)
% hObject    handle to duration_maximum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_edit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_edit_undo_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit_undo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = undo(handles);
save_for_undo = false;
plotCurrentSettingsData(handles,save_for_undo);
% plotData(handles);


% --------------------------------------------------------------------
function show_histogram_menu_Callback(hObject, eventdata, handles)
% hObject    handle to show_histogram_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(strcmp(get(hObject,'checked'),'on'))
    set(hObject,'checked','off');
    hideDistribution(handles);
    handles.user.distribution.show = false;
else
    set(hObject,'checked','on');
%     x = get(handles.user.distribution.vertical_linehandle,'xdata');
    
    handles.user.distribution.show = true;
    drawDistribution(handles);
end;
saveUndoState(handles);

updateLegend(handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function Main_Axes_contextmenu_Callback(hObject, eventdata, handles)
% hObject    handle to Main_Axes_contextmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function main_axes_cmenu_filter_menu_Callback(hObject, eventdata, handles)
% hObject    handle to main_axes_cmenu_filter_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmenu_filter_mean_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_filter_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(strcmp(get(hObject,'checked'),'on'))
    set(hObject,'checked','off');
    set(handles.cmenu_filter_median,'checked','on');
    handles.user.settings.median = true;
    handles.user.settings.mean = false;
else
    set(hObject,'checked','on');
    set(handles.cmenu_filter_median,'checked','off');
    handles.user.settings.median = false;
    handles.user.settings.mean = true;
end;

if(handles.user.settings.mean)
    set(get(handles.axes_main,'ylabel'),'string','Mean Power (\muV^2/Hz)','fontsize',13);
elseif(handles.user.settings.median)
    set(get(handles.axes_main,'ylabel'),'string','Median Power (\muV^2/Hz)','fontsize',13);
end;

showBusy(handles);

if(handles.user.distribution.show)
    handles.user.distribution.show = false;
    hideDistribution(handles);
    handles = clearSelection(handles);
    saveforundo = false;
    handles = plotCurrentSettingsData(handles,saveforundo);
    handles.user.distribution.show = true;
    saveUndoState(handles);

else
    handles = clearSelection(handles);
    handles = plotCurrentSettingsData(handles);
end
showReady(handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function cmenu_filter_median_Callback(hObject, eventdata, handles)
% hObject    handle to cmenu_filter_median (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(strcmp(get(hObject,'checked'),'on'))
    set(hObject,'checked','off');
    set(handles.cmenu_filter_mean,'checked','on');
    handles.user.settings.mean = true;
    handles.user.settings.median = false;
else
    set(hObject,'checked','on');
    set(handles.cmenu_filter_mean,'checked','off');
    handles.user.settings.mean = false;
    handles.user.settings.median = true;
end;

if(handles.user.settings.mean)
    set(get(handles.axes_main,'ylabel'),'string','Mean Power (\muV^2/Hz)','fontsize',13);
elseif(handles.user.settings.median)
    set(get(handles.axes_main,'ylabel'),'string','Median Power (\muV^2/Hz)','fontsize',13);
end;

showBusy(handles);
if(handles.user.distribution.show)
    handles.user.distribution.show = false;
    hideDistribution(handles);
    handles = clearSelection(handles);
    saveforundo = false;
    handles = plotCurrentSettingsData(handles,saveforundo);
    handles.user.distribution.show = true;
    saveUndoState(handles);

else
    handles = clearSelection(handles);
    handles = plotCurrentSettingsData(handles);
end
showReady(handles);
guidata(hObject,handles);
