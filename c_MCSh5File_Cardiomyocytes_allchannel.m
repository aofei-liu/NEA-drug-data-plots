%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NAME:   
% %     c_MCSh5File_Cardiomyocytes_allchannel
% % PURPOSE: 
% %     A data processing software package with a graphic user interface.
% % CATEGORY:
% %     Data processing.
% % CALLING SEQUENCE:
% %     c_MCSh5File_Cardiomyocytes_allchannel
% % INPUTS:
% %     Put in various parameters according to the intruction.
% % OUTPUS:
% %     Display data.
% %     Analyze action potential parameters.
% % COMMENTS:
% %     none.
% % HISTORY:
% %     Written by Bianxiao Cui on September 6th, 2019
% %     Updated by Bianxiao Cui on September 12th, 2019.
% %     Updated by Aofei on September 12th, 2019 (added 2 save functions)
% %     Updated by Aofei on September 19, 2019 (modified save functions)
% %     Updated by Bianxiao Cui on May 12th, 2020.
% %     Updated by Bianxiao Cui on July 11th, 2020. (updated plotting range)
% %     Updated by Bianxiao Cui on July 23rd, 2020. (Added automated peak finding)
% %     Updated by Bianxiao Cui on Oct. 30, 2020. (Improved
% autoPeakDetection).
% %     Updated by Ethan Foster on Nov 12th, 2020 (modified to accept h5 files)
% %     Updated by Ethan Foster on Nov 27th, 2020 (modified autoPeakDetection)
% %     Updated by Aofei on Feb 7th, 2022 (added batch analyze method)
% %     Updated by Aofei on Apr 5th, 2022 (modified peak detection methods)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = c_MCSh5File_Cardiomyocytes_allchannel(varargin)
% C_MCSH5FILE_CARDIOMYOCYTES_ALLCHANNEL M-file for c_MCSh5File_Cardiomyocytes_allchannel.fig
%      C_MCSH5FILE_CARDIOMYOCYTES_ALLCHANNEL, by itself, creates a new C_MCSH5FILE_CARDIOMYOCYTES_ALLCHANNEL or raises the existing
%      singleton*.
%
%      H = C_MCSH5FILE_CARDIOMYOCYTES_ALLCHANNEL returns the handle to a new C_MCSH5FILE_CARDIOMYOCYTES_ALLCHANNEL or the handle to
%      the existing singleton*.
%
%      C_MCSH5FILE_CARDIOMYOCYTES_ALLCHANNEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in C_MCSH5FILE_CARDIOMYOCYTES_ALLCHANNEL.M with the given input arguments.
%
%      C_MCSH5FILE_CARDIOMYOCYTES_ALLCHANNEL('Property','Value',...) creates a new C_MCSH5FILE_CARDIOMYOCYTES_ALLCHANNEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before c_cdf_process_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to c_MCSh5File_Cardiomyocytes_allchannel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help c_MCSh5File_Cardiomyocytes_allchannel

% Last Modified by GUIDE v2.5 07-Feb-2022 18:36:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @c_MCSh5File_Cardiomyocytes_allchannel_OpeningFcn, ...
                   'gui_OutputFcn',  @c_MCSh5File_Cardiomyocytes_allchannel_OutputFcn, ...
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


% --- Executes just before c_MCSh5File_Cardiomyocytes_allchannel is made visible.
function c_MCSh5File_Cardiomyocytes_allchannel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to c_MCSh5File_Cardiomyocytes_allchannel (see VARARGIN)

% % %     axes(handles.axes1);    %% make the existing axes1 the current axes.
% % %     cla;        %% clear current axes
% % %     set(gcf,'Units','pixels');
% % %     [x,map]=imread('yxz.tif');
% % %     imshow(x(1:512,1:512,:),map);

% Choose default command line output for c_MCSh5File_Cardiomyocytes_allchannel
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes c_MCSh5File_Cardiomyocytes_allchannel wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = c_MCSh5File_Cardiomyocytes_allchannel_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in Open.
function Open_Callback(hObject, eventdata, handles)
% hObject    handle to Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    [datafile,pathname]=uigetfile('*.h5','Choose an h5 file'); % read cui-lab data format files
    if datafile==0, return; end
    cd(pathname);
    
    file_parts = str2num(get(handles.num_files, 'String'));
    for i = 1:file_parts
        if i == 1
            filenames{i} = datafile;
        else
            name = datafile(1:strfind(datafile,'.')-1);
            ext = '.h5';
            filenames{i} = append(name, '000', num2str(i-1), ext);
        end
    end
    DataArray = [];
    for i = 1:file_parts
        try
            DataArray_temp=h5read(filenames{i},'/dataArray');
            DataArray = [DataArray,DataArray_temp];
            set(handles.FileName,'String',{filenames{i};...
                'Cyion File'});
            handles.fileType='Cyion';
        catch
            set(handles.FileName,'String',{filenames{i};...
                'MCS File'});
            handles.fileType='MCS';
            channelInfo=h5read(filenames{i},'/Data/Recording_0/AnalogStream/Stream_0/InfoChannel');
            DataArray_temp=h5read(filenames{i},'/Data/Recording_0/AnalogStream/Stream_0/ChannelData');
            DataArray = vertcat(DataArray, DataArray_temp);
        end
    end
    sz = size(DataArray);
    xsize = sz(1);
    ysize = sz(2);
    Nchannels = ysize;
    % Sampling time of DAQ card 
    SamplingInterval = double(channelInfo.Tick(1));
    sampling_time = SamplingInterval * 10 ^-6; %%SamplingInterval unit is microsecond
    sampling_frequency = 1./sampling_time;
    time_end =(sampling_time*xsize)- sampling_time;
    TimeArray = [0:sampling_time:time_end];
    TimeArray = permute(TimeArray,[2,1]);         
    set(handles.DisplayText,'String',{['sampling frequency (Hz):  ' num2str(sampling_frequency)];...
       ['  time duration(s):  ' num2str(time_end)]});

        
    handles.filename=filenames{1};
    handles.datapath=pathname;
    handles.DataArrayAllChannels=DataArray;
    handles.TimeArray=TimeArray;
    handles.channelInfo=channelInfo; 
    handles.sampling_time = sampling_time; %% s unit
    handles.sampling_frequency = sampling_frequency; 
   guidata(hObject,handles);

    
     % --- Executes on button press in PlotAllChannels.
function PlotAllChannels_Callback(hObject, eventdata, handles)
    % hObject    handle to PlotAllChannels (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    fig = figure(1);
    clf;
    channelInfo = handles.channelInfo;        
    TimeArray = handles.TimeArray;
    filename=handles.filename;

    startTime = str2num(get(handles.StartTime, 'String'));
    endTime = str2num(get(handles.EndTime, 'String')); 
    [temp, startIndex] = min(abs(TimeArray - startTime));
    [temp, endIndex] = min(abs(TimeArray - endTime));
    TimeArray_select = TimeArray(startIndex:endIndex);
    titleFilename = handles.FileName.String{1};
    sgtitle({strrep(append(titleFilename, '    ', string(startTime), ' sec - ', string(endTime), ' sec'), '_', '\_'), ''});   
    %%sgtitle({strrep(append('All channels:', '    ', string(startTime), ' sec - ', string(endTime), ' sec'), '_', '\_'), ''});

        for i = 1:60
            DataArray = double(handles.DataArrayAllChannels(:, i))*10^-3 * 2.3585; %%2.3585 factor determined by the company
            channelLabel = channelInfo.Label{i};
            subplotRow = str2num(channelLabel(2));
            subplotCol = str2num(channelLabel(1));
            subplotPos = (subplotRow - 1) * 8 + subplotCol;
            
            %determine if data is good for peak finding
            DataArray_select = DataArray(startIndex:endIndex);
            FilteredData=smooth(DataArray_select,5);
            stdData = movstd(smooth(DataArray_select, 5),20);
            spikeIntervalThreshold = 0.4; 
            stdThreshold = max(stdData)/3;
            [pks_temp,spike_locs_temp]=findpeaks(stdData,'MinPeakProminence',stdThreshold,...
            'MinPeakDistance',spikeIntervalThreshold*handles.sampling_frequency); 


            %% plotonechannel out the data
            hold on

            if subplotPos<8   
                %%if subplotPos == 0
                if subplotPos == 7  %%% this channel is always a bad channel and set to 0. 
                    subplot(8, 8, subplotPos)
                    plot(TimeArray(startIndex:endIndex),0*DataArray(startIndex:endIndex));
                    set( findobj(gca,'typehFig','line'), 'LineWidth', 1);
                    xticks([]);
                    ylim([-1, 1]);
                    t = title(append(string(channelLabel), '          '));

                else
                    subplot(8, 8, subplotPos)
                    plot(TimeArray(startIndex:endIndex),DataArray(startIndex:endIndex));
                    hold on;
                    % plot(TimeArray_select(spike_locs_temp), DataArray_select(spike_locs_temp), 'r*');
                    set( findobj(gca,'typehFig','line'), 'LineWidth', 1);
                    xticks([]);
                    ylim('auto');
                    temp=[min(DataArray(startIndex:endIndex)), max(DataArray(startIndex:endIndex))];
                    extended_yrange = 0.2*(temp(2)-temp(1));
                    extended_yrange = max(extended_yrange, 0.25);
                    ylim([temp(1)-extended_yrange temp(2)+extended_yrange]);
                    t = title(append(string(channelLabel), '          '));
                end
            elseif subplotPos>8 & subplotPos<57
                
                if subplotPos-2 == 13 | subplotPos-2 == 48  %%these are two pacing channels. 
                    subplot(8, 8, subplotPos)
                    plot(TimeArray(startIndex:endIndex),0*DataArray(startIndex:endIndex));
                    set( findobj(gca,'typehFig','line'), 'LineWidth', 1);
                    xticks([]);
                    ylim([-1, 1]);
                    t = title(append(string(channelLabel), '          '));
                else
                    subplot(8, 8, subplotPos)
                    plot(TimeArray(startIndex:endIndex),DataArray(startIndex:endIndex));
                    hold on;
                    % plot(TimeArray_select(spike_locs_temp), DataArray_select(spike_locs_temp), 'r*');
                    set( findobj(gca,'typehFig','line'), 'LineWidth', 1);
                    xticks([]);
                    ylim('auto');
                    temp=[min(DataArray(startIndex:endIndex)), max(DataArray(startIndex:endIndex))];
                    extended_yrange = 0.2*(temp(2)-temp(1));
                    extended_yrange = max(extended_yrange, 0.25);
                    ylim([temp(1)-extended_yrange temp(2)+extended_yrange]);
                    t = title(append(string(channelLabel), '          '));              
                end
                
            elseif subplotPos>57
                    subplot(8, 8, subplotPos)
                    plot(TimeArray(startIndex:endIndex),DataArray(startIndex:endIndex));
                    hold on;
                    % plot(TimeArray_select(spike_locs_temp), DataArray_select(spike_locs_temp), 'r*');
                    set( findobj(gca,'typehFig','line'), 'LineWidth', 1);
                    xticks([]);
                    ylim('auto');
                    temp=[min(DataArray(startIndex:endIndex)), max(DataArray(startIndex:endIndex))];
                    extended_yrange = 0.2*(temp(2)-temp(1));
                    extended_yrange = max(extended_yrange, 0.25);
                    ylim([temp(1)-extended_yrange temp(2)+extended_yrange]);
                    t = title(append(string(channelLabel), '          '));
            end
            t.Units = 'Normalize'; 
            t.Position(1) = 0;
            t.Position(2) = 1;
            t.HorizontalAlignment = 'left';  
        end
   set(handles.DisplayText,'String',{
            ['sampling frequency (Hz):  ' num2str(5000)];...
            ['  time duration(s):  ' num2str(TimeArray(length(TimeArray)))]});
        
   guidata(hObject,handles);
   
   % --- Executes on button press in PlotOneChannel.
function PlotOneChannel_Callback(hObject, eventdata, handles)
% hObject    handle to PlotOneChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    channelInfo = handles.channelInfo;        
    TimeArray = handles.TimeArray;
    filename=handles.filename;
% % %         channelLabel = str2num(get(handles.ChannelNumber,'String'));
% % %         channelIndex=find(str2num(cell2mat(channelInfo.Label)) == channelLabel);
% % %         channelID = channelIndex - 1;        
% % %         DataArray = double(handles.DataArrayAllChannels(:, channelIndex))*10^-3;
% % %         DataArray = DataArray * 2.3585; %%Scaling factor to match abf
% % %    set(handles.DisplayText,'String',{['sampling frequency (Hz):  ' num2str(5000)];...
% % %             ['  time duration(s):  ' num2str(TimeArray(length(TimeArray)))]});



        channelLabel = str2num(get(handles.ChannelNumber,'String'));
        channelIndex=find(str2num(cell2mat(channelInfo.Label)) == channelLabel);
        channelID = channelIndex - 1;        
        DataArray = double(handles.DataArrayAllChannels(:, channelIndex))*10^-3;
        DataArray = DataArray * 2.3585; %%Scaling factor to match abf
   set(handles.DisplayText,'String',{['sampling frequency (Hz):  ' num2str(5000)];...
            ['  time duration(s):  ' num2str(TimeArray(length(TimeArray)))]});
        

    %% plotonechannel out the data
    axes(handles.axes1);
    hold off;
     plot(TimeArray,DataArray);
        set( findobj(gca,'typehFig','line'), 'LineWidth', 1);
        set(gca, 'FontName', 'Arial')
        set(gca, 'FontSize', 8)
        xlabel('time (s)')
        ylabel('\DeltaV (mV)')
        
    handles.DataArray=DataArray;
    handles.DataSize = numel(DataArray);
    handles.axrange=axis(handles.axes1);
    handles.channelLabel=channelLabel;
    handles.channelIndex=channelIndex;
   guidata(hObject,handles);
 
   
% --- Executes on button press in AutoPeakDetection.
function AutoPeakDetection_Callback(hObject, eventdata, handles)
% hObject    handle to AutoPeakDetection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    DataArray = handles.DataArray;
    TimeArray = handles.TimeArray;
    sampling_time = handles.sampling_time;
    sampling_frequency = handles.sampling_frequency;

    axes(handles.axes1);
    xrange = xlim;
    starting_position = max(1,fix(xrange(1)/sampling_time));
    ending_position = min(handles.DataSize, fix(xrange(2)/sampling_time));
    DataArray_select = DataArray(starting_position:ending_position);
    TimeArray_select = TimeArray(starting_position:ending_position);
    NumPoints =  numel(DataArray_select);
    
    %%%%%%% develop a spike detection method %%%%%%%%%%%%
    %%%%%%% The algothrim automatically detect the peak and set the
    %%%%%%% threshold.  It uses a sliding window to adjust the threshold
    %%%%%%% in a time variable fashion. 
        slidingWindowInterval = str2num(get(handles.SlidingWindow,'String'));
        FilteredData=smooth(DataArray_select,slidingWindowInterval); %% 5-point average. The small range is to avoid averaging-induced shift in the maximum point.
        stdRange = 20; %% sliding window smoothing the data 
        stdData = movstd(FilteredData,stdRange);
        spikeIntervalThreshold = str2num(get(handles.SpikeIntervalThreshold,'String'));
        spikeIntervalThreshold = max(spikeIntervalThreshold, 0.1);  %% highest frequency 10
        stdThreshold = str2num(get(handles.StdThreshold,'String'));
        [pks_temp,spike_locs_temp]=findpeaks(stdData,'MinPeakProminence',stdThreshold,...
            'MinPeakDistance',spikeIntervalThreshold*sampling_frequency); 
        [max_std, max_std_loc] = max(stdData); % if porating, identifies time poration starts.

        slidingwindow_size = sampling_frequency*5; %%for every five seconds, re-calibrate the threshold
        spike_locs=[];
        for i=1:slidingwindow_size:numel(stdData)
            rightside_index = min(i+slidingwindow_size, numel(stdData));%%right side of the window
            maxstd_window = max(stdData(i:rightside_index));
            std_threshold_window = maxstd_window/3;
            spike_locs_window = find(spike_locs_temp>i & spike_locs_temp<rightside_index);
            pks_std_window = pks_temp(spike_locs_window);
            spike_locs_identified = spike_locs_window(find(pks_std_window>std_threshold_window));
            spike_locs = cat(1,spike_locs,spike_locs_temp(spike_locs_identified)); 
        end
        spike_locs = spike_locs(2:end); %%this line is to take care of corner cases when the selection starts at spike location of the AP.

    axes(handles.axes1);
    xrange = xlim;
    hold off;
    plot(TimeArray,DataArray,'b');
    hold on;
    plot(TimeArray_select,DataArray_select,'k');
    plot(TimeArray_select(spike_locs),DataArray_select(spike_locs),'*r');
    hold off;
    set( findobj(gca,'typehFig','line'), 'LineWidth', 1);
    set(gca, 'FontName', 'Arial')
    set(gca, 'FontSize', 8)
    xlabel('time (s)')
    ylabel('\DeltaV (mV)')

    axes(handles.axes2);
    hold off;
    plot(TimeArray_select,stdData,'k');
    hold on;
    plot(TimeArray_select(spike_locs),stdData(spike_locs),'r*');
    hold off;
    set( findobj(gca,'typehFig','line'), 'LineWidth', 1);
    set(gca, 'FontName', 'Arial')
    set(gca, 'FontSize', 8)
    xlabel('time (s)')
    ylabel('Standard Deviation')

    avg_frequency = (numel(spike_locs)-1)/...
        (TimeArray_select(spike_locs(end))- TimeArray_select(spike_locs(1)));
    set(handles.DisplayText,'String',{['# of spiking peaks:  ' num2str(numel(spike_locs))];...
       ['  for time (s) from:  ' num2str(TimeArray(starting_position))  '  to  '...
       num2str(TimeArray(ending_position))];
       ['  spiking frequency (Hz):  ' num2str(avg_frequency)];
       [' starting time(s) for poration (if applicable):' num2str(max_std_loc*sampling_time)]});

    handles.stdData_select=stdData;
    handles.DataArray_select = DataArray_select;
    handles.TimeArray_select = TimeArray_select;
    handles.spike_locs = spike_locs;
    handles.Data_select_baseline = DataArray_select*0;


   guidata(hObject,handles);
   
 
 
% --- Executes on button press in Analyze.
function Analyze_Callback(hObject, eventdata, handles)
% hObject    handle to Analyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    DataArray_select = handles.DataArray_select;
    TimeArray_select = handles.TimeArray_select;
    stdData_select=handles.stdData_select;
    TimeArray_select = handles.TimeArray_select;
    sampling_time = handles.sampling_time;
    spike_locs = handles.spike_locs;

        % baseline if available
        bas = isfield(handles, 'Data_select_baseline');
        if bas == 1
            DataArray_select = handles.DataArray_select - handles.Data_select_baseline;
        end

    AvgRange=fix(2*10^(-3)/handles.sampling_time); %%average over 2ms data
    FilteredData_select=smooth(DataArray_select,AvgRange); 

    clear max_locs max_vals min_locs min_vals start_locs start_vals max_min_ratio;
    clear APD100 APD90 APD50 APD10 SpikeVelocity SpikeInterval SpikeDuration;

    stdmedian = median(stdData_select);
    sortedarray=sort(stdData_select,'ascend');
    DataSize = numel(sortedarray);
    stdofstd = std(sortedarray(1:fix(DataSize*0.6)));

    for i=1:numel(spike_locs)-1
        %%finding the maximum and miniumum points position for each AP
        periodData = FilteredData_select(spike_locs(i):spike_locs(i+1));
        [value, locs] = max(periodData);
        max_locs(i) = spike_locs(i)+locs-1;
        max_vals(i) = value;
%%this part of the code is changed on 2022-04-05.  The replacement code to find the end of an action
% potential is after the APD calculations. The change is to take
% consideration of traces with low signal to noise ratio and for traces
% without clear minimum over potential. 
% %         [value, locs]=min(periodData);
% %         min_locs(i) = spike_locs(i)+locs-1;
% %         min_vals(i) = value;    

        %%finding the start point position for each AP
        start_pos=spike_locs(i); 
        while start_pos>=1 && stdData_select(start_pos)> stdmedian; %%replacing stdmedian+3*stdofstd; this is changed on 2022-04-02
            start_pos=start_pos-1;
        end
        start_locs(i)=start_pos;
                         %%start_vals(i)=FilteredData_select(start_pos);
                         %%%%modifed on 2022-04-02
        periodPoints = numel(periodData); %% number of data points in a period cycle
        avg_baseline = median(FilteredData_select(start_pos-fix(periodPoints*0.05):start_pos));
        start_vals(i)= avg_baseline; %%this is changed on 2022-04-02 to account for the slight variations in start_time identification. 

        %%Added on 2022-4-05. calculate the end of an action potential, approximately the
        %%minumum. The variable is called min_locs and min_vals. They work for traces that do not have a
        %%clear minimum.  
            Half_Amp = start_vals(i)+(max_vals(i)-start_vals(i))*0.1; %%APD90
            Half_Amp_points = numel(find(periodData>Half_Amp));  %%number of points with values greater than Half_amp
            if Half_Amp_points>1
                Range50ms = 50*10^(-3)/sampling_time;
                stdData50ms = movstd(periodData,Range50ms);
                range_start = Half_Amp_points;
                range_end = min(fix(Half_Amp_points*3), numel(periodData));
                med_stdData50ms = median(stdData50ms(range_start:range_end));
                    APend_pos = Half_Amp_points;
                %%APend_pos means the end of the repolarization phase
                while APend_pos>=1 && stdData50ms(APend_pos)> med_stdData50ms 
                    APend_pos=APend_pos+1;
                end
                min_loc_1 = spike_locs(i)+APend_pos-1;
                min_val_1=FilteredData_select(min_loc_1);
                [value, locs] = min(periodData);
                min_loc_2 = spike_locs(i)+locs-1;
                min_val_2 = value;
                if min_loc_2 < min_loc_1
                    min_vals(i) = min_val_2;
                    min_locs(i) = min_loc_2;
                else
                    min_vals(i) = min_val_1;
                    min_locs(i) = min_loc_1;
                end
            else
                min_locs(i) = spike_locs(i)+1; %%strange results, purposely set min to max
                min_vals(i)=FilteredData_select(min_locs(i));
            end  

        %%calculating the spike velocity and interval
        temp_data=FilteredData_select(spike_locs(i)-2:spike_locs(i)+2); %%unit mV
        temp_time=TimeArray_select(spike_locs(i)-2:spike_locs(i)+2); %% unit Sec
        temp_fit=polyfit(temp_time,temp_data,1);
        SpikeVelocity(i)=temp_fit(1);  %%SpikeVelocity is usually useless as the amplitude is not scaled. 
        SpikeInterval(i) = TimeArray_select(spike_locs(i+1))-TimeArray_select(spike_locs(i)); %%SpikeInterval is Cycle Length

        AmpSM(i) = max_vals(i)-start_vals(i); %%amplitude start to max
        AmpPP(i) = max_vals(i)-min_vals(i);   %%amplitude peak to peak
        PeriodData = FilteredData_select(start_locs(i):min_locs(i)); %%redefine the period data
        APD100(i) = numel(find(PeriodData > start_vals(i)))*sampling_time;
        APD90(i) = numel(find(PeriodData > (start_vals(i)+AmpSM(i)*0.1)))*sampling_time;
        APD80(i) = numel(find(PeriodData > (start_vals(i)+AmpSM(i)*0.2)))*sampling_time;
        APD70(i) = numel(find(PeriodData > (start_vals(i)+AmpSM(i)*0.3)))*sampling_time;
        APD60(i) = numel(find(PeriodData > (start_vals(i)+AmpSM(i)*0.4)))*sampling_time;
        APD50(i) = numel(find(PeriodData > (start_vals(i)+AmpSM(i)*0.5)))*sampling_time;
        APD40(i) = numel(find(PeriodData > (start_vals(i)+AmpSM(i)*0.6)))*sampling_time; 
        APD30(i) = numel(find(PeriodData > (start_vals(i)+AmpSM(i)*0.7)))*sampling_time; 
        APD20(i) = numel(find(PeriodData > (start_vals(i)+AmpSM(i)*0.8)))*sampling_time; 
        APD10(i) = numel(find(PeriodData > (start_vals(i)+AmpSM(i)*0.9)))*sampling_time; 

        %%calculating the time it takes from APD80 to APD20
        PeriodData = FilteredData_select(start_locs(i):max_locs(i)); %%redefine the period data
        t_rise(i) = numel(find(PeriodData > (start_vals(i)+AmpSM(i)*0.2) &...
            PeriodData < (start_vals(i)+AmpSM(i)*0.8)))*sampling_time;
        %%calculate triangulation of data (APD90-APD30 in the phase after
        %%the peak) and the max to min/cycle time ratio
        RepolarizationData = FilteredData_select(max_locs(i):min_locs(i)); %%define RepolarizationData
        RT_CL_ratio(i) = numel(RepolarizationData)*sampling_time/SpikeInterval(i);
        triangulation(i) = numel(find(RepolarizationData > (start_vals(i)+AmpSM(i)*0.1) &...
            RepolarizationData < (start_vals(i)+AmpSM(i)*0.7)))*sampling_time;

    end


    axes(handles.axes1);
    hold off;
    plot(TimeArray_select,DataArray_select,'k');
    hold on
    plot(TimeArray_select(spike_locs),DataArray_select(spike_locs),'r*');
    plot(TimeArray_select(max_locs),max_vals,'b*');
    plot(TimeArray_select(min_locs),min_vals,'c*');
    plot(TimeArray_select(start_locs),start_vals,'g*');
    hold off;
    set( findobj(gca,'typehFig','line'), 'LineWidth', 1);
    set(gca, 'FontName', 'Arial')
    set(gca, 'FontSize', 8)
    xlabel('time (s)')
    ylabel('\DeltaV (mV)')


    axes(handles.axes2);
    hold off;
    plot(TimeArray_select,stdData_select,'k');
    hold on;
    plot(TimeArray_select(spike_locs),stdData_select(spike_locs),'r*');
    plot(TimeArray_select(max_locs),stdData_select(max_locs),'b*');
    plot(TimeArray_select(min_locs),stdData_select(min_locs),'c*');
    plot(TimeArray_select(start_locs),stdData_select(start_locs),'g*');
    hold off;
    set( findobj(gca,'typehFig','line'), 'LineWidth', 1);
    set(gca, 'FontName', 'Arial')
    set(gca, 'FontSize', 8)
    xlabel('time (s)')
    ylabel('Standard Deviation')


    clear AnalyzedData;
    AnalyzedData(1).name = 'OrigData'; %%Original data, the selected time segment
    AnalyzedData(2).name = 'TimeData'; %%Time axins, the selected time segment
    AnalyzedData(3).name = 'StdData';  %%Standard deviation of the original data
    AnalyzedData(4).name = 'HighestPointLocation'; %%Array locations of the highest point
    AnalyzedData(5).name = 'LowestPointLocation';  %%Array locations of the lowest point
    AnalyzedData(6).name = 'StartPointLocation';   %%Array locations of the peak starting point
    AnalyzedData(7).name = 'SpikePointLocation';   %%Array locations of the spike point
    AnalyzedData(8).name = 'Amplitude'; %%from starting point to minimum
    AnalyzedData(9).name = 'AmplitudePeaktoPeak'; %%from maximum to minimum
    AnalyzedData(10).name = 'SpikeInterval'; %%Time intervals between neighboring peaks (s)
    AnalyzedData(11).name = 'SpikeFrequency'; %%1/SpikeInterval
    AnalyzedData(12).name = 'APD100';
    AnalyzedData(13).name = 'APD90';
    AnalyzedData(14).name = 'APD50';
    AnalyzedData(15).name = 'APD10';
    AnalyzedData(16).name = 'SpikeVelocity'; %%Slope at the maximum point
    AnalyzedData(17).name = 't-rise'; %%from APD80 to APD20
    AnalyzedData(18).name = 'APD30';
    AnalyzedData(19).name = 'APD40';
    AnalyzedData(20).name = 'APD60';
    AnalyzedData(21).name = 'APD20';
    AnalyzedData(22).name = 'APD70';
    AnalyzedData(23).name = 'APD80';
    AnalyzedData(24).name = 'triangulation';
    AnalyzedData(25).name = 'RepolarizationTime/CycleLength ratio';

    AnalyzedData(1).data = DataArray_select;
    AnalyzedData(2).data = TimeArray_select;
    AnalyzedData(3).data = stdData_select;
    AnalyzedData(4).data = max_locs;
    AnalyzedData(5).data = min_locs;
    AnalyzedData(6).data = start_locs;
    AnalyzedData(7).data = spike_locs;
    AnalyzedData(8).data = AmpSM; %%Amplitude from start to max
    AnalyzedData(9).data = AmpPP; %%Amplitude peak to peak
    AnalyzedData(10).data = SpikeInterval;
    AnalyzedData(11).data = 1./SpikeInterval;
    AnalyzedData(12).data = APD100;
    AnalyzedData(13).data = APD90;
    AnalyzedData(14).data = APD50;
    AnalyzedData(15).data = APD10;
    AnalyzedData(16).data = SpikeVelocity;
    AnalyzedData(17).data = t_rise;
    AnalyzedData(18).data = APD30;
    AnalyzedData(19).data = APD40;
    AnalyzedData(20).data = APD60;
    AnalyzedData(21).data = APD20;
    AnalyzedData(22).data = APD70;
    AnalyzedData(23).data = APD80;
    AnalyzedData(24).data = triangulation;
    AnalyzedData(25).data = RT_CL_ratio;

    APDratio = mean(APD50./APD90);
    Dip = AmpPP-AmpSM;
    DipAmpratio = mean(Dip/AmpSM);
    if APDratio>0.5
        celltype = 'ventricular';
    else
        celltype = 'atrial';
    end
    set(handles.DisplayText,'String',{['Avg APD50/APD90 ratio:  ' num2str(APDratio)];...
       ['  Avg Dip/Amp ratio:  ' num2str(DipAmpratio)]; ...
       ['  Cell type:  ' celltype]});

    handles.AnalyzedData = AnalyzedData;
     guidata(hObject,handles);
     
 
    
function DisplayText_Callback(hObject, eventdata, handles)
% hObject    handle to DisplayText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DisplayText as text
%        str2double(get(hObject,'String')) returns contents of DisplayText as a double


% --- Executes during object creation, after setting all properties.
function DisplayText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DisplayText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function FileName_Callback(hObject, eventdata, handles)
% hObject    handle to FileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FileName as text
%        str2double(get(hObject,'String')) returns contents of FileName as a double


% --- Executes during object creation, after setting all properties.
function FileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ChannelNumber_Callback(hObject, eventdata, handles)
% hObject    handle to SetMaxPeakWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SetMaxPeakWidth as text
%        str2double(get(hObject,'String')) returns contents of SetMaxPeakWidth as a double

% --- Executes during object creation, after setting all properties.
function ChannelNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChannelNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    

function StartTime_Callback(hObject, eventdata, handles)
% hObject    handle to StartTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StartTime as text
%        str2double(get(hObject,'String')) returns contents of StartTime as a double


% --- Executes during object creation, after setting all properties.
function StartTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StartTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EndTime_Callback(hObject, eventdata, handles)
% hObject    handle to EndTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EndTime as text
%        str2double(get(hObject,'String')) returns contents of EndTime as a double


% --- Executes during object creation, after setting all properties.
function EndTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EndTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%%%%%%%%%%%%%%%%%%%%%% menubar items %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------------------------------------------------------------------
function Tools_Callback(hObject, eventdata, handles)
% hObject    handle to Tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function CopyFigure_Callback(hObject, eventdata, handles)
% hObject    handle to CopyFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handlles    structure with handles and user data (see GUIDATA)
fig = gcf;
print(fig, '-clipboard','-dmeta');


% --------------------------------------------------------------------
function SaveFigure_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename = uiputfile('.fig');
savefig(filename);

% --------------------------------------------------------------------
function Zoom_Callback(hObject, eventdata, handles)
% hObject    handle to Zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%zoom on;
zoom on;

% --------------------------------------------------------------------
function DataCursor_Callback(hObject, eventdata, handles)
% hObject    handle to DataCursor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

datacursormode toggle;
  
% --------------------------------------------------------------------
function Colormap_Callback(hObject, eventdata, handles)
% hObject    handle to Colormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function original_Callback(hObject, eventdata, handles)
% hObject    handle to original (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
colormap(original);

% --------------------------------------------------------------------
function jet_Callback(hObject, eventdata, handles)
% hObject    handle to jet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
colormap(jet);

% --------------------------------------------------------------------
function hot_Callback(hObject, eventdata, handles)
% hObject    handle to hot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
colormap(hot);

% --------------------------------------------------------------------
function gray_Callback(hObject, eventdata, handles)
% hObject    handle to gray (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
colormap(gray);

    
%%%%%%%%%%%%%%%%%%%%%%%% End of menu items %%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on selection change in PlotOut.
function PlotOut_Callback(hObject, eventdata, handles)
% hObject    handle to PlotOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PlotOut contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PlotOut

AnalyzedData = handles.AnalyzedData;
TimeArray_select = AnalyzedData(2).data;
max_locs = AnalyzedData(4).data;
spike_locs = AnalyzedData(7).data;

axes(handles.axes3);
val = get(hObject,'Value');
switch val

    case 2  %Amplitude
        axes(handles.axes3);
        cla reset;
        yyaxis left
        plot(TimeArray_select(max_locs),AnalyzedData(8).data,'*-');
        ylabel('AmplitudeSP (mV)');
        yyaxis right
        plot(TimeArray_select(max_locs),AnalyzedData(9).data,'*-');
        ylabel('AmplitudePP (mV)');
        xlabel('Time (s)');
        set(handles.DisplayText,'String',{['AvgAmp start-to-Peak (mV):  ' num2str(mean(AnalyzedData(8).data))];...
            ['  AvgAmp peak-to-peak (mV):  ' num2str(mean(AnalyzedData(9).data))]}); 
        
    case 3  %TimeInterval
        axes(handles.axes3);
        cla reset;
        plot(TimeArray_select(max_locs),AnalyzedData(10).data,'*-');
        ylabel('CycleLength (s)');  
        xlabel('Time (s)');
        set(handles.DisplayText,'String',{['Avg cycle time (s):  ' num2str(mean(AnalyzedData(10).data))];...
            ['  Avg max-to-min ratio:  ' num2str(mean(AnalyzedData(25).data))]});     
        
    case 4 % t_rise (20% to 80% APD)
        axes(handles.axes3);
        cla reset;
        plot(TimeArray_select(spike_locs(1:end-1)),AnalyzedData(17).data,'*-');
        ylim([0 0.05]);
        ylabel('Duration APD20-80 (s)');
        xlabel('Depol Time (s)');
        set(handles.DisplayText,'String',{['Avg spiking velocity (mv/s):  ' num2str(mean(AnalyzedData(16).data))];...
            ['  Avg Spiking duration (s):  ' num2str(mean(AnalyzedData(17).data))]}); 
% %         axes(handles.axes3);
% %         cla reset;
% %         yyaxis left;
% %         plot(TimeArray_select(spike_locs(1:end-1)),AnalyzedData(16).data,'*-');
% %         ylabel('Spiking velocity (mV/s)');
% %         yyaxis right;
% %         plot(TimeArray_select(spike_locs(1:end-1)),AnalyzedData(17).data,'*-');
% %         ylabel('Duration APD10-90 (s)');
% %         xlabel('Time (s)');
% %         set(handles.DisplayText,'String',{['Avg spiking velocity (mv/s):  ' num2str(mean(AnalyzedData(16).data))];...
% %             ['  Avg Spiking duration (s):  ' num2str(mean(AnalyzedData(17).data))]});     
        
    case 5 %APD 90, 50, 10
        axes(handles.axes3);
        cla reset;
        hold off;
        plot(TimeArray_select(max_locs), AnalyzedData(13).data, '*-');
        hold on;
        plot(TimeArray_select(max_locs), AnalyzedData(14).data, '*-');
        plot(TimeArray_select(max_locs), AnalyzedData(15).data, '*-');
        hold off;
        ylabel('APDs (s)');
        xlabel('Time (s)');
        set(handles.DisplayText,'String',{['APD90 (s): ',num2str(mean(AnalyzedData(13).data))];...
            ['APD50 (s): ', num2str(mean(AnalyzedData(14).data))];...
            ['APD10 (s): ',num2str(mean(AnalyzedData(15).data))]});
        
     case 6 %Triangulation
        axes(handles.axes3);
        cla reset;

        %%%
        plot(TimeArray_select(max_locs), AnalyzedData(24).data, '*-');
        ylabel('Triangulation/s');
        xlabel('Time (s)');
        set(handles.DisplayText,'String',{['Mean triangulation value',num2str(mean(AnalyzedData(24).data))]});

     case 7 %Overpotential
        axes(handles.axes3);
        cla reset;
        AmpPP = AnalyzedData(9).data;
        AmpSP = AnalyzedData(8).data;
        overpot = (AmpPP-AmpSP)./AmpPP;
        plot(TimeArray_select(max_locs), overpot, '*-');
        ylabel('Overpot/Amp ratio');
        xlabel('Time (s)');
        set(handles.DisplayText,'String',{['Mean overpotential ratio',num2str(mean(overpot))]});
        
% %      case 8 %Triangulation (ratio method)
% %         axes(handles.axes3);
% %         cla reset;
% %         %%% this is the old ratio method
% %         APD90 = AnalyzedData(13).data;
% %         APD50 = AnalyzedData(14).data;
% %         APD10 = AnalyzedData(15).data;
% %         tri = (APD90-APD50)./(APD50-APD10);
% %         plot(TimeArray_select(max_locs), tri, '*-');
% %         ylabel('Triangulation ratio');
% %         xlabel('Time (s)');
% %         set(handles.DisplayText,'String',{['Mean triangulation ratio',num2str(mean(tri))]});
        
        case 8  %Replolarization Time to Cycle Length ratio
        axes(handles.axes3);
        cla reset;
        plot(TimeArray_select(max_locs),AnalyzedData(25).data,'*-');
        ylabel('RT/CL ratio');  
        xlabel('Time (s)');
        set(handles.DisplayText,'String',{['RT/CL ratio:  ' num2str(mean(AnalyzedData(25).data))];...
            ['  Avg max-to-min ratio:  ' num2str(mean(AnalyzedData(25).data))]});        
end

xlabel('time (s)')
set( findobj(gca,'typehFig','line'), 'LineWidth', 1);
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 8)
handles.val = val;
guidata(hObject,handles);
% --- Executes during object creation, after setting all properties.
function PlotOut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlotOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveData.
function saveData_Callback(hObject, eventdata, handles)
% hObject    handle to saveData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% saves all data in .mat struct
savefilename = strcat('Analyzed_', handles.filename(1:end - 3), '_', string(handles.channelLabel), '.mat');
[filename,pathname]=uiputfile(savefilename);
save_file = fullfile(pathname, filename);
AnalyzedData=handles.AnalyzedData;
save(save_file,'AnalyzedData');
guidata(hObject,handles);

% --- Executes on button press in openAnalyzed.
function openAnalyzed_Callback(hObject, eventdata, handles)
% hObject    handle to openAnalyzed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Saves currently displayed plotted worked up data only, as a .txt csv.

    fclose all;
    [datafile,pathname]=uigetfile('Analyzed_*.mat','Choose an analyzed data file'); 
    if datafile==0, return; end
    cd(pathname);
    filename=datafile;
    file = fullfile( pathname , filename);
    clear AnalyzedData;
    load(file);
    handles.AnalyzedData=AnalyzedData;
    guidata(hObject,handles);



% --- Executes on button press in baseline.
function baseline_Callback(hObject, eventdata, handles)
% hObject    handle to baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

DataArray_select = handles.DataArray_select;
TimeArray_select = handles.TimeArray_select;
sampling_time = handles.sampling_time;
sampling_frequency = handles.sampling_frequency;
lambda = 2*10^13*(sampling_frequency/5000)^3;

numiter = 3; %%numer of iteraters
% als settings: p = 0.015
Data_select_baseline = als(DataArray_select, lambda, 0.015, numiter);

axes(handles.axes1);
cla reset;
hold on;
plot(TimeArray_select,DataArray_select,'b');
plot(TimeArray_select,DataArray_select-Data_select_baseline,'r');
hold off;
set( findobj(gca,'typehFig','line'), 'LineWidth', 1);
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 8)
xlabel('time (s)')
ylabel('V/(mV)')

handles.Data_select_baseline = Data_select_baseline;
guidata(hObject,handles);

function z = als(y, lambda, p, numiter)
% Estimate baseline with asymmetric least squares
if isrow(y)
    y=y';
end
m = length(y);
D = diff(speye(m), 2);w = ones(m, 1);
for i = 1:numiter
    W = spdiags(w, 0, m, m);
    C = chol(W + lambda * (D' * D));
    z = C \ (C' \ (w .* y));
    w = p * (y > z) + (1 - p) * (y < z);
end



function SlidingWindow_Callback(hObject, eventdata, handles)
% hObject    handle to SlidingWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SlidingWindow as text
%        str2double(get(hObject,'String')) returns contents of SlidingWindow as a double


% --- Executes during object creation, after setting all properties.
function SlidingWindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SlidingWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function StdThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to StdThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StdThreshold as text
%        str2double(get(hObject,'String')) returns contents of StdThreshold as a double


% --- Executes during object creation, after setting all properties.
function StdThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StdThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SpikeIntervalThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to SpikeIntervalThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SpikeIntervalThreshold as text
%        str2double(get(hObject,'String')) returns contents of SpikeIntervalThreshold as a double


% --- Executes during object creation, after setting all properties.
function SpikeIntervalThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SpikeIntervalThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_files_Callback(hObject, eventdata, handles)
% hObject    handle to num_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_files as text
%        str2double(get(hObject,'String')) returns contents of num_files as a double


% --- Executes during object creation, after setting all properties.
function num_files_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in analyze_all_channels.
function analyze_all_channels_Callback(hObject, eventdata, handles)
% hObject    handle to analyze_all_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    DataArrayAll = handles.DataArrayAllChannels;
    channelInfo = handles.channelInfo;
    TimeArray = handles.TimeArray;
    sampling_time = handles.sampling_time;
    sampling_frequency = handles.sampling_frequency;
    filename = handles.filename;
    save_path = append(handles.datapath, filename(1:strfind(filename,'.')-1)); %save into subfolder with filename name
    mkdir(save_path);
    %grab all times to analyze
    stim_time = str2num(get(handles.stim_start, 'String'));
    t_interval = str2num(get(handles.time_skip, 'String'));
    end_time = str2num(get(handles.end_time, 'String'));
    start_times = (stim_time:t_interval:end_time);
    
    % determine window size for analysis; default set to 10s
    window_size = str2num(get(handles.window_size, 'String'));
    % check if data should start with offset time
    offset_start = get(handles.offset_start, 'Value');
   
    if offset_start
        start_times(1) = stim_time + 60; %timepoint 1 at 1 min instead of 0
    end
    
    for tm =1:numel(start_times)
        set(handles.DisplayText,'String',{['going through timepoint: ',num2str(tm),...
            ' of ', num2str(numel(start_times))]});
        % find the start and end of the period to analyze
        startIndex = max(1,fix(start_times(tm)/sampling_time));
        endIndex = min(length(TimeArray), fix((start_times(tm)+window_size)/sampling_time));
        TimeArray_select = TimeArray(startIndex:endIndex);
        labels = channelInfo.Label;
        
        APD90_avg = [];
        APD50_avg = [];
        APD10_avg = [];
        t_rise_avg = [];
        SP_Amp_avg = [];
        PP_Amp_avg = [];
        frequency = [];
        triangle = [];
        overpotential = [];
        max_min_time = [];

        fig = figure(1);
        clf(fig);
        num_good = 0;
        hold on
        for i = 1:60 % iterate through number of channels
            DataArray = double(DataArrayAll(:, i))*10^-3 * 2.3585; %%2.3585 factor determined by the company
            channelLabel = channelInfo.Label{i};
            subplotRow = str2num(channelLabel(2));
            subplotCol = str2num(channelLabel(1));
            subplotPos = (subplotRow - 1) * 8 + subplotCol;
            DataArray_select = DataArray(startIndex:endIndex);
            
            %find stdev
            %%average over 1ms data
            AvgRange=fix(1*10^(-3)/sampling_time); 
            FilteredData=smooth(DataArray_select,AvgRange);
            stdData = movstd(FilteredData,AvgRange*4);
            spikeIntervalThreshold = 0.4; %threshold of 0.4
            %std threshold set at min of max(stdData/2) or 0.05
            maxStd = max(stdData);
            stdThreshold = min(maxStd/2, 0.05);
            sortedStd=sort(stdData,'ascend');
            stdSize = numel(sortedStd);
            %find peaks with rough algorithm first, then filter down later
            [pks_temp,spike_locs_temp]=findpeaks(stdData,'MinPeakProminence',stdThreshold,...
                'MinPeakDistance',spikeIntervalThreshold*sampling_frequency); 
            spike_locs = spike_locs_temp(pks_temp>(maxStd/3));
            pks = pks_temp(pks_temp>(maxStd/3));
            spike_locs = spike_locs(2:numel(spike_locs)); %%this line is to take care of corner cases when the selection starts at spike location of the AP.
            pks = pks(2:numel(pks));
            %spiking frequency
            
            %perform rest of analysis
            stdmedian = median(stdData);
            %stdofstd = std(sortedStd(1:fix(stdSize*0.6)));
            clear max_locs max_vals min_locs min_vals start_locs start_vals;
            clear SP_Amp PP_Amp APD90 APD50 APD10 t_rise triangulation max_min;
            if numel(spike_locs) > 2 && numel(spike_locs) < window_size*3 %shouldn't be beating > 3hz, and need at least 3 for frequency
                for j=1:numel(spike_locs)-1
                    %%finding the maximum and miniumum points position for each AP
                    periodData = FilteredData(spike_locs(j):spike_locs(j+1));
                    [value, locs] = max(periodData);
                    max_locs(j) = spike_locs(j)+locs-1;
                    max_vals(j) = value;

                    %%finding the start point position for each AP
                    %%as of 2022-04-04, use stdmedian and median baseline
                    %%value instead.
                    start_pos=spike_locs(j); 
                    while start_pos>1 && stdData(start_pos)>stdmedian %%replacing stdmedian+3*stdofstd
                        start_pos=start_pos-1;
                    end
                    
                    start_locs(j)=start_pos;
                    periodPoints = numel(periodData); %% number of data points in a period cycle
                    if start_pos - fix(periodPoints*0.05) < 1
                        start_vals(j)= periodData(start_pos);
                    else
                        start_vals(j) = median(FilteredData(start_pos-fix(periodPoints*0.05):start_pos)); 
                    end
                    
                    %%Added on 2022-4-05. calculate the end of an action potential, approximately the
                    %%minumum. The variable is called min_locs and min_vals. They work for traces that do not have a
                    %%clear minimum.  
                    Half_Amp = start_vals(j)+(max_vals(j)-start_vals(j))*0.1; %%APD90
                    Half_Amp_points = numel(find(periodData>Half_Amp));  %%number of points with values greater than Half_amp
                    if Half_Amp_points > 1
                        Range50ms = 50*10^(-3)/sampling_time;
                        stdData50ms = movstd(periodData,Range50ms);
                        range_start = Half_Amp_points;
                        range_end = min(fix(Half_Amp_points*3), numel(periodData));
                        med_stdData50ms = median(stdData50ms(range_start:range_end));
                        APend_pos = Half_Amp_points;
                        %%APend_pos means the end of the repolarization phase
                        while APend_pos>=1 && stdData50ms(APend_pos)> med_stdData50ms 
                            APend_pos=APend_pos+1;
                        end
                        min_loc_1 = spike_locs(j)+APend_pos-1;
                        min_val_1=FilteredData(min_loc_1);
                        [value, locs] = min(periodData);
                        min_loc_2 = spike_locs(j)+locs-1;
                        min_val_2 = value;
                        if min_loc_2 < min_loc_1
                            min_vals(j) = min_val_2;
                            min_locs(j) = min_loc_2;
                        else
                            min_vals(j) = min_val_1;
                            min_locs(j) = min_loc_1;
                        end
                    else
                        min_locs(j) = spike_locs(j)+1; %%strange results, purposely set min to max
                        min_vals(j)=FilteredData(min_locs(j));
                    end

                    SP_Amp(j) = max_vals(j)-start_vals(j); %%amplitude start to max
                    PP_Amp(j) = max_vals(j)-min_vals(j);   %%amplitude peak to peak
                    periodData_APD = FilteredData(start_locs(j):min_locs(j)); %%redefine the period data
                    APD90(j) = numel(find(periodData_APD > (start_vals(j)+SP_Amp(j)*0.1)))*sampling_time;
                    APD50(j) = numel(find(periodData_APD > (start_vals(j)+SP_Amp(j)*0.5)))*sampling_time;
                    APD10(j) = numel(find(periodData_APD > (start_vals(j)+SP_Amp(j)*0.9)))*sampling_time;

                    %t_rise is the time from APD80 to APD20
                    periodData_trise = FilteredData(start_locs(j):max_locs(j)); %%redefine the period data
                    t_rise(j) = numel(find(periodData_trise > (start_vals(j)+SP_Amp(j)*0.2) &...
                                    periodData_trise < (start_vals(j)+SP_Amp(j)*0.8)))*sampling_time;
                    %%calculate triangulation of data (APD90-APD30 in the phase after
                    %%the peak) and the duration of the entire phase
                    periodData_tri = FilteredData(max_locs(j):min_locs(j)); %%redefine the period data
                    max_min(j) = numel(periodData_tri)*sampling_time;
                    triangulation(j) = numel(find(periodData_tri > (start_vals(j)+SP_Amp(j)*0.1) &...
                                            periodData_tri < (start_vals(j)+SP_Amp(j)*0.7)))*sampling_time;
                end
            else %data is bad
                APD50 = [];
                APD90 = [];
                APD10 = [];
                t_rise = [];
                SP_Amp = [];
                PP_Amp = [];
                max_min = [];
            end
                
            %determine if data is good by consistency of APD and t-rise
            apd_value = mean(APD90)>0.1; %good APD value should be >0.1
            %if values fluctuate a lot it's not good
            apd_fluc = std(APD50)/mean(APD50) < 0.2;
            max_min_fluc = std(max_min)/mean(max_min) < 0.2;
            %amp should not fluctuate too much
            ampSP_fluc = std(SP_Amp)/mean(SP_Amp) < 0.2;
            ampPP_fluc = std(PP_Amp)/mean(PP_Amp) < 0.2;
            %trise should not fluctuate too much
            trise_fluc = std(t_rise)/mean(t_rise) < 0.5;
            data_good = apd_value && apd_fluc && max_min_fluc && ampSP_fluc && trise_fluc && ampPP_fluc;
            if subplotPos == 7 || subplotPos == 15 || subplotPos == 50 % special case
                data_good = false;
            end
            
            %store relevant data
            if data_good
                num_good = num_good + 1;
                APD90_avg(i) = mean(APD90);
                APD50_avg(i) = mean(APD50);
                APD10_avg(i)= mean(APD10);
                t_rise_avg(i) = mean(t_rise);
                SP_Amp_avg(i) = mean(SP_Amp);
                PP_Amp_avg(i) = mean(PP_Amp);
                triangle(i) = mean(triangulation);
                max_min_time(i) = mean(max_min);
                
                % calculate overpotential from PP/SP ratio               
                overpotential(i) = mean(SP_Amp)/(mean(PP_Amp)-mean(SP_Amp));
                % calculate frequency from spike_locs
                first_pk_time = TimeArray_select(spike_locs(1));
                last_pk_time = TimeArray_select(spike_locs(numel(spike_locs)));
                frequency(i) = (numel(spike_locs)-1)/(last_pk_time - first_pk_time); 
                
            else
                APD90_avg(i) = NaN;
                APD50_avg(i) = NaN;
                APD10_avg(i)= NaN;
                t_rise_avg(i) = NaN;
                SP_Amp_avg(i) = NaN;
                PP_Amp_avg(i) = NaN;
                frequency(i) = NaN;
                triangle(i) = NaN;
                overpotential(i) = NaN;
                max_min_time(i) = NaN;
            end
                
            
            %plot timepoint
            if subplotPos == 7 || subplotPos == 15 || subplotPos == 50 %%% this channel is always a bad channel and set to 0. 
                subplot(8, 8, subplotPos)
                plot(TimeArray_select,0*DataArray_select);
                set( findobj(gca,'typehFig','line'), 'LineWidth', 1);
                xticks([]);
                ylim([-1, 1]);
                t = title(append(string(channelLabel), ' ', string(data_good),'    '));

            else
                subplot(8, 8, subplotPos)
                plot(TimeArray_select,DataArray_select);
                hold on;
                plot(TimeArray_select(spike_locs), DataArray_select(spike_locs), 'r*');
                if data_good
                    plot(TimeArray_select(start_locs), DataArray_select(start_locs), 'b*');
                    plot(TimeArray_select(min_locs), DataArray_select(min_locs), 'c*');
                end
                set( findobj(gca,'typehFig','line'), 'LineWidth', 1);
                xticks([]);
                ylim('auto');
                temp=[min(DataArray_select), max(DataArray_select)];
                extended_yrange = 0.2*(temp(2)-temp(1));
                extended_yrange = max(extended_yrange, 0.25);
                ylim([temp(1)-extended_yrange temp(2)+extended_yrange]);
                t = title(append(string(channelLabel), ' ', string(data_good),'    '));
            end
            t.Units = 'Normalize'; 
            t.Position(1) = 0;
            t.Position(2) = 1;
            t.HorizontalAlignment = 'left'; 
        end
        
        sgtitle({strrep(append(filename, ' open=', string(num_good), ' ', string(start_times(tm)), ' sec - ',...
                string(start_times(tm)+window_size), ' sec'), '_', '\_'), ''});
        % store data at this timepoint
        APD90_avg = APD90_avg';
        APD50_avg = APD50_avg';
        APD10_avg = APD10_avg';
        t_rise_avg = t_rise_avg';
        SP_Amp_avg = SP_Amp_avg';
        PP_Amp_avg = PP_Amp_avg';
        frequency = frequency';
        triangle = triangle';
        overpotential = overpotential';
        max_min_time = max_min_time';
        table_data = table(labels,APD90_avg,APD50_avg,APD10_avg,t_rise_avg,SP_Amp_avg,PP_Amp_avg, ...
                           frequency,triangle,overpotential,max_min_time);
        
        % check for autosave
        autosave = get(handles.autosave_check, 'Value');
        if autosave
            fname = append('data_t', string(round((tm-1)*t_interval/60)), 'min.txt');
            writetable(table_data, fullfile(save_path, fname));
            %imgname = append('image', string(round((tm-1)*t_interval/60)), 'min.fig');
            %saveas(fig, fullfile(save_path, imgname));
            imgname_2 = append('image', string(round((tm-1)*t_interval/60)), 'min.png');
            saveas(fig, fullfile(save_path, imgname_2));
        else
            [filename,pathname]=uiputfile('output.txt');
            imgname_2 = append(filename(1:strfind(filename,'.')-1), '.png');
            writetable(table_data, fullfile(pathname, filename));
            saveas(fig, fullfile(pathname, imgname_2));

        end
    end

   guidata(hObject,handles);



function stim_start_Callback(hObject, eventdata, handles)
% hObject    handle to stim_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stim_start as text
%        str2double(get(hObject,'String')) returns contents of stim_start as a double


% --- Executes during object creation, after setting all properties.
function stim_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function time_skip_Callback(hObject, eventdata, handles)
% hObject    handle to time_skip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_skip as text
%        str2double(get(hObject,'String')) returns contents of time_skip as a double


% --- Executes during object creation, after setting all properties.
function time_skip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_skip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function end_time_Callback(hObject, eventdata, handles)
% hObject    handle to end_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of end_time as text
%        str2double(get(hObject,'String')) returns contents of end_time as a double


% --- Executes during object creation, after setting all properties.
function end_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to end_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in offset_start.
function offset_start_Callback(hObject, eventdata, handles)
% hObject    handle to offset_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of offset_start



function window_size_Callback(hObject, eventdata, handles)
% hObject    handle to window_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of window_size as text
%        str2double(get(hObject,'String')) returns contents of window_size as a double


% --- Executes during object creation, after setting all properties.
function window_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to window_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in autosave_check.
function autosave_check_Callback(hObject, eventdata, handles)
% hObject    handle to autosave_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autosave_check
