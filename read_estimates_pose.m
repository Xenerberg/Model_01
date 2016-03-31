%% Import data from text file.
% Script for importing data from the following text file:
%
%    C:\Users\Iseberg\Documents\MATLAB\Model_01\estimates.txt
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2016/01/21 17:06:28
clear segment_data;
clear data;
clear filtered_segment_data;
%% Initialize variables.
filename = './estimates.txt';
delimiter = ' ';

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5,6,7,8,9,10,11,12]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
    
    
    
    
end


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
estimates = raw;

[r,w] = size(estimates);
row_counter = 0;
for i_Count = 1:7:r
    row_counter = row_counter + 1;
    data{row_counter,1} = estimates{i_Count, 2};%Frame number
    data{row_counter,2} = estimates{i_Count + 1, 1};%Time stamp
    data{row_counter,3} = estimates{i_Count + 2, 1};%estimation flag

    for i_Count_2 = 1:4
        temp = cell2mat(estimates(i_Count + i_Count_2+2,1:12));%Trans TCP to grasp
        data{row_counter,3+2*i_Count_2-1} = [temp(1:3);temp(5:7);temp(9:11)];
        data{row_counter,3+2*i_Count_2} = [temp(4);temp(8);temp(12)];
    end    
    
end    
reject = cell2mat(data(cell2mat(data(:,3))~=0));
%plot data characteristics
b_plot = 0;
if (b_plot == 1)
    figure;
    subplot(2,2,1);

    plot(cell2mat(data(:,2)));
    title('time-stamps w.r.t samples');
    subplot(2,2,2);
    stairs(cell2mat(data(:,3)),'Linewidth',3);
    title('flag state');
    subplot(2,2,3);
    hold all;
    title('position vector');
    vector = reshape(cell2mat(data(:,9)),3,110)';
    stem(cell2mat(data(:,2)),vector(:,1),'LineStyle','-.','MarkerEdgeColor','black','Marker','o','MarkerFaceColor','red');
    stem(cell2mat(data(:,2)),vector(:,2),'LineStyle','-.','MarkerEdgeColor','black','Marker','v','MarkerFaceColor','red');
    stem(cell2mat(data(:,2)),vector(:,3),'LineStyle','-.','MarkerEdgeColor','black','Marker','*','MarkerFaceColor','red');
end
%Check for data continuity
% difference_time = diff(cell2mat(data(:,2)));
% logical_index = difference_time < 0.09 | difference_time > 0.11;
% index_discont = find(logical_index==1);
% 
% last_count = 1;
% for i_Count = 1:length(index_discont)
%    segment_data{i_Count} = {data(last_count:index_discont(i_Count),:)};
%    last_count = index_discont(i_Count)+1;
% end
segment_data{1} = {data(1:end,:)};
%Check for flag and remove -1 elements
for i_Count = 1:length(segment_data)
   filtered_segment_data{i_Count} = segment_data{i_Count}{1}(cell2mat(segment_data{i_Count}{1}(:,3)) ~= -1,:);
end

%At present, use only 2 of the segment
data = filtered_segment_data{1};

% data(reject,4) = {nan(3,3)};
% data(reject,6) = {nan(3,3)};
% data(reject,8) = {nan(3,3)};
% data(reject,10) = {nan(3,3)};
% data(reject,5) = {nan(3,1)};
% data(reject,7) = {nan(3,1)};
% data(reject,9) = {nan(3,1)};
% data(reject,11) = {nan(3,1)};
%Fetch data quality
v_quality = cell2mat(data(:,3));
%Assign quality factor with respect to the time-series object
v_quality(v_quality == -1) = -128;
v_quality(v_quality == 0) = 127;
%Create the time vector
time_vector = cell2mat(data(:,2));
%Convert time vector into an epoch measurement that starts at zero
time_vector = time_vector - (time_vector(1));
% signal_1 = zeros(4);
% signal_2 = zeros(4);
counter = 1;
row_counter = length(data);
for i_Count = 4:2:11   
    trans = reshape(cell2mat(data(:,i_Count))',[3,3,row_counter]);
    signal(counter) = {timeseries(trans,time_vector,v_quality,'Name',strcat('Signal_Transf_',num2str(floor((i_Count-3)/2))))};
    signal{counter}.TimeInfo.UserData = 'Sampling time for camera';
    trans = reshape(cell2mat(data(:,i_Count+1))',[3,1,row_counter]);
    signal(counter+1) = {timeseries(trans,time_vector,v_quality,'Name',strcat('Signal_vector_',num2str(floor((i_Count-3)/2))))};
    signal{counter}.TimeInfo.UserData = 'Sampling time for camera';
    counter = counter + 2;   
end
ts_coll = tscollection(signal,'name','Signal_Collection');
filtered_data = data(reject,:);

%% Clear temporary variable
clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;