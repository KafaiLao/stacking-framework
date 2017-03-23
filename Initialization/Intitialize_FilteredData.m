% This function is used to read the raw data and concanate data from
% different subjects into one variable 
clear;clc;addpath('Function');
%%%%%%%%%%%%%%%%%%%%%%% Setting and parameters %%%%%%%%%%%%%%%%%%%%%%%
passBand = [7 70]; % Refer to Chen et al (2015) paper
bandPassFilter = true; % Apply the bandpass filter or not
CAR = false; % Apply the common average reference or not
time = 5.5;
%%%%%%%%%%%%%% Read value of stimulus frequency and phase %%%%%%%%%%%%%%
addpath('C:\SSVEP\Tsinghua database 2016');
temp = load('Phase.mat');
stimuFreq = temp.freqs;
[stimuFreq,nidx] = sort(stimuFreq);
fsample = 250;
numSubject = 35;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create data name that are going to read
dataName = cell(numSubject,1);
for subject = 1:numSubject
    dataName{subject} = ['s' num2str(subject) '.mat'];
end
fprintf('Stimulus frequency: \n');
disp(stimuFreq);
fprintf('Sampling frequency: \n');
disp(fsample);
numSubject = length(dataName);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter design
if bandPassFilter
    % cheby1chev type I
    chebyord(
    bpFilt = designfilt('bandpassiir','DesignMethod','cheby1','FilterOrder',6, ...
        'PassbandFrequency1',passBand(1),'PassbandFrequency2',passBand(2), ...
        'PassbandRipple',0.5,'SampleRate',fsample);
    % Butterworth
%     bpFilt = designfilt('bandpassiir','FilterOrder',10, ...
%         'HalfPowerFrequency1',passBand(1),'HalfPowerFrequency2',passBand(2), ...
%         'SampleRate',fsample);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove the redundant time segment and channels
delay = 0.5*fsample; %0.5 record for pre-stimulus onset
sChannel = [48 54:58 61:63]; %(Pz, PO5, PO3, POz, PO4, PO6, O1, Oz, and O2)
% Data preprocessing (Read -> reshape -> concatenate)
allData = cell(numSubject,1);

for targetSubject = 1:numSubject
    temp = load(strtrim(dataName{targetSubject}));
    data = temp.data;
    % [Channels, Samples, Freqencies, Trials]
    data = data(sChannel,delay+1:delay+fsample*time,:,:);
    % [Trials, Freqencies, Samples, Channels]
    data = permute(data,[4 3 2 1]);
    ndata = data;
    
    for trial = 1:size(data,1)
        for freq = 1:size(data,2)
            x = squeeze(ndata(trial,nidx(freq),:,:));
            if bandPassFilter && CAR
                data(trial,freq,:,:) = filtfilt(bpFilt,bsxfun(@minus,x,mean(x,2)));
            elseif bandPassFilter
                data(trial,freq,:,:) = filtfilt(bpFilt,x);
            elseif CAR
                data(trial,freq,:,:) = bsxfun(@minus,x,mean(x,2));
            else
                data(trial,freq,:,:) = x;
            end
        end
    end
    allData{targetSubject} = data;
end
dataSize = size(data);
rmpath('C:\SSVEP\Tsinghua database 2016');
str = sprintf('benchmark_prefilter%d_CAR%d',bandPassFilter,CAR);
save(['ProcessedData\' str '.mat'],'allData','stimuFreq','fsample','dataSize');