% This function is used to read the raw data and concanate data from
% different subjects into one variable 
clear;clc;addpath('..\Function');
%%%%%%%%%%%%%%%%%%%%%%%%% Setting and parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
passBand = [7 70]; % Refer to Chen et al (2015) paper
stopBand = [6 71];
Rp = 0.5; Rs = 100;
bandPassFilter = true; % Apply the bandpass filter or not
CAR = false; % Apply the common average reference or not
time = 6; % Keep the whole signal including the pre-stimulus
%%%%%%%%%%%%%%%%% Read value of stimulus frequency and phase %%%%%%%%%%%%%%
addpath('C:\SSVEP\Tsinghua database 2016');
temp = load('Phase.mat');
stimuFreq = temp.freqs;
[stimuFreq,nidx] = sort(stimuFreq);
fsample = 250;
numSubject = 35;
trialLength = 6;
freqLength = 40;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create data name that are going to load in workspace
dataName = cell(numSubject,1);
for subject = 1:numSubject
    dataName{subject} = ['s' num2str(subject) '.mat'];
end
fprintf('Stimulus frequency: \n');
disp(stimuFreq);
fprintf('Sampling frequency: \n');
disp(fsample);
numSubject = length(dataName);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter design (cheby1chev type I) for extract SSVEP/Alpha/Beta related signal
if bandPassFilter
    ssvep_filter = genBPFilter(passBand,stopBand,Rp,Rs,fsample);
end
alpha_filter = genBPFilter([8 13],[7.5 13.5],Rp,Rs,fsample);
lowbeta_filter = genBPFilter([13 16],[12.5 16.5],Rp,Rs,fsample);
beta_filter = genBPFilter([16 20],[15.5 20.5],Rp,Rs,fsample);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove the redundant time segment and channels
% delay = 0.5*fsample; %0.5 record for pre-stimulus onset
delay = 0; % (23/03/2017): Since the pre-stimulus signal may be useful, keep it
sChannel = [48 54:58 61:63]; %(Pz, PO5, PO3, POz, PO4, PO6, O1, Oz, and O2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intialize the variable
preL = round(0.7*fsample); % Segement for pre-stimulus
allData = zeros(trialLength,freqLength,time*fsample,length(sChannel),numSubject);
prestimulus_alpha = zeros(trialLength,freqLength,preL,length(sChannel),numSubject);
prestimulus_lowbeta = zeros(size(prestimulus_alpha));
prestimulus_beta = zeros(size(prestimulus_alpha));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data preprocessing (Read -> reshape -> concatenate)
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
                data(trial,freq,:,:) = filtfilt(ssvep_filter(1,:),ssvep_filter(2,:),bsxfun(@minus,x,mean(x,2)));
            elseif bandPassFilter
                data(trial,freq,:,:) = filtfilt(ssvep_filter(1,:),ssvep_filter(2,:),x);
            elseif CAR
                data(trial,freq,:,:) = bsxfun(@minus,x,mean(x,2));
            else
                data(trial,freq,:,:) = x;
            end
            % if CAR is needed, may need to rearrange the below code
            prestimulus_alpha(trial,freq,:,:,targetSubject) = filtfilt(alpha_filter(1,:),alpha_filter(2,:),x(1:preL,:));
            prestimulus_lowbeta(trial,freq,:,:,targetSubject) = filtfilt(lowbeta_filter(1,:),lowbeta_filter(2,:),x(1:preL,:));
            prestimulus_beta(trial,freq,:,:,targetSubject) = filtfilt(beta_filter(1,:),beta_filter(2,:),x(1:preL,:));
        end
    end
    allData(:,:,:,:,targetSubject) = data;
end
dataSize = size(allData);
rmpath('C:\SSVEP\Tsinghua database 2016');
str = sprintf('benchmark_ssvep_prefilter%d_CAR%d',bandPassFilter,CAR);
save(['..\ProcessedData\' str '.mat'],'allData','stimuFreq','fsample','dataSize');
str = sprintf('benchmark_prestimulus_rhythm');
save(['..\ProcessedData\' str '.mat'],'prestimulus_alpha','prestimulus_lowbeta',...
    'prestimulus_beta','stimuFreq','fsample','preL');