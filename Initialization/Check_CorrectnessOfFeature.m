clear;clc;addpath('..\Function');
%% Testing parameter
nameDataset = 'Benchmark';
method = 'ECCA4';
%% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load data to workspace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% allData: SSVEP data in cell format
% allData{j} contains SSVEP data from j^th subject
prefilter = true; % Use the data that have been filtered (for the whole 5s signal)
name = sprintf('benchmark_ssvep_prefilter%d_CAR0.mat',prefilter);
load(['..\ProcessedData\' name]);
% warning('off','stats:canoncorr:NotFullRank');
startTime = 0.14+0.5; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Basic info of data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trialLength = dataSize(1); %Number of recorded EEG response for each stimulus frequency
freqLength = dataSize(2); %Number of visual stimulus 
sampleLength = dataSize(3); %Number of time points in each record
channelLength = dataSize(4); %Number of channels used in each experiment
numSubject = dataSize(5); %Number of subjects in the data set
trialSeq = 1:trialLength;
subjectSeq = 1:numSubject;
startIdx = round(fsample*startTime);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('..\Feature\ECCA4_Benchmark_Nh5_time100_filter1_cv.mat');

% Check random trial and and random frequency and random source data and
% random target template
sinTemplate = genSinTemplate(stimuFreq,fsample,1,5);
for ch = 1:10
    ranTrial = randsample(1:6,1);
    ranFreq = randsample(1:40,1);
    ranSource = randsample(1:35,1);
    ranTarget = randsample(1:35,1);
    
    cvTrial = trialSeq ~= ranTrial;
    
    Xnew = squeeze(allData(ranTrial,ranFreq,startIdx+1:startIdx+floor(fsample*1),:,ranTarget));
    Xtemplate = squeeze(mean(allData(cvTrial,:,startIdx+1:startIdx+floor(fsample*1),:,ranSource)));
    [~,pvec] = ccaExtend(Xnew,Xtemplate,sinTemplate,'Combination4');
    
    result = squeeze(trainFeatureSet(ranTarget,(ranTrial-1)*freqLength+ranFreq,ranSource,:));
    if all((pvec(:,1) - result) < 1e-9)
        fprintf('Reuslt no problem, go ahead\n');
    else
        fprintf('Oh no!\n');
    end
end

fprintf('\n\n\n');

for ch = 1:10
    ranTrial = randsample(1:6,1);
    ranFreq = randsample(1:40,1);
    ranSource = randsample(1:35,1);
    ranTarget = randsample(1:35,1);
    while ranSource == ranTarget
        ranTarget = randsample(1:35,1);
    end
    feataureSelector = find(subjectSeq ~= ranTarget);
    cvTrial = trialSeq ~= ranTrial;
    
    Xnew = squeeze(all_ssvep(ranTrial,ranFreq,startIdx+1:startIdx+floor(fsample*1),:,ranTarget));
    Xtemplate = squeeze(mean(all_ssvep(:,:,startIdx+1:startIdx+floor(fsample*1),:,ranSource)));
    [~,pvec] = ccaExtend(Xnew,Xtemplate,sinTemplate,'Combination4');
    
    result = squeeze(targetFeatureSet(ranTarget,(ranTrial-1)*freqLength+ranFreq,find(feataureSelector == ranSource),:));
    if all((pvec(:,1) - result) < 1e-9)
        fprintf('Reuslt no problem, go ahead\n');
    else
        fprintf('Oh no!\n');
    end
end

