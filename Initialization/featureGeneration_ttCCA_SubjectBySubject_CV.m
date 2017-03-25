clear;clc;addpath('..\Function\');
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

% One-hot encoded label for each data
label = zeros(numSubject,trialLength*freqLength,freqLength);
for i = 1:numSubject
    label(i,:,:) = repmat(eye(freqLength),trialLength,1);
end
%% Loop through different time length and Nh
NhSeq = 5;
timeSeq = 0.25:0.25:4;

for Nhidx = 1:length(NhSeq)
    Nh = NhSeq(Nhidx);
    for tidx = 1:length(timeSeq);
        time = timeSeq(tidx);
        % Artifical(sinusoidal template)
        sinTemplate = genSinTemplate(stimuFreq,fsample,time,Nh);
        %% Extract features for all data (including test data)
        % featureSet: [testSubject,size of data, views, features]
        % Each subject provide nF views on the input data + one view from CCA with sinusodial template
        trainFeatureSet = zeros(numSubject,trialLength*freqLength,numSubject*3+1,freqLength);
        wlen = floor(time*fsample);
        % Extract feature through CCA with artifical/subject-specific template
        for testSubject = 1:numSubject
            for cvTestTrial = 1:trialLength
                % Take trial 2 - 6 for training and trial 1 for validation
                cvTrainSeq = find(trialSeq ~= cvTestTrial);
                train_ssvep = allData(cvTrainSeq,:,startIdx+1:startIdx+wlen,:,:);
                train_template = squeeze(mean(train_ssvep));
                valid_ssvep = squeeze(allData(cvTestTrial,:,startIdx+1:startIdx+wlen,:,testSubject));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                tempFeatureGroup1 = zeros(freqLength,numSubject,freqLength);
                tempFeatureGroup2 = zeros(freqLength,numSubject,freqLength);
                tempFeatureGroup3 = zeros(freqLength,numSubject,freqLength);
                tempFeatureGroupCCA = zeros(freqLength,numSubject,freqLength);
                
                parfor ithTemplate = 1:numSubject
                    tempFeature1 = zeros(freqLength,freqLength);
                    tempFeature2 = zeros(freqLength,freqLength);
                    tempFeature3 = zeros(freqLength,freqLength);
                    tempFeatureCCA = zeros(freqLength,freqLength);
                    Xtemplate = squeeze(train_template(:,:,:,ithTemplate));
                    for freq = 1:freqLength
                        Xnew = squeeze(valid_ssvep(freq,:,:));
                        [~,pvec] = ccaExtend(Xnew,Xtemplate,sinTemplate,'Combination4');
                        tempFeature1(freq,:) = pvec(:,1);
                        tempFeature2(freq,:) = pvec(:,3);
                        tempFeature3(freq,:) = pvec(:,4);
                        tempFeatureCCA(freq,:) = pvec(:,2);
                    end
                    tempFeatureGroup1(:,ithTemplate,:) = tempFeature1;
                    tempFeatureGroup2(:,ithTemplate,:) = tempFeature2;
                    tempFeatureGroup3(:,ithTemplate,:) = tempFeature3;
                    tempFeatureGroupCCA(:,ithTemplate,:) = tempFeatureCCA;
                end
                insertIdx = (cvTestTrial-1)*freqLength+1:cvTestTrial*freqLength;
                trainFeatureSet(testSubject,insertIdx,1:numSubject,:) = tempFeatureGroup1;
                trainFeatureSet(testSubject,insertIdx,numSubject+1:numSubject*2,:) = tempFeatureGroup2;
                trainFeatureSet(testSubject,insertIdx,numSubject*2+1:numSubject*3,:) = tempFeatureGroup3;
                trainFeatureSet(testSubject,insertIdx,end,:) = squeeze(tempFeatureGroupCCA(:,1,:));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end              
            fprintf(sprintf('Finish s%d\n',testSubject));
        end
        save(['..\Feature\' sprintf('%s_%s_Nh%d_time%d_filter%d_cv.mat',method,nameDataset,Nh,time*100,1)],'trainFeatureSet','label');
        fprintf(sprintf('Have saved %s_%s_Nh%d_time%d_filter%d_cv.mat\n',method,nameDataset,Nh,time*100,1));
    end
end