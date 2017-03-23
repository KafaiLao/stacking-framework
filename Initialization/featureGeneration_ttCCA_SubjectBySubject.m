clear;clc;addpath('Function');
%% Testing parameter
nameDataset = 'Benchmark';
method = 'ECCA4';
%% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load data to workspace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% allData: SSVEP data in cell format
% allData{j} contains SSVEP data from j^th subject
prefilter = true; % Use the data that have been filtered (for the whole 5s signal)
name = sprintf('benchmark_prefilter%d_CAR0.mat',prefilter);
load(['ProcessedData\' name]);
% warning('off','stats:canoncorr:NotFullRank');
startTime = 0.14; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Basic info of data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trialLength = dataSize(1); %Number of recorded EEG response for each stimulus frequency
freqLength = dataSize(2); %Number of visual stimulus 
sampleLength = dataSize(3); %Number of time points in each record
channelLength = dataSize(4); %Number of channels used in each experiment
numSubject = length(allData); %Number of subjects in the data set
trialSeq = 1:trialLength;
subjectSeq = 1:numSubject;
startIdx = round(fsample*startTime);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the average EEG response template for each subject
numTemplate = numSubject + 1;
allTemplate = cell(numSubject,1);
for subject = 1:numSubject
    allTemplate{subject} = squeeze(mean(allData{subject}));
end
allTemplateMat = cat(4,allTemplate{:});

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
        featureSet = zeros(numSubject,trialLength*freqLength,numTemplate*3+1,freqLength);
        % Extract feature through CCA with artifical/subject-specific template
        for testSubject = 1:numSubject
            ssvep = allData{testSubject};
            sourceSeq = find(subjectSeq ~= testSubject);
            allTemplate{numSubject+1} = squeeze(mean(allTemplateMat(:,:,:,sourceSeq),4));

            tempFeature = zeros(freqLength*trialLength,numTemplate*3+1,freqLength);
            tempFeatureGroup1 = zeros(freqLength*trialLength,numTemplate,freqLength);
            tempFeatureGroup2 = zeros(freqLength*trialLength,numTemplate,freqLength);
            tempFeatureGroup3 = zeros(freqLength*trialLength,numTemplate,freqLength);
            tempFeatureGroupCCA = zeros(freqLength*trialLength,numTemplate,freqLength);
            parfor ithTemplate = 1:numTemplate
                tempFeature1 = zeros(freqLength*trialLength,freqLength);
                tempFeature2 = zeros(freqLength*trialLength,freqLength);
                tempFeature3 = zeros(freqLength*trialLength,freqLength);
                tempFeatureCCA = zeros(freqLength*trialLength,freqLength);
                Xtemplate = allTemplate{ithTemplate};
                Xtemplate = Xtemplate(:,startIdx+1:startIdx+floor(time*fsample),:);
                for trial = 1:trialLength
                    for freq = 1:freqLength
                        Xnew = squeeze(ssvep(trial,freq,startIdx+1:startIdx+floor(time*fsample),:));
                        [~,pvec] = ccaExtend(Xnew,Xtemplate,sinTemplate,'Combination4');
                        tempFeature1((trial-1)*freqLength+freq,:) = pvec(:,1);
                        tempFeature2((trial-1)*freqLength+freq,:) = pvec(:,3);
                        tempFeature3((trial-1)*freqLength+freq,:) = pvec(:,4);
                        tempFeatureCCA((trial-1)*freqLength+freq,:) = pvec(:,2);
                    end
                end
                tempFeatureGroup1(:,ithTemplate,:) = tempFeature1;
                tempFeatureGroup2(:,ithTemplate,:) = tempFeature2;
                tempFeatureGroup3(:,ithTemplate,:) = tempFeature3;
                tempFeatureGroupCCA(:,ithTemplate,:) = tempFeatureCCA;
            end
            featureSet(testSubject,:,1:numTemplate,:) = tempFeatureGroup1;
            featureSet(testSubject,:,numTemplate+1:numTemplate*2,:) = tempFeatureGroup2;
            featureSet(testSubject,:,numTemplate*2+1:numTemplate*3,:) = tempFeatureGroup3;
            featureSet(testSubject,:,end,:) = squeeze(tempFeatureGroupCCA(:,1,:));
            
            fprintf(sprintf('Finish s%d\n',testSubject));
        end
        save(['Feature\' sprintf('%s_%s_Nh%d_time%d_filter%d.mat',method,nameDataset,Nh,time*100,1)],'featureSet','label');
        fprintf(sprintf('Have saved %s_%s_Nh%d_time%d_filter%d.mat\n',method,nameDataset,Nh,time*100,1));
    end
end