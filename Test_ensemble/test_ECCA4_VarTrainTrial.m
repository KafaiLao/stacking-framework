clearvars -except list method;clc;addpath('Function');
%% Testing parameter
method = 'ECCA4';
nameDataset = 'Benchmark';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
% Load data to workspace
% allData: SSVEP data in cell format
% allData{j} contains SSVEP data from j^th subject
prefilter = true; % Use the data that have been filtered (for the whole 5s signal)
name = sprintf('benchmark_prefilter%d_CAR0.mat',prefilter);
load(['ProcessedData\' name]);
fprintf('Stimulus frequency:\n');
disp(stimuFreq);
%%%%%%%%%%%%%%%%%% Specific starting time (i.e., remove some initial data points %%%%%%%%%%%%%%%%%%%%%
startTime = 0.14; 
startIdx = round(fsample*startTime);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic info of data
trialLength = dataSize(1); %Number of recorded EEG response for each stimulus frequency
freqLength = dataSize(2); %Number of visual stimulus 
sampleLength = dataSize(3); %Number of time points in each record
channelLength = dataSize(4); %Number of channels used in each experiment
numSubject = length(allData); %Number of subjects in the data set
trialSeq = 1:trialLength;
subjectSeq = 1:numSubject;
%% Recognition
timeSeq = 0.25:0.25:4;
NhSeq = 5;
repetition = 30;
for numTrainTrial = 1:trialLength - 2
    rec_confusion = zeros(length(NhSeq),length(timeSeq),numSubject,freqLength,freqLength);
    for Nhidx = 1:length(NhSeq)
        Nh = NhSeq(Nhidx);
        acc = zeros(numSubject+2,length(timeSeq));
        for tidx = 1:length(timeSeq)
            tic;
            time = timeSeq(tidx);
            % Artifical(sinusoidal template)
            sinTemplate = genSinTemplate(stimuFreq,fsample,time,Nh);
            rec = zeros(numSubject,1);
            parfor targetSubject = 1:numSubject
                ssvep = allData{targetSubject};
                rec_subjectSpecific = zeros(freqLength);
                for repeatCount = 1:repetition
                    for trial = 1:trialLength
                        trainTestSeq = circshift(trialSeq,[0 -trial]);
                        trainIdx = randsample(trainTestSeq(1:end-1),numTrainTrial);
                        if numTrainTrial ~= 1
                            XTemplate = squeeze(mean(ssvep(trainIdx,:,startIdx+1:startIdx+floor(time*fsample),:)));
                        else
                            XTemplate = squeeze(ssvep(trainIdx,:,startIdx+1:startIdx+floor(time*fsample),:));
                        end
                        for freq = 1:freqLength
                            Xnew = squeeze(ssvep(trial,freq,startIdx+1:startIdx+floor(time*fsample),:));
                            score = ccaExtend(Xnew,XTemplate,sinTemplate,'Combination4');
                            [~,maxLoc] = max(score);
                            if maxLoc == freq, rec(targetSubject) = rec(targetSubject) + 1; end;
                            rec_subjectSpecific(freq,maxLoc) = rec_subjectSpecific(freq,maxLoc) + 1;
                        end
                    end
                end
                rec_confusion(Nhidx,tidx,targetSubject,:,:) = rec_subjectSpecific;
                fprintf('Accuracy of S%d: %.2f\n',targetSubject,rec(targetSubject)*100/(trialLength*freqLength)/repetition);
            end
            acc(1:numSubject,tidx) = rec/(trialLength*freqLength)*100/repetition;
            acc(numSubject+1,tidx) = mean(acc(1:numSubject,tidx));
            acc(numSubject+2,tidx) = std(acc(1:numSubject,tidx));
            fprintf('Finish %.1fs simulation\n',time);
            fprintf('Mean accuracy %.2f \n',acc(numSubject+1,tidx));
            toc
        end
        %%  Save result to xls
        col_header = strsplit(num2str(timeSeq));     %Row cell array (for column labels)
        col_header = strcat(col_header,'s');
        row_header = cell(numSubject + 2,1);
        for i = 1:numSubject
            row_header{i}=['s' num2str(i)];     %Column cell array (for row labels)
        end
        row_header{numSubject+1} = 'Mean';
        row_header{numSubject+2} = 'Std';
        str = sprintf('%s_%s_startfrom%dms_filtered_numTrainTrial%d.xlsx',method,nameDataset,startTime*1000,numTrainTrial);
        filename = ['Result\' str];
        xlswrite(filename,acc,sprintf('Nh = %d',Nh),'B2');     %Write data
        xlswrite(filename,col_header,sprintf('Nh = %d',Nh),'B1');     %Write column header
        xlswrite(filename,row_header,sprintf('Nh = %d',Nh),'A2');      %Write row header
        fprintf('Finish N = %d simulation\n',Nh);
    end
    fprintf('Finish numTrainTrial = %d simulation\n',numTrainTrial);
    str = sprintf('%s_%s_startfrom%dms_filtered_numTrainTrial%d',method,nameDataset,startTime*1000,numTrainTrial);
    save(['Result\' str '_confusion.mat'],'rec_confusion');
end