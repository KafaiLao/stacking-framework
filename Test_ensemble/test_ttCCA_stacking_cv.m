clear;clc;addpath('..\Function');
%% Testing parameter
nameDataset = 'Benchmark';
method = 'ECCA4';
%% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load data to workspace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startTime = 0.14;
timeSeq = 0.25:0.25:4;
NhSeq = 5;
Tempature = 1;
numSubject = 35;
freqLength = 40;
trialLength = 6;
rec_confusion = zeros(length(NhSeq),length(timeSeq),numSubject,freqLength,freqLength);
rec_Beta = zeros(length(NhSeq),length(timeSeq),numSubject,numSubject*3+1,freqLength);
rec_Constant = zeros(length(NhSeq),length(timeSeq),numSubject,freqLength);
rec_largeModelOutput = zeros(length(NhSeq),length(timeSeq),numSubject,freqLength*trialLength,freqLength);
for Nhidx = 1:length(NhSeq)
    Nh = NhSeq(Nhidx);
    acc = zeros(numSubject+2,length(timeSeq));
    for tidx = 1:length(timeSeq)
        time = timeSeq(tidx);
        load(['..\Feature\' sprintf('%s_%s_Nh%d_time%d_filter%d_level1_Data.mat',method,nameDataset,Nh,time*100,1)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Convert the correlation coefficient to probability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        trainFeatureSetProb = zeros(size(trainFeatureSet));
        for i = 1:size(trainFeatureSet,1)
            for j = 1:size(trainFeatureSet,2)
                for k = 1:size(trainFeatureSet,3)
                    trainFeatureSetProb(i,j,k,:) = softmaxfun(trainFeatureSet(i,j,k,:),Tempature);
                end
            end
        end
        
        testFeatureSetProb = zeros(size(testFeatureSet));
        for i = 1:size(testFeatureSet,1)
            for j = 1:size(testFeatureSet,2)
                for k = 1:size(testFeatureSet,3)
                    testFeatureSetProb(i,j,k,:) = softmaxfun(testFeatureSet(i,j,k,:),Tempature);
                end
            end
        end        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        templateTags = [repmat(1:numSubject,1,3) 0];
        sourceSeq = 1:numSubject;
        trainLabel = reshape(label(1:numSubject-1,:,:),(numSubject-1)*size(label,2),size(label,3));
        testLabel = squeeze(label(1,:,:));
        rec = zeros(numSubject,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        for testSubject = 1:numSubject
            takeIdx = find(templateTags ~= testSubject);
            % Take out the test (target subject) data
            testData = squeeze(testFeatureSetProb(testSubject,:,takeIdx,:));
            % Take out the train (source subjects) data
            trainData = squeeze(trainFeatureSetProb(sourceSeq ~= testSubject,:,takeIdx,:));
            trainData = reshape(trainData,size(trainData,1)*size(trainData,2),size(trainData,3),size(trainData,4));
            
            Beta = zeros(size(trainData,2)+1,freqLength);
            
            % Solve L linear regression problems
            for lProblem = 1:freqLength
                lTrainData = squeeze(trainData(:,:,lProblem));
                lTrainLabel = squeeze(trainLabel(:,lProblem));
                sw = ones(length(lTrainLabel),1);
                sw(logical(lTrainLabel)) = 1;
                sw = sw/sum(sw);
                Beta(:,lProblem) = glmfit(lTrainData,lTrainLabel,'normal','weights',sw);
            end
            rec_Beta(Nhidx,tidx,testSubject,templateTags ~= testSubject,:) = Beta(2:end,:);
            rec_Constant(Nhidx,tidx,testSubject,:) = Beta(1,:);
            
            score = zeros(freqLength,1);
            % Recognize testData
            for trial = 1:size(testData,1)
                for lProblem = 1:freqLength
                    feature = squeeze(testData(trial,:,lProblem));
                    score(lProblem) = glmval(Beta(:,lProblem),feature,'identity');
%                     score(lProblem) = feature*Beta(:,lProblem);
                end
                rec_largeModelOutput(Nhidx,tidx,testSubject,trial,:) = score;
                [~,prediction] = max(score);
                [~,groundTruth] = max(testLabel(trial,:));
                if prediction == groundTruth, rec(testSubject) = rec(testSubject) + 1; end
                rec_confusion(Nhidx,tidx,testSubject,groundTruth,prediction) = rec_confusion(Nhidx,tidx,testSubject,groundTruth,prediction) + 1;
            end
            fprintf('Accuracy of S%d at time %.2fs = %.2f\n',testSubject,time,rec(testSubject)/size(testData,1)*100);
        end
        acc(1:numSubject,tidx) = rec/size(testData,1)*100;
        acc(numSubject+1,tidx) = mean(acc(1:numSubject,tidx));
        acc(numSubject+2,tidx) = std(acc(1:numSubject,tidx));
        fprintf('Finish %.1fs simulation\n',time);
        fprintf('Mean accuracy %.2f \n',acc(numSubject+1,tidx));
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
    str = sprintf('ttCCA_stack_%s_%s_startfrom%dms_filtered_CV.xlsx',method,nameDataset,startTime*1000);
    filename = ['Result\' str];
    xlswrite(filename,acc,sprintf('Nh = %d',Nh),'B2');     %Write data
    xlswrite(filename,col_header,sprintf('Nh = %d',Nh),'B1');     %Write column header
    xlswrite(filename,row_header,sprintf('Nh = %d',Nh),'A2');      %Write row header
    fprintf('Finish N = %d simulation\n',Nh);
end
str = sprintf('ttCCA_stack_%s_%s_startfrom%dms_filtered_CV',method,nameDataset,startTime*1000);
save(['..\Result\' str '_confusion.mat'],'rec_confusion','rec_Beta','rec_Constant','rec_largeModelOutput');
