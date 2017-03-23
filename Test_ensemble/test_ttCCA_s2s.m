clear;clc;addpath('Function');
%% Testing parameter
nameDataset = 'Benchmark';
method = 'ECCA4';
%% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load data to workspace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startTime = 0.14;
timeSeq = 0.25:0.25:4;
NhSeq = 5;
numSubject = 35;
freqLength = 40;
rec_confusion = zeros(length(NhSeq),length(timeSeq),numSubject,numSubject,freqLength,freqLength);
rec_Beta = zeros(length(NhSeq),length(timeSeq),numSubject,(numSubject+1)*3+1,freqLength);
rec_Constant = zeros(length(NhSeq),length(timeSeq),numSubject,freqLength);

Nhidx = 1;
Nh = NhSeq(Nhidx);
acc = zeros(numSubject+3,numSubject);
for tidx = 1:length(timeSeq)
    time = timeSeq(tidx);
    load(['Feature\' sprintf('%s_%s_Nh%d_time%d_filter%d.mat',method,nameDataset,Nh,time*100,1)],'featureSet','label');
    
    templateTags = [repmat(1:numSubject+1,1,3) 0];
    sourceSeq = 1:numSubject;
    trainLabel = reshape(label(1:numSubject-1,:,:),(numSubject-1)*size(label,2),size(label,3));
    testLabel = squeeze(label(1,:,:));
    rec = zeros(numSubject,numSubject);
    
    for testSubject = 1:numSubject
        for sourceSubject = 1:numSubject
            takeIdx = templateTags == sourceSubject;
            takeIdx(end) = true;
            testData = squeeze(featureSet(testSubject,:,takeIdx,:));
            
            % Recognize testData
            for trial = 1:size(testData,1)
                X = squeeze(testData(trial,:,:));
                score = mean(sign(X).*X.^2);
                [~,prediction] = max(score);
                [~,groundTruth] = max(testLabel(trial,:));
                if prediction == groundTruth, rec(sourceSubject,testSubject) = rec(sourceSubject,testSubject) + 1; end
                rec_confusion(Nhidx,tidx,testSubject,sourceSubject,groundTruth,prediction) = ...
                    rec_confusion(Nhidx,tidx,testSubject,sourceSubject,groundTruth,prediction) + 1;
            end
            fprintf('Accuracy of S%d at time %.2fs = %.2f\n',testSubject,time,rec(sourceSubject,testSubject)/size(testData,1)*100);
        end
    end
    acc(1:numSubject,1:numSubject) = rec/size(testData,1)*100;
    temp = acc(1:numSubject,1:numSubject).*~logical(eye(numSubject));
    acc(numSubject+1,1:numSubject) = max(temp);
    acc(numSubject+2,1:numSubject) = mean(acc(1:numSubject,1:numSubject));
    acc(numSubject+3,1:numSubject) = std(acc(1:numSubject,1:numSubject));
    fprintf('Finish %.1fs simulation\n',time);
    %%  Save result to xls
    row_header = cell(numSubject + 2,1);
    for i = 1:numSubject
        row_header{i}=['s' num2str(i)];     %Column cell array (for row labels)
    end
    col_header = row_header(1:numSubject)';
    row_header{numSubject+1} = 'MaxS';
    row_header{numSubject+2} = 'Mean';
    row_header{numSubject+3} = 'Std';
    str = sprintf('ttCCA_s2s_%s_%s_startfrom%dms_filtered.xlsx',method,nameDataset,startTime*1000);
    filename = ['Result\' str];
    xlswrite(filename,acc,sprintf('t = %.3fs',time),'B2');     %Write data
    xlswrite(filename,col_header,sprintf('t = %.3fs',time),'B1');     %Write column header
    xlswrite(filename,row_header,sprintf('t = %.3fs',time),'A2');      %Write row header
    fprintf('Finish N = %d simulation\n',Nh);
end
str = sprintf('ttCCA_s2s_%s_%s_startfrom%dms_filtered',method,nameDataset,startTime*1000);
save(['Result\' str '_confusion.mat'],'rec_confusion');
