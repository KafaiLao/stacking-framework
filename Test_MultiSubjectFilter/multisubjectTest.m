clear;clc;close all;
addpath('C:\Users\controllab\Dropbox\02Academic Work\01SSVEP\Toolbox\eeglab13_5_4b');
addpath('Function\');
% Input: X (1375x8x35) one-trial EEG signal from each subject
load('ProcessedData\benchmark_prefilter1_CAR0.mat')
ssvep_all = cat(5,allData{:});

startIdx = round(0.14*fsample);
tSeq = 0.25:0.25:4;

sinTemplate = genSinTemplate(stimuFreq,fsample,4,5);
% for subject = 1:size(ssvep_all,5)
%     ssvep = squeeze(ssvep_all(:,:,startIdx+1:startIdx+1000,:,subject));
%     template = squeeze(mean(ssvep));
%     SF = canoncorr(squeeze(template(1,:,:)),squeeze(sinTemplate(1,:,:)));
%     W(:,subject) = SF(:,1);
% end
% 
% for tidx = 1:length(tSeq)
% % time = 1;
% time = tSeq(tidx);
% sinTemplate = genSinTemplate(stimuFreq,fsample,time,5);
% %% First have a pliot test
% % figure;
% rec = zeros(size(ssvep_all,5),1);
% 
% for subject = 1:size(ssvep_all,5)
%     ssvep = squeeze(ssvep_all(:,:,startIdx+1:startIdx+floor(time*fsample),:,subject));
% %     template = squeeze(mean(ssvep));
% %     SF = canoncorr(squeeze(template(1,:,:)),squeeze(sinTemplate(1,:,:)));
% %     w = SF(:,1);
%     w = W(:,subject);
% %     subplot(6,6,subject);topoplot([zeros(19,1);w], 'C:\SSVEP\Tsinghua database 2016\61-channels.loc','electrodes','off');title(sprintf('S%d',subject));
%     for trial = 1:size(ssvep_all,1)
%         for freq = 2:size(ssvep_all,2)
%             Xnew = squeeze(ssvep(trial,freq,:,:));
%             p = zeros(size(ssvep_all,2),1);
%             for tFreq = 1:size(ssvep_all,2)
%                 [~,~,r] = canoncorr(Xnew*w,squeeze(sinTemplate(tFreq,:,:)));
%                 p(tFreq) = r(1);
%             end
%             [~,loc] = max(p);
%             if loc == freq,rec(subject) = rec(subject) + 1;end;
%         end
%     end
% end
% meanacc_samew(tidx) = mean(rec/2.34);
% end
% stop
%% Take out one
figure;
train_freq = 1;
train_ssvep = squeeze(ssvep_all(:,train_freq,startIdx+1:startIdx+floor(time*fsample),:,:));
Y = squeeze(sinTemplate(train_freq,:,:));
for subject = 1:size(ssvep_all,5)
    Zall(:,:,subject) = squeeze(mean(train_ssvep(:,:,:,subject)));   
    [Wtemp,Vtemp] = canoncorr(squeeze(Zall(:,:,subject)),Y);
%     Wall(:,subject) = Wtemp(:,1)/sqrt(size(train_ssvep,2)-1);
%     Vall(:,subject) = Vtemp(:,1)/sqrt(size(train_ssvep,2)-1);
    Wall(:,subject) = Wtemp(:,1);
%     Vall(:,subject) = Vtemp(:,1);
%     subplot(6,6,subject);topoplot([zeros(19,1);Wall(:,subject)], 'C:\SSVEP\Tsinghua database 2016\61-channels.loc','electrodes','off');
end
% Wall = 1e-4*randn(9,35);
% Vall = 1e-4*randn(10,35);
v = 1e-4*randn(10,1);
w0 = mean(Wall,2);
% subplot(6,6,36);topoplot([zeros(19,1);w0], 'C:\SSVEP\Tsinghua database 2016\61-channels.loc','electrodes','off');
stop
betaSeq = 0.7;
for idx = 1:length(betaSeq)
beta = betaSeq(idx);
alpha = 0.001;
mu = 1;
ls = randn(35,1);
ps = randn(1,1);
for i = 1:1000
    if i == 200, alpha = alpha/10;end;
    output(i) = objFun_fixV(Zall,Y,Wall,v,w0,beta);
    [Wall,v,w0,ls,ps,dws,dvs,dw0] = updateGrad_fixV(Zall,Y,Wall,v,w0,beta,alpha,mu,ls,ps);
%     lambda = 1.01*lambda;
end
plot(output)

rec = zeros(size(ssvep_all,5),1);
for subject = 1:size(ssvep_all,5)
    ssvep = squeeze(ssvep_all(:,:,startIdx+1:startIdx+floor(time*fsample),:,subject));
    w = Wall(:,subject);
    for trial = 1:size(ssvep_all,1)
        for freq = 1:size(ssvep_all,2)
            Xnew = squeeze(ssvep(trial,freq,:,:));
            for tFreq = 1:size(ssvep_all,2)
                [~,~,r] = canoncorr(Xnew*w,squeeze(sinTemplate(tFreq,:,:)));
                p(tFreq) = r(1);
            end
            [~,loc] = max(p);
            if loc == freq,rec(subject) = rec(subject) + 1;end;
        end
    end
end
meanacc(idx) = mean(rec/2.4);
end
