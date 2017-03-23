function template = genSinTemplate(stimuFreq, fsample, time, Nharmonic)
% Generate artifical sinosudial template for conventional CCA
% [Input]
% stimuFreq [number of stimulus frequency x 1]: visual stimulus frequencies
% fsample [scalar]: sampling rate of the EEG signal
% time [scalar]: Time length (in second)of the output template
% Nharmonic [scalar]: Number of harmonic we consider in the template
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
template = zeros(length(stimuFreq),floor(fsample*time),2*Nharmonic);
timeSeq = 0:1/fsample:time-1/fsample;
for freqIdx = 1:length(stimuFreq)
    for i = 1:Nharmonic
        template(freqIdx,:,(i-1)*2+1) = sin(2*pi*i*stimuFreq(freqIdx)*timeSeq);
        template(freqIdx,:,i*2) = cos(2*pi*i*stimuFreq(freqIdx)*timeSeq);
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Last modified 23/03/2017 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Ka Fai Lao, University of Macau %%%%%%%%%%%%%%%%%%%%%