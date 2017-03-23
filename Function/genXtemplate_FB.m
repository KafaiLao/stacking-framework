function XtemplateFB = genXtemplate_FB(Xtemplate,filterBank)
% Generate filtered template for conventional filter bank apporaches
% [Input]
% Xtemplate [freqLength x sampleLength x channelLength]
% filterBank [cell array]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XtemplateFB = zeros([length(filterBank) size(Xtemplate)]);
for sb = 1:length(filterBank)
    for freq = 1:size(Xtemplate,1)
        XtemplateFB(sb,freq,:,:) = filtfilt(filterBank{sb},squeeze(Xtemplate(freq,:,:)));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Last modified 23/03/2017 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Ka Fai Lao, University of Macau %%%%%%%%%%%%%%%%%%%%%
end