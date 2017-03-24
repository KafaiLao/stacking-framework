function initialize_EEGLAB()
username=getenv('USERNAME');
str = sprintf('C:\\Users\\%s\\Dropbox\\02Academic Work\\01SSVEP\\Toolbox\\eeglab13_5_4b\\functions',username);
addpath(genpath(str));
end