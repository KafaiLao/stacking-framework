clear;clc;
timeSeq = 0.25:0.25:4;
OrgTemplateTags = [repmat(1:36,1,3) 0];

for tidx = 1:length(timeSeq)
    targetFeatureSet = zeros(35,240,103,40);
    str = sprintf('ECCA4_Benchmark_Nh5_time%d_filter1',timeSeq(tidx)*100);
    load(['Feature\' str]);
    %%%%%%%%% Take out the target subject feature (going to classify)%%%%%
    for targetSubject = 1:35
        takeIdx = find(OrgTemplateTags ~= targetSubject & OrgTemplateTags ~= 36);
        targetFeatureSet(targetSubject,:,:,:) = featureSet(targetSubject,:,takeIdx,:);
    end
    %%%%%%%%% Take out the target subject feature (going to classify)%%%%%
    str = sprintf('ECCA4_Benchmark_Nh5_time%d_filter1_cv',timeSeq(tidx)*100);
    load(['Feature\' str]);
    sourceFeatureSet = featureSet;
    
    str = sprintf('ECCA4_Benchmark_Nh5_time%d_filter1_final.mat',timeSeq(tidx)*100);
    save(['Feature\' str],'sourceFeatureSet','targetFeatureSet','label');
    fprintf('Finsh time: %.2fs\n',timeSeq(tidx));
end