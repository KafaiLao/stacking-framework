clear;clc;
timeSeq = 0.25:0.25:4;

for tidx = 1:length(timeSeq)
    str = sprintf('ECCA4_Benchmark_Nh5_time%d_filter1',timeSeq(tidx)*100);
    load(['..\Feature\' str]);
    
    str = sprintf('ECCA4_Benchmark_Nh5_time%d_filter1_cv',timeSeq(tidx)*100);
    load(['..\Feature\' str]);
    
    str = sprintf('ECCA4_Benchmark_Nh5_time%d_filter1_level1_Data.mat',timeSeq(tidx)*100);
    save(['..\Feature\' str],'trainFeatureSet','testFeatureSet','label');
    fprintf('Finsh time: %.2fs\n',timeSeq(tidx));
end