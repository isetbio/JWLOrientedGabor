%% s_combine_eccbasedcoverage_fovea_accuracy

dataPth     = fullfile(ogRootPath,'data','classification',expName);

d = dir(fullfile(dataPth, '100trials_trainHighContrast', '*eccen0.00*'));

[~, reindex] = sort( str2double( regexp( {d.name}, '\d+', 'match', 'once' )));
d = d(reindex);

P =[];

for ii = 1:size(d,1)
    accuracy = load(fullfile(d(ii).folder, d(ii).name));
    
    P = [P accuracy.P];
    
end


fName   = sprintf('Classify_coneOutputs_contrast0.100_pa0_eye11_eccen0.00_defocus0.00_noise-random_sf4.00.mat');
save(fullfile(dataPth, '100trials_trainHighContrast',fName),'P');