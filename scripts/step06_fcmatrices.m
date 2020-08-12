path_proj = '/projects/students/TYWong';
addpath(genpath(fullfile(path_proj, 'utilities', 'functions')));
addpath(fullfile(path_proj, 'utilities', 'spm12'));

path_wd = fullfile(path_proj, '2020_pnc_fc_extbeh');

load(fullfile(path_wd, 'data', 'preprocessed', 'data_completed.mat'));
data = data_complete;

% change the SubjID to string
subjlist       = table2array(data(:, 1));
path_output    = fullfile(path_wd, 'outputs');
path_brainmask = fullfile(path_wd, 'masks', 'brainmask_JHU.nii');
cohort         = 'PNC';

seed = 'adhd_mcc';
mcc_fcmaps = one_extract_fcmaps(subjlist, seed, cohort, path_output, path_brainmask);
mcc_fcmaps_shen = one_extract_fcmaps_shen(subjlist, seed, cohort, path_output, path_brainmask);
%mcc_fcmaps_adjusted = one_fcmaps_adjust4covs(mcc_fcmaps, seed, cohort, covs, path_output, scale_idx, cat_idx);

seed = 'adhd_preSMA';
presma_fcmaps = one_extract_fcmaps(subjlist, seed, cohort, path_output, path_brainmask);
presma_fcmaps_shen = one_extract_fcmaps_shen(subjlist, seed, cohort, path_output, path_brainmask);
csvwrite('presma_fcmaps_shen.csv', presma_fcmaps_shen);
%presma_fcmaps_adjusted = one_fcmaps_adjust4covs(mcc_fcmaps, seed, cohort, covs, path_output, scale_idx, cat_idx);

seed = 'adhd_striatum';
stiatum_fcmaps = one_extract_fcmaps(subjlist, seed, cohort, path_output, path_brainmask);
striatum_fcmaps_shen = one_extract_fcmaps_shen(subjlist, seed, cohort, path_output, path_brainmask);
csvwrite('striatum_fcmaps_shen.csv', striatum_fcmaps_shen);
%stiatum_fcmaps_adjusted = one_fcmaps_adjust4covs(mcc_fcmaps, seed, cohort, covs, path_output, scale_idx, cat_idx);

seed = 'dbd_amygdala';
amy_fcmaps = one_extract_fcmaps(subjlist, seed, cohort, path_output, path_brainmask);
amy_fcmaps_shen = one_extract_fcmaps_shen(subjlist, seed, cohort, path_output, path_brainmask);
csvwrite('amy_fcmaps_shen.csv', amy_fcmaps_shen);
%amy_fcmaps_adjusted = one_fcmaps_adjust4covs(mcc_fcmaps, seed, cohort, covs, path_output, scale_idx, cat_idx);

seed = 'dbd_IPL';
ipl_fcmaps = one_extract_fcmaps(subjlist, seed, cohort, path_output, path_brainmask);
ipl_fcmaps_shen = one_extract_fcmaps_shen(subjlist, seed, cohort, path_output, path_brainmask);
csvwrite('ipl_fcmaps_shen.csv', ipl_fcmaps_shen);
%ipl_fcmaps_adjusted = one_fcmaps_adjust4covs(mcc_fcmaps, seed, cohort, covs, path_output, scale_idx, cat_idx);









