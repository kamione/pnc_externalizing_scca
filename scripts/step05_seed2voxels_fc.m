%% Environment Setup ------------------------------------------------------
path_proj = '/projects/students/TYWong';
addpath(genpath(fullfile(path_proj, 'utilities', 'functions')));
addpath(fullfile(path_proj, 'utilities', 'spm12'));

path_wd        = fullfile(path_proj, '2020_pnc_fc_extbeh');
path_brainmask = fullfile(path_wd, 'masks', 'brainmask_JHU.nii');
path_output    = fullfile(path_wd, 'outputs');
cohort         = 'PNC';

%% Data I/O ---------------------------------------------------------------
% load table
path_data  = fullfile(path_wd, 'data', 'preprocessed', 'pnc_n1149_data.csv');
data = readtable(path_data);
data.SUBJID = string(data.SUBJID);
subjlist   = table2array(data(:, 1));
n_subj     = length(subjlist);

% check existence of files
flag = ones(n_subj, 1);
for ith_subj = 1:n_subj
    path_cohort_part1 = '/projects/PNC/data';
    path_cohort_part2  = 'restingfMRI/Preprocess_BBR/10_fMRI2Atlas';
    preprocessing_pipeline = '6motion_3tissue';
    
    subj_id = subjlist{ith_subj};
    path_vol = fullfile(path_cohort_part1, subj_id, path_cohort_part2, ...
                        preprocessing_pipeline, ...
                        [ subj_id '_4d_final_inAtlas_FS.img']);
             
    if ~exist(path_vol)
       flag(ith_subj) = 0;
    end
end


data_complete = data(logical(flag), :);
subjlist_complete = subjlist(logical(flag));

% save complete data to .mat
save(fullfile(path_wd, 'data', 'preprocessed', 'data_completed.mat'), 'data_complete');


%% Seed to Voxel Connectivity ---------------------------------------------
% seed Amygdala
seed         = 'dbd_amygdala';
path_roimask = fullfile(path_wd, 'outputs', 'ALEseedsJHU', sprintf('%s_MNI2JHU.img', seed));
amy_info = one_seed2voxel_corr(subjlist_complete, cohort, path_output, seed, path_brainmask);
save(fullfile(path_output, 'PNC_seed2voxelfc_summary', sprintf('%s_%s_info.mat', cohort, seed)), 'amy_info');

% seed IPL
seed         = 'dbd_IPL';
path_roimask = fullfile(path_wd, 'outputs', 'ALEseedsJHU', sprintf('%s_MNI2JHU.img', seed));
ipl_info = one_seed2voxel_corr(subjlist_complete, cohort, path_output, seed, path_brainmask);
save(fullfile(path_output, 'PNC_seed2voxelfc_summary',  sprintf('%s_%s_info.mat', cohort, seed)), 'ipl_info');

% seed MCC
seed         = 'adhd_mcc';
path_roimask = fullfile(path_wd, 'outputs', 'ALEseedsJHU', sprintf('%s_MNI2JHU.img', seed));
mcc_info = one_seed2voxel_corr(subjlist_complete, cohort, path_output, seed, path_brainmask);
save(fullfile(path_output, 'PNC_seed2voxelfc_summary', sprintf('%s_%s_info.mat', cohort, seed)), 'mcc_info');

% seed pre-SMA
seed         = 'adhd_preSMA';
path_roimask = fullfile(path_wd, 'outputs', 'ALEseedsJHU', sprintf('%s_MNI2JHU.img', seed));
presma_info = one_seed2voxel_corr(subjlist_complete, cohort, path_output, seed, path_brainmask);
save(fullfile(path_output, 'PNC_seed2voxelfc_summary', sprintf('%s_%s_info.mat', cohort, seed)), 'presma_info');

% seed Striatum
seed         = 'adhd_striatum';
path_roimask = fullfile(path_wd, 'outputs', 'ALEseedsJHU', sprintf('%s_MNI2JHU.img', seed));
striatum_info = one_seed2voxel_corr(subjlist_complete, cohort, path_output, seed, path_brainmask);
save(fullfile(path_output, 'PNC_seed2voxelfc_summary', sprintf('%s_%s_info.mat', cohort, seed)), 'striatum_info');
