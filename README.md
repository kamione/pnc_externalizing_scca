## Information

- Authors: Ting-Yat Wong
- Emails: tywong@gmail.com
- Last Updates: 12-Aug-2020

## Folder Structure

### scripts
+ Step01: roi_mni2jhu.sh
    + transform ALE seeds from a MNI to JHU format
+ Step02: seed2voxels_fc.m
    + calculate individual seed-based FC maps to each voxel
+ Step03: fcmatrices.m
    + prepare FC matrices of each seed for SCCA analyses
+ Step04: prep_data4fmri.R
    + grab the subjects from raw files
+ Step05: combine_behbrain.R
    + combine behavior and brain data for SCCA analyses to R data
+ Step06: scca_adhd_rois.R
    + run SCCA on ADHD-related brain hubs and generate relevant figures of results
        1. regressing out variables
        2. split into discovery and replication
        3. SCCA main analysis
        4. permutation analysis for no. of modes
        5. boostrapping analysis for stability
+ Step07: scca_dbd_rois.R
    + run SCCA on DBD-related brain hubs and generate relevant figures of results
        1. regressing out variables
        2. split into discovery and replication
        3. SCCA main analysis
        4. permutation analysis for no. of modes
        5. boostrapping analysis for stability

### utilities
+ Folder: R
    + R functions used for SCCA analyses
