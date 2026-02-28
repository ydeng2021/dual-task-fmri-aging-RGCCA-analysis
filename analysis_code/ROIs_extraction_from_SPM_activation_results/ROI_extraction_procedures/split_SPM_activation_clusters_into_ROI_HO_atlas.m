%% Split SPM Brain Activation Mask in ROIs by Harvard Oxford (HO) atlas 
% Yan Deng
% 08.05.2025 Version 1.0
% 31.05.2025, Version 2.0 
% 12.11.2025, Version 3.0, Single Go, NoGo and Motor

% Reslice atlas to match ROI
data_path = 'Z:\01_COGMOT_Project\02_GroupLevel\After_fMRIPrep_02012025\All_subjects_afterQC\After_fmriprep_one-sample_t_test_AllSubjects_QC_con_0005\';
atlas_path = 'Z:\01_COGMOT_Project\02_GroupLevel\After_fMRIPrep_02012025\All_subjects_afterQC\After_fmriprep_one-sample_t_test_AllSubjects_QC_con_0006\';

cd(data_path); % Change to the data folder

% Set filenames
maskfilename = 'Single_Motor_all_clusters.nii'; % From step 1
ref_image = [data_path, maskfilename];      % reference = our ROI Mask
src_image = [atlas_path 'atlas.nii'];       % source = Harvard Oxford atlas

% Create file list for reslicing
reslice_flags = struct( ...
    'interp', 0, ...           % Nearest neighbor (preserves integer labels)
    'wrap', [0 0 0], ...
    'mask', 0, ...
    'which', 1, ...
    'mean', 0 );

spm_reslice({ref_image, src_image}, reslice_flags);
% SPM12: spm_reslice (v7141)     
% This generates: ratlas.nii in the atlas folder. "r" stands for "resliced"
% ratlas.nii has the same voxel size, orientation, and dimensions as source image.

%% Step 2: Use the resliced AAL 
% Load ROI and HO atlas
roi_vol = spm_vol(ref_image);
roi_data = spm_read_vols(roi_vol);

% Apply Harvard-Oxford atlas to split the ROIs mask
atlas_vol = spm_vol([atlas_path 'ratlas.nii']); % Use the resliced atlas
atlas_data = spm_read_vols(atlas_vol);

% Import atlas labels, short form
labels = importdata([atlas_path 'short_labels.txt']);   % 132x1 cell

unique_ids = unique(atlas_data(:));
unique_ids(unique_ids==0) = [];  % remove background
% 132 labels

% Create output dir, the same name as the maskfile name
[~, output_dir, ~] = fileparts(maskfilename);
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

min_voxels = 20;  % or 10/30 depending on your judgment

% Loop through each overlapping HO region
for i = 1:length(unique_ids)
    id = unique_ids(i);
    mask = (roi_data > 0) & (atlas_data == id);
    if nnz(mask) > min_voxels
        out = roi_vol;
        safe_name = strrep(labels{id}, ' ', '_');
        out.fname = fullfile(output_dir, ['SM_ROI_' safe_name '.nii']);
        spm_write_vol(out, double(mask));
        fprintf('Saved: %s\n', out.fname);
    end
end


%% Count Voxels per ROI Mask
% Folder where ROI.nii files are stored
roi_dir = output_dir;  % change if needed
roi_files = dir(fullfile(roi_dir, '*.nii'));

roi_names = cell(1,length(roi_files));
voxel_counts = zeros(1,length(roi_files));

for i = 1:length(roi_files)
    fname = fullfile(roi_dir, roi_files(i).name);
    V = spm_vol(fname);
    Y = spm_read_vols(V);    
    roi_names{i} = roi_files(i).name;
    voxel_counts(i) = nnz(Y);  % count non-zero voxels
end

% Display summary
T = table(roi_names', voxel_counts', 'VariableNames', {'ROI', 'NumVoxels'});
disp(T);

% Save to CSV
writetable(T, fullfile(roi_dir, [output_dir '_into_ROIs.csv']));

% Done at 15:20, 09.05.2025  Split all clusters into HO atlas ROIs