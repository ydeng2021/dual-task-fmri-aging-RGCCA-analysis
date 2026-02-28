%% 2nd_level Group Analysis Batch
% Version 3.0 for group one-sample t-test
% Test on All Subjects
% Yan Deng
% 30.05.2025 set K threshold = 20 after fMRIPrep 

%-----------------------------------------------------------------------

clear;
% Add the path
%addpath 'W:\01_COGMOT_Project\03_programs\2nd-level_analysis_24'

% Specify subject to process
data_path = 'W:\01_COGMOT_Project\01_fmriprep_single';
subjects_dir = dir(data_path);
subjects = {subjects_dir(3:end).name};
% {'sub-01'}    {'sub-02'}    {'sub-03'}    {'sub-04'}    {'sub-05'}    {'sub-06'}    {'sub-07'}

% Output Directory for 2nd_level analysis:
analysis_path ='W:\01_COGMOT_Project\02_GroupLevel\After_fMRIPrep';
contrast_path = fullfile(analysis_path, 'All_subjects_afterQC'); 

contrast = {'con_0001','con_0002', 'con_0003','con_0004','con_0005','con_0006', 'con_0007', 'con_0008','con_0009', 'con_0010'};
contrast_name = {'Single_Go', '-Single_NoGo','Single GoNoGo', 'Dual>Single Motor', 'Single Motor',...
   'Dual Go', 'Dual GoNoGO', 'Dual_Go > Single_Go', 'Dual > Two Singles', 'Dual > Two Singles (-1, -1)'};
% contrast 1 = 'Single_Go';
% contrast 2 = '- Single_NoGo';
% contrast 3 = 'Single GoNoGo';
% contrast 4 = 'Dual>Single_Motor';
% contrast 5 = 'Motor_single';
% contrast 6 = 'Dual_Go';
% contrast 7 = 'Dual GoNoGO';
% contrast 8 = 'Dual_Go >Single_Go';
% contrast 9 = 'Dual> Two Singles (-0.5, -0.5)';
% contrast 10 = 'Dual> Two Singles (-1, -1)';

% Import subjects summary table
% 13.04.2025 update this file further, exclude sub-08, 14, 28, 29
csvfile = 'subjects_summary_64_young_old_QC_NC.csv';
Subjects_summary= import_csv(csvfile);

ii = 1; % Choose the contrast

%% Initiating SPM
%___________________
spm('Defaults','fMRI');
spm_jobman('initcfg');

contrast_folder = ['After_fmriprep_one-sample_t_' contrast{ii}];

    clear matlabbatch
    % Output Directory for 2nd_level analysis: one-sample t-test
    %--------------------------------------------------------------------------  
    
    matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = cellstr(contrast_path);
    matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = contrast_folder;
    
    spm_jobman('run',matlabbatch);

%-----------------------------------------------------------------------
% One-sample t-test design
cnt =1; 
for nsub = 1:length(subjects)
    TF = (Subjects_summary{nsub,6} == 'Old');
    FT = (Subjects_summary{nsub,6} == 'Young');
    if  FT == 1  % ||   TF == 1      % Include both old and young subjects
        condition1_imgs{cnt,1}=fullfile(data_path,subjects{nsub},['GLM_fMRIPrep_' subjects{nsub}],[contrast{ii} '.nii,1']);
        cnt = cnt +1;
    end
end

%%
% Initialize batch
clear matlabbatch

matlabbatch{1}.spm.stats.factorial_design.dir = {fullfile(contrast_path,contrast_folder)};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = condition1_imgs;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

%-----------------------------------------------------------------------
% Estimation
%------------------------------------------------------------------------
matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(contrast_path,contrast_folder, 'SPM.mat')};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% Inference
%---------------------------------------------------------------------
matlabbatch{3}.spm.stats.con.spmmat =  {fullfile(contrast_path,contrast_folder, 'SPM.mat')};
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = contrast_name{ii};
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1; % To test for positive effects
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = ['- ' contrast_name{ii}];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = -1; % To test for negative effects
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;  % p-value threshold
matlabbatch{4}.spm.stats.results.conspec.extent = 20; % change from 0 to 20
matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{4}.spm.stats.results.units = 1;
matlabbatch{4}.spm.stats.results.export{1}.ps = true;

% run the batch
spm_jobman('run',matlabbatch);
