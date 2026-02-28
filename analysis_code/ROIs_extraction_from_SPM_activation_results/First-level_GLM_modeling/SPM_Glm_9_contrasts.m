%% First-Level GLM analysis
%  Data analysis with SPM toolbox
%  Author: Yan Deng: 09.2024
%  Version: 5.0  Further make the script shorter. 
%  Version 6.0, add the Inference results, so it can show it automatically


clear;
% pwd

%% Data Path Settings
data_path = 'W:\01_COGMOT_Project\01_Single'; % windows system
data_dir = dir(data_path);
subjects = {data_dir(3:end).name};

%% Initiating SPM
%___________________
spm('Defaults','fMRI');
spm_jobman('initcfg');

for nsub = 1:length(subjects)      

  % The subject path:
    subject_path = fullfile(data_path, subjects{nsub});  
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLM SPECIFICATION, ESTIMATION, INFERENCE, RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
func_dataFolder = 'func';
f = spm_select('FPList', fullfile(subject_path, func_dataFolder), '^swrf.*');
rp_f = spm_select('FPList', fullfile(subject_path, func_dataFolder), '^rp_f.*');
% Should use FPlist, not List.  List can not find the right file.

    onsets_path = 'W:\Your_path\First-level_modeling\onsets'; 
    onsets_name = ['onsets_','sub-', num2str(nsub, '%02d'),'_GoNoGo_Test_7.mat']; % Corrected Onsets_IV
    onsets_file = fullfile(onsets_path,onsets_name);
    
clear matlabbatch
% Output Directory for GLM analysis
%--------------------------------------------------------------------------
glm_path = [subjects{nsub} '_GLM_'];  % This folder for all three tasks
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = cellstr(subject_path);
matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = glm_path ;

spm_jobman('run',matlabbatch);

clear matlabbatch
% Model Specification
%-------------------------------------------------------------------------
matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(fullfile(subject_path,glm_path));
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 0.85; %  RT = 850ms
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;  % 16 slices
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;  % half of the slices
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(f);

%%
matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('names', {}, 'onsets', {}, 'durations', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = cellstr(onsets_file);
matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = cellstr(rp_f);
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

% Model Estimation
%--------------------------------------------------------------------------
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% Contrasts
%--------------------------------------------------------------------------
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Go_single';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 0 0 0 0 0 0 0 0 0];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'NoGo_single';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 1 0 0 0 0 0 0 0 0];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'Single Nogo > Go';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [-1 1 0 0 0 0 0 0 0 0];
matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'Motor Dual>Single';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 0 -1 1 0 0 0 0 0 0];
matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'Motor_single';
matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0 0 1 0 0 0 0 0 0 0];
matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'Go_dual';
matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [0 0 0 1 0 0 0 0 0 0];
matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'Dual NoGO > Go';
matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [0 0 0 0 -1 1 0 0 0 0];
matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'Go Dual>Single';
matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [-1 0 0 1 0 0 0 0 0 0];
matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = 'Dual> Two Single';
matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [-1 0 -1 2 0 0 0 0 0 0];
matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0; % This line was gone before.

% Inference Results
%--------------------------------------------------------------------------
matlabbatch{4}.spm.stats.results.spmmat = cfg_dep('Contrasts: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{4}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{4}.spm.stats.results.conspec.extent = 0;
matlabbatch{4}.spm.stats.results.print = true;

% Inference Results_2
%--------------------------------------------------------------------------

matlabbatch{5}.spm.stats.results.spmmat = cfg_dep('Contrasts: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec.contrasts = 3;
matlabbatch{5}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{5}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{5}.spm.stats.results.conspec.extent = 0;
matlabbatch{5}.spm.stats.results.print = true; % YD added the results_2

% Inference Results_3
%--------------------------------------------------------------------------

matlabbatch{6}.spm.stats.results.spmmat = cfg_dep('Contrasts: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{6}.spm.stats.results.conspec.contrasts = 5;
matlabbatch{6}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{6}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{6}.spm.stats.results.conspec.extent = 0;
matlabbatch{6}.spm.stats.results.print = true; % YD added the results_3

% Inference Results_4
%--------------------------------------------------------------------------

matlabbatch{7}.spm.stats.results.spmmat = cfg_dep('Contrasts: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{7}.spm.stats.results.conspec.contrasts = 6;
matlabbatch{7}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{7}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{7}.spm.stats.results.conspec.extent = 0;
matlabbatch{7}.spm.stats.results.print = true; % YD added the results_4

% Inference Results_5
%--------------------------------------------------------------------------

matlabbatch{7}.spm.stats.results.spmmat = cfg_dep('Contrasts: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{7}.spm.stats.results.conspec.contrasts = 8;
matlabbatch{7}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{7}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{7}.spm.stats.results.conspec.extent = 0;
matlabbatch{7}.spm.stats.results.print = true; % YD added the results_5

spm_jobman('run', matlabbatch);

end