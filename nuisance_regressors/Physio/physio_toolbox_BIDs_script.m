% Batch PhysIO (RETROICOR) regressor generation for BIDs fMRI data
%
% Uses the TAPAS PhysIO Toolbox (via SPM12) to generate physiological
% noise regressors from cardiac (PPU) and respiratory recordings for a
% set of subjects, sessions, and runs.
%
% To use:
%   1. Set the paths and subject lists below
%   2. Choose which subject list to run (standard = single session, fixed runs; special = multiple sessions or a variable number of runs)
%   3. Run the script in MATLAB with SPM12 + PhysIO on the path
% 
% Note: with vendor = 'BIDs' and empty sampling_interval / relative_start_acquisition, a JSON sidecar is required. 
% The .json file must sit in the same folder as its .tsv.gz and share the same filename. 
% To use one shared JSON instead, set log_files.scan_timing to its path, or set the sampling interval and start time manually.
%
% You would need: MATLAB, SPM12, TAPAS PhysIO Toolbox (tested with PhysIO R2022a-v8.1.0)

clear;
clc;

base_path = '/path/to/your/BIDS_dataset';   % root of the BIDS dataset
spm_path  = '/path/to/spm12';               % SPM12 folder

% Standard subjects: single session, 6 runs each
standard_subjects = {'0001', '0002', '0003'};

% Special subjects: e.g. two fMRI sessions, or more than 6 runs
special_subjects  = {'0004', '0005'};

% Select which list to process
subjects = standard_subjects;
% subjects = special_subjects;

task_label   = 'flk';   % BIDS task label
acq_label    = 'hcp';   % BIDS acquisition label
n_slices     = 58;      % number of slices
TR           = 1.8;     % repetition time (s)
onset_slice  = 29;      % reference slice (typically n_slices / 2)
n_runs_std   = 6;       % number of runs for standard subjects

addpath(spm_path);
spm('defaults', 'FMRI');

for subj_idx = 1:length(subjects)
    subject_id = subjects{subj_idx};
    fprintf('Processing subject %s\n', subject_id);

    is_special = ismember(subject_id, special_subjects);

    if is_special
        session_dirs = dir(fullfile(base_path, sprintf('sub-%s', subject_id), 'ses-*'));

        if isempty(session_dirs)
            fprintf('No session directories for subject %s\n', subject_id);
            continue;
        end

        for sess_idx = 1:length(session_dirs)
            session_folder = fullfile(session_dirs(sess_idx).folder, session_dirs(sess_idx).name);
            func_folder = fullfile(session_folder, 'func');

            base_out_dir = fullfile(base_path, 'physio_regressors', sprintf('ST%s_%s', subject_id, session_dirs(sess_idx).name));

            cardiac_files = dir(fullfile(func_folder, sprintf('sub-%s_%s_task-%s_acq-%s_run-*_recording-cardiac_physio.tsv.gz', subject_id, session_dirs(sess_idx).name, task_label, acq_label)));

            fprintf('Processing session %s with %d runs\n', session_dirs(sess_idx).name, length(cardiac_files));

            if ~exist(base_out_dir, 'dir')
                mkdir(base_out_dir);
            end

            for run_idx = 1:length(cardiac_files)
                run_match = regexp(cardiac_files(run_idx).name, 'run-(\d+)', 'tokens');
                if ~isempty(run_match)
                    run_str = run_match{1}{1};
                else
                    run_str = sprintf('%02d', run_idx);
                end

                cardiac_file = fullfile(cardiac_files(run_idx).folder, cardiac_files(run_idx).name);
                resp_file = strrep(cardiac_file, 'recording-cardiac_physio.tsv.gz', 'recording-respiratory_physio.tsv.gz');
                bold_file = sprintf('%s/sub-%s/%s/func/sub-%s_%s_task-%s_acq-%s_run-%s_bold.nii.gz', base_path, subject_id, session_dirs(sess_idx).name, subject_id, session_dirs(sess_idx).name, task_label, acq_label, run_str);
                output_dir = fullfile(base_out_dir, ['run-' run_str]);

                fprintf('Processing run %s\n', run_str);

                clear nvols
                temp_dir = tempdir;
                temp_files = gunzip(bold_file, temp_dir);
                temp_nii = temp_files{1};

                vol_info = spm_vol(temp_nii);
                nvols = length(vol_info);

                delete(temp_nii);

                fprintf('Run %s: Found %d volumes\n', run_str, nvols);

                if ~exist(output_dir, 'dir')
                    mkdir(output_dir);
                end

                clear matlabbatch
                matlabbatch = {};
                matlabbatch{1}.spm.tools.physio.save_dir = {output_dir};
                matlabbatch{1}.spm.tools.physio.log_files.vendor = 'BIDs';
                matlabbatch{1}.spm.tools.physio.log_files.cardiac = {cardiac_file};
                matlabbatch{1}.spm.tools.physio.log_files.respiration = {resp_file};
                matlabbatch{1}.spm.tools.physio.log_files.scan_timing = {''};
                matlabbatch{1}.spm.tools.physio.log_files.sampling_interval = []; % Takes from json
                matlabbatch{1}.spm.tools.physio.log_files.relative_start_acquisition = [];
                matlabbatch{1}.spm.tools.physio.log_files.align_scan = 'last'; % Typically, aligning the last scan to the end of the logfile is beneficial, since start of logfile and scans might be shifted due to pre-scans
                matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nslices = n_slices;
                matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.NslicesPerBeat = [];
                matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.TR = TR;
                matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Ndummies = 0;
                matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nscans = nvols;
                matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.onset_slice = onset_slice; % Nslices/2
                matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.time_slice_to_slice = [];
                matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nprep = [];
                matlabbatch{1}.spm.tools.physio.scan_timing.sync.nominal = struct([]);
                matlabbatch{1}.spm.tools.physio.preproc.cardiac.modality = 'PPU';
                matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.min = 0.4;
                matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.file = sprintf('initial_cpulse_kRpeakfile_run%s.mat', run_str);
                matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.max_heart_rate_bpm = 90;
                matlabbatch{1}.spm.tools.physio.preproc.cardiac.posthoc_cpulse_select.off = struct([]);
                matlabbatch{1}.spm.tools.physio.preproc.respiratory.filter.passband = [0.01 2]; % Hz
                matlabbatch{1}.spm.tools.physio.preproc.respiratory.despike = false;
                matlabbatch{1}.spm.tools.physio.model.output_multiple_regressors = sprintf('multiple_regressors_run%s.txt', run_str);
                matlabbatch{1}.spm.tools.physio.model.output_physio = sprintf('physio_run%s.mat', run_str);
                matlabbatch{1}.spm.tools.physio.model.orthogonalise = 'none';
                matlabbatch{1}.spm.tools.physio.model.censor_unreliable_recording_intervals = true;
                matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.c = 3;
                matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.r = 4;
                matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.cr = 1;
                matlabbatch{1}.spm.tools.physio.model.rvt.no = struct([]);
                matlabbatch{1}.spm.tools.physio.model.hrv.no = struct([]);
                matlabbatch{1}.spm.tools.physio.model.noise_rois.no = struct([]);
                matlabbatch{1}.spm.tools.physio.model.movement.no = struct([]);
                matlabbatch{1}.spm.tools.physio.model.other.no = struct([]);
                matlabbatch{1}.spm.tools.physio.verbose.level = 2;
                matlabbatch{1}.spm.tools.physio.verbose.fig_output_file = sprintf('physio_run%s.jpeg', run_str);
                matlabbatch{1}.spm.tools.physio.verbose.use_tabs = false;

                % Run the batch
                try
                    spm_jobman('run', matlabbatch);
                    fprintf('Subject %s, Session %s, Run %s processed successfully\n', subject_id, session_dirs(sess_idx).name, run_str);
                catch ME
                    fprintf('ERROR processing Subject %s, Session %s, Run %s: %s\n', subject_id, session_dirs(sess_idx).name, run_str, ME.message);
                end
            end
        end

    else % standard subjects
        base_out_dir = fullfile(base_path, 'physio_regressors', sprintf('ST%s', subject_id));

        for run = 1:n_runs_std
            run_str = sprintf('%02d', run);

            cardiac_file = sprintf('%s/sub-%s/ses-01/func/sub-%s_ses-01_task-%s_acq-%s_run-%s_recording-cardiac_physio.tsv.gz', base_path, subject_id, subject_id, task_label, acq_label, run_str);
            resp_file = sprintf('%s/sub-%s/ses-01/func/sub-%s_ses-01_task-%s_acq-%s_run-%s_recording-respiratory_physio.tsv.gz', base_path, subject_id, subject_id, task_label, acq_label, run_str);
            bold_file = sprintf('%s/sub-%s/ses-01/func/sub-%s_ses-01_task-%s_acq-%s_run-%s_bold.nii.gz', base_path, subject_id, subject_id, task_label, acq_label, run_str);
            output_dir = fullfile(base_out_dir, ['run-' run_str]);

            fprintf('Processing run %s\n', run_str);

            clear nvols
            temp_dir = tempdir;
            temp_files = gunzip(bold_file, temp_dir);
            temp_nii = temp_files{1};

            vol_info = spm_vol(temp_nii);
            nvols = length(vol_info);

            delete(temp_nii);

            fprintf('Run %s: Found %d volumes\n', run_str, nvols);

            clear matlabbatch
            matlabbatch = {};
            matlabbatch{1}.spm.tools.physio.save_dir = {output_dir};
            matlabbatch{1}.spm.tools.physio.log_files.vendor = 'BIDs';
            matlabbatch{1}.spm.tools.physio.log_files.cardiac = {cardiac_file};
            matlabbatch{1}.spm.tools.physio.log_files.respiration = {resp_file};
            matlabbatch{1}.spm.tools.physio.log_files.scan_timing = {''};
            matlabbatch{1}.spm.tools.physio.log_files.sampling_interval = []; % Takes from json
            matlabbatch{1}.spm.tools.physio.log_files.relative_start_acquisition = [];
            matlabbatch{1}.spm.tools.physio.log_files.align_scan = 'last'; % Typically, aligning the last scan to the end of the logfile is beneficial, since start of logfile and scans might be shifted due to pre-scans
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nslices = n_slices;
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.NslicesPerBeat = [];
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.TR = TR;
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Ndummies = 0;
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nscans = nvols;
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.onset_slice = onset_slice; % Nslices/2
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.time_slice_to_slice = [];
            matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nprep = [];
            matlabbatch{1}.spm.tools.physio.scan_timing.sync.nominal = struct([]);
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.modality = 'PPU';
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.min = 0.4;
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.file = sprintf('initial_cpulse_kRpeakfile_run%s.mat', run_str);
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.max_heart_rate_bpm = 90;
            matlabbatch{1}.spm.tools.physio.preproc.cardiac.posthoc_cpulse_select.off = struct([]);
            matlabbatch{1}.spm.tools.physio.preproc.respiratory.filter.passband = [0.01 2];
            matlabbatch{1}.spm.tools.physio.preproc.respiratory.despike = false;
            matlabbatch{1}.spm.tools.physio.model.output_multiple_regressors = sprintf('multiple_regressors_run%s.txt', run_str);
            matlabbatch{1}.spm.tools.physio.model.output_physio = sprintf('physio_run%s.mat', run_str);
            matlabbatch{1}.spm.tools.physio.model.orthogonalise = 'none';
            matlabbatch{1}.spm.tools.physio.model.censor_unreliable_recording_intervals = true;
            matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.c = 3;
            matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.r = 4;
            matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.cr = 1;
            matlabbatch{1}.spm.tools.physio.model.rvt.no = struct([]);
            matlabbatch{1}.spm.tools.physio.model.hrv.no = struct([]);
            matlabbatch{1}.spm.tools.physio.model.noise_rois.no = struct([]);
            matlabbatch{1}.spm.tools.physio.model.movement.no = struct([]);
            matlabbatch{1}.spm.tools.physio.model.other.no = struct([]);
            matlabbatch{1}.spm.tools.physio.verbose.level = 2;
            matlabbatch{1}.spm.tools.physio.verbose.fig_output_file = sprintf('physio_run%s.jpeg', run_str);
            matlabbatch{1}.spm.tools.physio.verbose.use_tabs = false;

            try
                spm_jobman('run', matlabbatch);
                fprintf('Subject %s, Run %s processed successfully\n', subject_id, run_str);
            catch ME
                fprintf('ERROR processing Subject %s, Run %s: %s\n', subject_id, run_str, ME.message);
            end
        end
    end
end

fprintf('All subjects processed\n');
