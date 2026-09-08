function preprocessing_SPM12_univariate(subject_ids)
% Preprocessing pipeline in SPM12 for fMRI univariate analysis.
%
% Runs the following steps:
%   1. DICOM to NIfTI conversion (dcm2niix)
%   2. Slice timing correction
%   3. Fieldmap estimation (VDM file creation)
%   4a. Realign & Unwarp (with fieldmap)
%   4b. Realign only (without fieldmap)
%   5. Spatial smoothing (multiple kernel sizes)
%
% Based on a script by David Wisniewski:
%   https://github.com/CCN-github/fMRI-preprocessing-SPM12/blob/master/preprocessing_BIDS_SPM12.m
%
% USAGE:
%   preprocessing_SPM12_univariate([101, 102, 103])
%
% INPUT:
%   subject_ids  - numeric array of subject IDs to preprocess

%% FLAGS — set to 1 to run a step, 0 to skip
% Run in order: slice time → fieldmap → realign&unwarp → smooth
do_slice_time     = 0;   % slice time correction
do_fieldmap       = 0;   % calculate fieldmap (VDM file)
do_realign_unwarp = 0;   % realignment + unwarping (requires fieldmap)
do_realign_only   = 0;   % realignment only (use when no fieldmap available)
do_smooth         = 1;   % spatial smoothing

%% USER SETTINGS
% --- Paths ------------------------------------------------------------------
SPM_PATH      = '/path/to/spm12';        % SPM12 folder
DCM2NIIX_PATH = '/path/to/dcm2niix';     % dcm2niix
BASE_PATH     = '/path/to/your/data';    % root data directory
% Expected structure under BASE_PATH:
% <BASE_PATH>/<subject_id>/scans/        — raw DICOM folders
% <BASE_PATH>/<subject_id>/preprocessed/ — created by this script

% Acquisition parameters
nruns = 6;   % number of functional runs

% Slice timing
param.st.nslices  = 58; % number of slices per volume
param.st.tr       = 1.8; % TR in seconds
param.st.ta       = 1.8 - (1.8 / 58); % TA = TR - (TR / nslices)
param.st.so       = [0 0.9225 0.0625 0.9825 0.1225 1.045 0.185 1.105 0.2475 1.1675 ...
                     0.3075 1.2275 0.37 1.29 0.43 1.35 0.4925 1.4125 0.5525 1.475 ...
                     0.615 1.535 0.675 1.5975 0.7375 1.6575 0.8 1.72 0.86 0 0.9225 ...
                     0.0625 0.9825 0.1225 1.045 0.185 1.105 0.2475 1.1675 ...
                     0.3075 1.2275 0.37 1.29 0.43 1.35 0.4925 1.4125 0.5525 1.475 ...
                     0.615 1.535 0.675 1.5975 0.7375 1.6575 0.8 1.72 0.86] * 1000;
                   % Interleaved slice acquisition order (from .json)
param.st.refslice = max(param.st.so) / 2;  % reference slice (middle of acquisition)

% Smoothing kernels (mm) — one smoothed output is produced per kernel
param.smooth.kernels = [2, 3, 4];

addpath(SPM_PATH);
spm('defaults', 'fmri');
spm_jobman('initcfg');

for subject_idx = 1:length(subject_ids)
    current_id   = subject_ids(subject_idx);
    subject_path = fullfile(BASE_PATH, num2str(current_id));
    fprintf('\n=== Preprocessing subject: %d ===\n', current_id);

    %% ------------------------------------------------------------------
    %  PART 1: DICOM → NIfTI conversion
    %  Converts raw DICOM files to compressed NIfTI (.nii.gz) and then
    %  gunzips them so SPM can read them
    % -------------------------------------------------------------------
    output_base = fullfile(subject_path, 'preprocessed');
    if ~exist(output_base, 'dir'), mkdir(output_base); end

    for d = {'func', 'fmap_mag', 'fmap_phase', 'anat'}
        dir_path = fullfile(output_base, d{1});
        if ~exist(dir_path, 'dir'), mkdir(dir_path); end
    end

    % Functional runs
    func_dirs = dir(fullfile(subject_path, 'scans', '*-func_bold_acq_hcp_task_flk_run_*'));
    func_dirs = func_dirs(~contains({func_dirs.name}, 'WIP_PMU'));  % remove PMU noise folders

    if ~isempty(func_dirs)
        % Sort runs by scan number (first numeric token in folder name)
        scan_numbers = zeros(length(func_dirs), 1);
        for i = 1:length(func_dirs)
            parts = split(func_dirs(i).name, '-');
            scan_numbers(i) = str2double(parts{1});
        end
        [~, sort_idx] = sort(scan_numbers);
        func_dirs = func_dirs(sort_idx);

        for i = 1:length(func_dirs)
            dicom_path = fullfile(subject_path, 'scans', func_dirs(i).name, 'resources', 'DICOM', 'files');
            run_name   = sprintf('run-%d', i);
            cmd = sprintf('"%s" -z y -w 1 -f "%s" -o "%s" "%s"', ...
                DCM2NIIX_PATH, run_name, fullfile(output_base, 'func'), dicom_path);
            [status, result] = system(cmd);
            if status ~= 0
                warning('dcm2niix failed for run %d, subject %d:\n%s', i, current_id, result);
            end
        end
    end

    % Fieldmaps
    has_fieldmaps = false;
    fmap_dirs = dir(fullfile(subject_path, 'scans', '*-fmap_acq_boldGRE'));

    if ~isempty(fmap_dirs)
        has_fieldmaps = true;
        scan_numbers = zeros(length(fmap_dirs), 1);
        for i = 1:length(fmap_dirs)
            parts = split(fmap_dirs(i).name, '-');
            scan_numbers(i) = str2double(parts{1});
        end
        [~, sort_idx] = sort(scan_numbers);
        fmap_dirs = fmap_dirs(sort_idx);

        for i = 1:length(fmap_dirs)
            dicom_path = fullfile(subject_path, 'scans', fmap_dirs(i).name, 'resources', 'DICOM', 'files');
            if ~exist(dicom_path, 'dir'), continue; end
            if i == 1
                out_dir = fullfile(output_base, 'fmap_mag');
            else
                out_dir = fullfile(output_base, 'fmap_phase');
            end
            cmd = sprintf('"%s" -z y -w 1 -f "fmap_%d" -o "%s" "%s"', ...
                DCM2NIIX_PATH, i, out_dir, dicom_path);
            [status, result] = system(cmd);
            if status ~= 0
                warning('dcm2niix failed for fieldmap %d, subject %d:\n%s', i, current_id, result);
            end
        end
    else
        fprintf('No fieldmaps found for subject %d — will use realign-only.\n', current_id);
    end

    % Anatomical (T1w)
    anat_dirs = dir(fullfile(subject_path, 'scans', '*-anat_T1w_acq_mprage'));
    if ~isempty(anat_dirs)
        dicom_path = fullfile(subject_path, 'scans', anat_dirs(1).name, 'resources', 'DICOM', 'files');
        if exist(dicom_path, 'dir')
            cmd = sprintf('"%s" -z y -w 1 -f "T1w" -o "%s" "%s"', ...
                DCM2NIIX_PATH, fullfile(output_base, 'anat'), dicom_path);
            system(cmd);
        end
    end

    % Gunzip all .nii.gz files produced by dcm2niix
    gz_files = dir(fullfile(output_base, '**', '*.nii.gz'));
    for i = 1:length(gz_files)
        gunzip(fullfile(gz_files(i).folder, gz_files(i).name), gz_files(i).folder);
    end


    %% ------------------------------------------------------------------
    %  PART 2: SLICE TIMING CORRECTION
    %  Interleaved acquisition
    % -------------------------------------------------------------------
    if do_slice_time
        clear matlabbatch;
        func_dir = fullfile(output_base, 'func');

        matlabbatch{1}.spm.temporal.st.scans    = {cellstr(spm_select('FPList', func_dir, '^run-.*\.nii$'))};
        matlabbatch{1}.spm.temporal.st.nslices  = param.st.nslices;
        matlabbatch{1}.spm.temporal.st.tr       = param.st.tr;
        matlabbatch{1}.spm.temporal.st.ta       = param.st.ta;
        matlabbatch{1}.spm.temporal.st.so       = param.st.so;
        matlabbatch{1}.spm.temporal.st.refslice = param.st.refslice;
        matlabbatch{1}.spm.temporal.st.prefix   = 'a';

        spm_jobman('run', matlabbatch);
        clear matlabbatch;
    end


    %% ------------------------------------------------------------------
    %  PART 3: FIELDMAP ESTIMATION (VDM file)
    %  Uses magnitude and phase fieldmap images to estimate a voxel
    %  displacement map (VDM)
    % -------------------------------------------------------------------
    if has_fieldmaps && do_fieldmap
        successVDM = createVDM_SPM12(current_id, BASE_PATH);
        if ~successVDM
            warning('VDM creation failed for subject %d. Falling back to realign-only.', current_id);
            has_fieldmaps = false;
        end
    end


    func_dir = fullfile(BASE_PATH, num2str(current_id), 'preprocessed', 'func');

    %% ------------------------------------------------------------------
    %  PART 4A: REALIGN & UNWARP (with fieldmap)
    % -------------------------------------------------------------------
    if has_fieldmaps && do_realign_unwarp
        clear matlabbatch;
        fmap_dir = fullfile(BASE_PATH, num2str(current_id), 'preprocessed', 'fmap_phase');
        vdm_file = cellstr(spm_select('FPList', fmap_dir, '^vdm5.*\.nii$'));

        for run = 1:nruns
            filter = sprintf('^arun-%d.nii$', run);
            matlabbatch{1}.spm.spatial.realignunwarp.data(run).scans  = cellstr(spm_select('FPList', func_dir, filter));
            matlabbatch{1}.spm.spatial.realignunwarp.data(run).pmscan = vdm_file;
        end

        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality  = 0.9;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep      = 2;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm     = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm      = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp  = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap    = [0 0 0];
        matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight   = '';

        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn    = [12 12];
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder  = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda    = 1e5;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm        = 0;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot       = [4 5];
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot       = [];
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm    = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem       = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi       = 5;
        matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround  = 'Average';

        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich   = [2 1];
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp   = 4;
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap      = [0 0 0];
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask      = 1;
        matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix    = 'u';

        spm_jobman('run', matlabbatch);
        clear matlabbatch;


    %% ------------------------------------------------------------------
    %  PART 4B: REALIGN ONLY (no fieldmap)
    %  Standard motion correction when no fieldmap is available
    % -------------------------------------------------------------------
    elseif do_realign_only
        fprintf('Realign-only (no fieldmap) for subject %d\n', current_id);
        clear matlabbatch;

        for run = 1:nruns
            filter = sprintf('^arun-%d.nii$', run);
            matlabbatch{1}.spm.spatial.realign.estwrite.data{run} = ...
                cellstr(spm_select('FPList', func_dir, filter));
        end

        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep     = 2;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm    = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm     = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp  = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap    = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight  = '';
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which   = [2 1];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp  = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap    = [0 0 0];
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask    = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix  = 'r';

        spm_jobman('run', matlabbatch);
        clear matlabbatch;
    end


    %% ------------------------------------------------------------------
    %  PART 5: SPATIAL SMOOTHING
    %  File prefix: s<kernel_size> (e.g. s3 for 3 mm FWHM)
    % -------------------------------------------------------------------
    if do_smooth
        img_files = cellstr(spm_select('FPList', func_dir, '^uarun-.*\.nii$'));
        if isempty(img_files) || isempty(img_files{1})
            error('No unwarped files (uarun-*) found in %s.\nIf you ran realign-only, update the filter to ^rarun-.* or ^arun-.*.', func_dir);
        end
        fprintf('Found %d volumes for smoothing.\n', length(img_files));

        for k = 1:length(param.smooth.kernels)
            clear matlabbatch;
            current_kernel = param.smooth.kernels(k);

            matlabbatch{1}.spm.spatial.smooth.data   = img_files;
            matlabbatch{1}.spm.spatial.smooth.fwhm   = repmat(current_kernel, 1, 3);
            matlabbatch{1}.spm.spatial.smooth.dtype  = 0;
            matlabbatch{1}.spm.spatial.smooth.im     = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = sprintf('s%d', current_kernel);

            fprintf('Smoothing at %d mm FWHM for subject %d...\n', current_kernel, current_id);
            spm_jobman('run', matlabbatch);
            clear matlabbatch;
        end

        fprintf('=== Done: subject %d ===\n', current_id);
    end

end 

end 
