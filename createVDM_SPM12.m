function successVDM = createVDM_SPM12(current_id, base_path)
% Creates a Voxel Displacement Map (VDM) for fieldmap-based distortion
% correction using SPM12's FieldMap toolbox.
%
% The VDM encodes the B0 field inhomogeneity at each voxel and is used
% by Realign & Unwarp (Part 4A of the preprocessing pipeline) to correct
% susceptibility-induced geometric distortions in EPI images.
%
% INPUTS:
%   current_id  - numeric subject ID (e.g. 101)
%   base_path   - root data directory (same BASE_PATH as the main script)
%
% OUTPUT:
%   successVDM  - logical; true if VDM file was created successfully
%
% DEPENDENCIES:
%   SPM12 with FieldMap toolbox
%
% -------------------------------------------------------------------------
% USER SETTINGS — replace placeholder paths with your own
% -------------------------------------------------------------------------

SPM_PATH      = '/path/to/spm12';                          % SPM12 root
FIELDMAP_PATH = fullfile(SPM_PATH, 'toolbox', 'FieldMap'); % FieldMap toolbox
FIELDMAP_T1   = fullfile(FIELDMAP_PATH, 'T1.nii');         % template for brain masking
% BEST comment everis

% -------------------------------------------------------------------------

addpath(SPM_PATH);
addpath(FIELDMAP_PATH);
spm('defaults', 'fmri');
spm_jobman('initcfg');

% Directory layout created by the main preprocessing script
preproc_dir    = fullfile(base_path, num2str(current_id), 'preprocessed');
fmap_phase_dir = fullfile(preproc_dir, 'fmap_phase');
fmap_mag_dir   = fullfile(preproc_dir, 'fmap_mag');
func_dir       = fullfile(preproc_dir, 'func');

% Input images
%   Phase image  : difference-of-phase map (second fieldmap echo, _ph suffix)
%   Magnitude    : magnitude image at shorter TE (better signal)
%   EPI reference: first volume of run-1, used to register VDM to functional space
phase_img = cellstr(spm_select('FPList', fmap_phase_dir, '^fmap_2_e2_ph\.nii$'));
mag_img   = cellstr(spm_select('FPList', fmap_mag_dir,   '^fmap_1_e1\.nii$'));
epi_ref   = cellstr(spm_select('FPList', func_dir,       '^arun-1.*\.nii$'));
% ^ uses run-1 only; regex '^arun-1\.' avoids matching run-10 etc.

if isempty(phase_img{1}) || isempty(mag_img{1}) || isempty(epi_ref{1})
    warning('createVDM: missing input files for subject %d. Check directory contents.', current_id);
    successVDM = false;
    return;
end

matlabbatch = {};

% Phase & magnitude images
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase     = phase_img;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = mag_img;

matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et        = [4.92 7.38]; % echo times (ms): short TE, long TE
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 1;            % skull-strip fieldmap
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir   = -1;           % phase-encode blip direction
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert      = 19.97;        % total EPI readout time (ms)
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm     = 0;            % fieldmap NOT acquired with EPI
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm       = 0;            % do not apply Jacobian modulation

% Phase unwrapping options (Mark3D method)
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm   = 10;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad    = 0;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws     = 1;

% Brain masking options
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {FIELDMAP_T1};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm     = 5;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode   = 2;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate  = 4;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh   = 0.5;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg      = 0.02;

% EPI reference for VDM registration
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session.epi  = epi_ref;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm     = 1;  % match VDM to EPI space
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname     = 'session';
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0; % don't write unwarped EPI
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat         = '';
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat    = 0;

spm_jobman('run', matlabbatch);
clear matlabbatch;

% Check that the VDM file was actually created
vdm_file = dir(fullfile(fmap_phase_dir, 'vdm5_*.nii'));
if isempty(vdm_file)
    warning('createVDM: VDM file not found after running FieldMap for subject %d.', current_id);
    successVDM = false;
else
    fprintf('VDM created successfully for subject %d: %s\n', current_id, vdm_file(1).name);
    successVDM = true;
end

end
