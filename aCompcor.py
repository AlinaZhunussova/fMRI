import os
from nipype.algorithms.confounds import ACompCor

IDs = [
    "sub-01", "sub-02", "sub-03", "sub-04", "sub-05",
    "sub-06", "sub-07", "sub-08", "sub-09", "sub-10",
    # ...
]

parpath = "/mypath/data_drive/"
outdir = os.path.join(parpath, "AcompCor_regressors_fmriprep_style")
os.makedirs(outdir, exist_ok=True)

successful_runs = 0
failed_runs = []

for subject_id in IDs:
    print(f"\nProcessing subject {subject_id}")
    
    func_dir = f'/mypath/preprocessed_multisession/{subject_id}/func/'
    mask_dir = os.path.join(parpath, f'compcor_masks_in_func_space_fmriprepv_rigidt/{subject_id}/')

    csf_mask = os.path.join(mask_dir, 'c3T1w_func_masked_bin.nii.gz')        # CSF mask
    wm_mask = os.path.join(mask_dir, 'acompcor_wm_func_masked_bin.nii.gz')   # WM mask
    combined_mask = os.path.join(mask_dir, 'acompcor_wmcsf_func_masked_bin.nii.gz')  # Combined WM+CSF
    
    os.makedirs(os.path.join(outdir, subject_id), exist_ok=True)
    
    # Process 6 runs
    prefix = "uarun"
    
    for run_num in range(1, 7):
        fmripath = os.path.join(func_dir, f"{prefix}-{run_num}.nii")
        
        if not os.path.exists(fmripath):
            print(f"  Run {run_num}: File not found")
            continue
        
        try:
            print(f"  Processing Run {run_num}")
            
            ccinterface = ACompCor()
            ccinterface.inputs.realigned_file = fmripath
            ccinterface.inputs.mask_files = [csf_mask, wm_mask, combined_mask]
            ccinterface.inputs.merge_method= 'none'
            ccinterface.inputs.num_components = 6  #change to a different value or also could be variance-based like ccinterface.inputs.variance_threshold = 0.5
            ccinterface.inputs.pre_filter = False # polynomial is dealing with nonlinear, low-frequency artifacts in fMRI time series
            ccinterface.inputs.repetition_time = 1.8
            ccinterface.inputs.components_file = os.path.join(outdir, subject_id, f"{subject_id}_run-{run_num:02d}_compcor.txt")
            ccinterface.run()
            
            successful_runs += 1
            print(f"     Success")
            
        except Exception as e:
            print(f"     ERROR: {str(e)}")
            failed_runs.append((subject_id, f"run-{run_num:02d}", str(e)))
            
print(f"\n{'='*70}")
print(f"SUMMARY: {successful_runs} successful, {len(failed_runs)} failed")
if failed_runs:
    print(f"\nFAILED RUNS:")
    for subject, run, reason in failed_runs:
        print(f"   {subject} {run}: {reason}") 