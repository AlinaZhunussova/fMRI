# According to fmriprep pipeline, where they do acc_msk_brain = pe.MapNode(ApplyMask(), name='acc_msk_brain', iterfield=['in_file']) 
# I mask my masks in BOLD space with the mean functional brain mask using fslmaths -mas:use (following image > 0) to mask the current image.
# According to fmriprep pipeline, where they do acc_msk_bin = pe.MapNode(Binarize(thresh_low=0.99), name='acc_msk_bin', iterfield=['in_file']) 
# I am binarizing the mask and setting the threshold to 0.99 using fslmaths with the -thr 0.99 -bin options.

PARTICIPANTS=(sub-01 sub-02 sub-03 sub-04 sub-05
              sub-06 sub-07 sub-08 sub-09 sub-10
              sub-11 sub-12 sub-13 sub-14 sub-15
              sub-16 sub-17 sub-18 sub-19 sub-20
              sub-21 sub-22 sub-23 sub-24 sub-25
              sub-26 sub-27 sub-28 sub-29 sub-30
              sub-31 sub-32 sub-33 sub-34 sub-35
              sub-36 sub-37 sub-38 sub-39 sub-40
              sub-41 sub-42 sub-43 sub-44 sub-45
              sub-46 sub-47 sub-48 sub-49 sub-50)

ID=${PARTICIPANTS[$SLURM_ARRAY_TASK_ID]}

mask_dir="/mypath/compcor_masks_in_func_space_fmriprepv_rigidt/${ID}"
func_path="/mypath/states_fmri_data/${ID}/preprocessed_multisession/func"
output_dir="${mask_dir}"

brain_mask="${func_path}/meanuarun-1_brain_mask.nii.gz"

missing=0
if [[ ! -f "${brain_mask}" ]]; then
    echo "Missing brain mask: ${brain_mask}"
    missing=1
fi

for mask in "${mask_dir}/acompcor_wm_func.nii.gz" \
            "${mask_dir}/acompcor_wmcsf_func.nii.gz" \
            "${mask_dir}/c3T1w_func.nii.gz"; do
    if [[ ! -f "${mask}" ]]; then
        echo "Missing file: ${mask}"
        missing=1
    fi
done

if [ ${missing} -eq 1 ]; then
    echo "  Skipping ${ID} due to missing files"
    exit 1
fi

# Mask WM
fslmaths "${mask_dir}/acompcor_wm_func.nii.gz" \
    -mas "${brain_mask}" \
    "${output_dir}/acompcor_wm_func_masked.nii.gz"

# Mask WM+CSF
fslmaths "${mask_dir}/acompcor_wmcsf_func.nii.gz" \
    -mas "${brain_mask}" \
    "${output_dir}/acompcor_wmcsf_func_masked.nii.gz"

# Mask c3T1w
fslmaths "${mask_dir}/c3T1w_func.nii.gz" \
    -mas "${brain_mask}" \
    "${output_dir}/c3T1w_func_masked.nii.gz"

# Binarize WM
fslmaths "${output_dir}/acompcor_wm_func_masked.nii.gz" \
    -thr 0.99 -bin \
    "${output_dir}/acompcor_wm_func_masked_bin.nii.gz"

# Binarize WM+CSF
fslmaths "${output_dir}/acompcor_wmcsf_func_masked.nii.gz" \
    -thr 0.99 -bin \
    "${output_dir}/acompcor_wmcsf_func_masked_bin.nii.gz"

# Binarize c3T1w
fslmaths "${output_dir}/c3T1w_func_masked.nii.gz" \
    -thr 0.99 -bin \
    "${output_dir}/c3T1w_func_masked_bin.nii.gz"

echo "All masks are masked and binarized"