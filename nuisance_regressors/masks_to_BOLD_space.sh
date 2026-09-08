#According to fmriprep pipeline, where they do inverse transform and Gaussian interpolation,
#I bring BOLD to T1w using antsRegistrationSyN.sh -d 3 -t r, and then coregister masks to BOLD space using antsApplyTransforms -n Gaussian  

PARTICIPANTS=(sub-01 sub-02 sub-03 sub-04 sub-05 sub-06 sub-07 sub-08 sub-09 sub-10
              sub-11 sub-12 sub-13 sub-14 sub-15 sub-16 sub-17 sub-18 sub-19 sub-20
              sub-21 sub-22 sub-23 sub-24 sub-25 sub-26 sub-27 sub-28 sub-29 sub-30
              sub-31 sub-32 sub-33 sub-34 sub-35 sub-36 sub-37 sub-38 sub-39 sub-40
              sub-41 sub-42 sub-43 sub-44 sub-45 sub-46 sub-47 sub-48 sub-49 sub-50)

ID=${PARTICIPANTS[$SLURM_ARRAY_TASK_ID]}

compcor_path="/mypath/compcor_masks_after_fmriprep/${ID}"
mean_func_path="/mypath/states_fmri_data/${ID}/preprocessed_multisession/func"
anat_path="/mypath/states_fmri_data/${ID}/preprocessed_multisession/anat"
output_dir="/mypath/compcor_masks_in_func_space_fmriprepv_rigidt/${ID}"

mkdir -p "${output_dir}"

func_ref="${mean_func_path}/meanuarun-1_stripped.nii"
t1_brain="${anat_path}/T1w_bias_corrected_stripped.nii"

echo "Checking input files..."
missing=0
for f in "${func_ref}" "${t1_brain}" \
         "${compcor_path}/acompcor_wm.nii.gz" \
         "${compcor_path}/acompcor_wmcsf.nii.gz" \
         "${compcor_path}/c3T1w.nii"; do
    if [[ ! -f "$f" ]]; then
        echo "Missing file: $f"
        missing=1
    fi
done

if [ ${missing} -eq 1 ]; then
    echo "  Skipping ${ID} due to missing files"
    exit 1
fi

echo "All input files found"

# Register functional to T1w space (like fmriprep: boldref2anat)

antsRegistrationSyN.sh -d 3 -t r \
    -m "${func_ref}" \
    -f "${t1_brain}" \
    -o "${output_dir}/reg_SSmeanfunc2SST1w_"

affine_file="${output_dir}/reg_SSmeanfunc2SST1w_0GenericAffine.mat"

if [[ ! -f "${affine_file}" ]]; then
    echo "Registration failed for subject ${ID}"
    exit 1
fi

echo "Registration completed"

# Transform masks from T1w to functional space

antsApplyTransforms \
    -d 3 -n Gaussian -v 0 \
    -i "${compcor_path}/acompcor_wm.nii.gz" \
    -r "${func_ref}" \
    -t [${affine_file}, 1] \
    -o "${output_dir}/acompcor_wm_func.nii.gz"

antsApplyTransforms \
    -d 3 -n Gaussian -v 0 \
    -i "${compcor_path}/acompcor_wmcsf.nii.gz" \
    -r "${func_ref}" \
    -t [${affine_file}, 1] \
    -o "${output_dir}/acompcor_wmcsf_func.nii.gz"

antsApplyTransforms \
    -d 3 -n Gaussian -v 0 \
    -i "${compcor_path}/c3T1w.nii" \
    -r "${func_ref}" \
    -t [${affine_file}, 1] \
    -o "${output_dir}/c3T1w_func.nii.gz"

echo "All transformations completed for subject ${ID}"
