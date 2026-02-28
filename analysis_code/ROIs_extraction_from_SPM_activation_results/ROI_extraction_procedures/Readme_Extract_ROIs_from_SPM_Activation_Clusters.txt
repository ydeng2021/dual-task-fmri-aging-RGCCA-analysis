README: Extract ROIs from SPM Activation Results

After obtaining the second-level results for each contrast in SPM, 
you can export a mask of your region(s) of interest (ROI) for further brain connectivity analysis.

Step 1: Export a Cluster Mask from SPM

1. Open the SPM Results window and click Results.
2. Select your desired contrast and enter the threshold settings.
3. In the Results -> Display window, click Save.
4. You will see several saving options, such as:
save...
thresholded SPM
all clusters (binary)
all clusters (n-ary)
current clusters
5. Choose “all clusters (binary)”, then enter an output filename, such as 
"Single_Motor_all_clusters". 
SPM will automatically save the file as "Single_Motor_all_clusters.nii" in the current working directory.

This binary file masks significant clusters are with 1 (inside the cluster) and non-significant voxels with 0 (outside the cluster).

Step 2: Split SPM clusters into ROIs using the Harvard-Oxford atlas

Run the MATLAB script:"split_SPM_activation_clusters_into_ROI_HO_atlas.m"
Setup
1. Edit the script to specify your path on line 8 and line 9:
data_path = 'your_data_directory';
atlas_path = 'your_atlas_directory';
2. Make sure the following three atlas files (from the CONN toolbox) are included in your atlas_path:
atlas.txt
atlas.nii
short_lables.txt
3. On line 14, specify the mask output file created in step 1, for example:
maskfilename = 'Single_Motor_all_clusters.nii'; 
Run the script in MATLAB.
It will automatically split the significant clusters into ROIs aligned with the Harvard–Oxford atlas space.
Each ROI will be saved as an individual file for further analysis.

Step 3: Import the Extracted ROIs into the CONN Toolbox

The extracted ROI files can be directly imported into the CONN toolbox for ROI-to-ROI functional connectivity analysis, 
as proposed in the manuscript.
Open your CONN project.
Go to Setup → ROIs.
Click New ROI and select the ROI NIfTI files generated from Step 2.
Once imported, these ROIs can be used to define source and target regions for ROI-to-ROI connectivity analyses.
Continue with the First-level and Second-level analyses in CONN as usual.