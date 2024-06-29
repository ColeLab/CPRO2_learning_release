# CPRO2_learning_release

Public code release for Mill & Cole 2023: "Neural representation dynamics reveal computational principles of cognitive task learning"  

Author: Ravi Mill (rdm146@rutgers.edu)

Includes the Matlab (v2019a) function quantify_rule_type_subject.m, which quantifies the strength of compositional and conjunctive neural representations in block-to-block fMRI regional activations evoked by the C-PRO2 paradigm. The remaining files included with the repository are required for the function to run (see details below).

Top of the main function provides a detailed description of the function, as well as input and output variables. 

DEPENDENCIES:
1. Connectome workbench must be installed and accessible to your system path, including the command line version (wb_command), https://www.humanconnectome.org/software/get-connectome-workbench
2. (Included with the repo) Relevant files from Matlab Gifti toolbox https://github.com/gllmflndn/gifti, allowing for opening (ciftiopen) and saving (ciftisave and ciftisavereset) HCP CIFTI files
3. (Included with the repo) CAB-NP network partition .dlabel file (CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii) from https://github.com/ColeLab/ColeAnticevicNetPartition that links vertex/voxel data to affiliated Glasser atlas regions
4. Example data to test code execution: i) Behavioral data for one subject is provided in the data/ subdirectory, ii) fMRI vertexwise task activation data for one subject, which also needs to be stored in the data/ subdirectory, but needs to be downloaded first from this access link https://rutgers.box.com/s/8abw0n2ydeau2688rdozozzrqgq5hf37.

Note that running the function on the supplied example data for one subject on the Rutgers Amarel HPC cluster (2x Intel Xeon Gold 6338 Processors (48MB cache, 2.0GHz), 1 CPU requested, 80GB memory) takes ~114 seconds for execution. 
