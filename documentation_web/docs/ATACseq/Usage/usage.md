---
sidebar_position: 3
displayed_sidebar: documentationSidebar

---

# Configuration

In order to use PIPPA, two main configuration files need to be set up. For examples of these config files, go to ATACSeq/config
1. config.txt: This is the main configuration file that will source the pipeline.
  - This file should be found in the main directory of the project/dataset to be analysed.
  - The main variables to be changed in this configuration file are:
    - `PROJECT`: name of the project's folder.
    - `MAIN_DIR`: full path to the project's folder.
    - `REFERENCES_DIR`: full path to directory with all references used throughout the pipeline.
    - Specify the modules versions to be loaded and the full paths to the softwares used in the pipeline (e.g. picard).
  - An alphabetically ordered array of cell types to which samples belong to: `CELLTYPES`. 
  - Cutoff values for peak calling should be specified by the user. If not specified, default will be used (0.05). 
    - `macs_sample_thres` : Sample-level peak calling q-value cutoff
    - `macs_group_thres` : Cell fraction-level peak calling q-value cutoff
    - `macs_chr_thres` : Peak calling on sex chromosomes q-value cutoff
2. config.r: This is the configuration file for running R scripts.
  - The `mainDir` variable should be changed to be the full path to the project's directory. This should match `MAIN_DIR` variable in the config.txt file.
  - Other parameters or threshold should be changed as required.
  - An alphabetically ordered list of cell types to which samples belong to needs to be specified as *cellTypes*.

Additional source files of R packages that are installed manually will also be located in this folder. A temporary directory can be either specified in the system's _.bashrc_ file or in the *config.txt* as an additional variable `TMPDIR`. Please make sure you have full permissions for the directory specified. For each main bash script, make sure the SBATCH details are filled in/removed if not apply.