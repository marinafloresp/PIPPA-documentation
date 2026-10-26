---
sidebar_position: 2
displayed_sidebar: documentationSidebar

---

# Software 

Several software are needed to run the pipeline. Some of these can be used from a conda. The set up script (0_setUpPipeline.sh) will create a conda environment where some of these packages will be installed by conda. As explained in *Requisites*, the condaEnv.yml file can be used to set up the conda environment. This can be done using either Anaconda or Miniconda, which local source directory needs to be specified in the *config.txt* file (`CONDA`). In order to download and install refer to: [Anaconda](https://docs.anaconda.com/anaconda/install/) or [Miniconda](https://docs.anaconda.com/miniconda/miniconda-install/)

Important software that will be installed in this environment are: MACS3 (3.0.2), SAMstats, SAMtools, BEDTools, Bowtie2, MultiQC, R, VerifyBamID, Python3, Pandoc, GATK, Picard. If you already have a conda environment, please specify its name or path in the *config.txt* file. 
  - Python version in this conda environment should be >= 3.12, as this is required by MACS3. For further details about MACS3 requirements: [MACS3 documentation](https://macs3-project.github.io/MACS/docs/INSTALL.html)
  - R version should be >= 4.3. Most packages are installed through the conda environment. A couple of packages would be installed when running the *0_setUpPipeline.sh* script. Some packages are no longer supported in the regular CRAN repository and need to be installed from source. To do this, download the source file of the package [ptest R package](https://cran.r-project.org/src/contrib/Archive/ptest/), keep this in the *config* folder and this will be installed from source when running the *0_setUpPipeline.sh* script.
