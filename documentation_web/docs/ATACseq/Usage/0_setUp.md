---
sidebar_position: 4
displayed_sidebar: documentationSidebar

---

# 0. Set up pipeline

This is an optional step of the pipeline. It will check that required files, packages or libraries that are as expected by the pipeline.

`sbatch ./jobSubmission/0_setUp.sh (project directory) `

##### -scripts executed-
- ./Rscripts/installLibraries.r : checks for already installed R libraries and install if not found.

##### -parameters-
- `(project-name)` project's directory.

##### -requires-
- condaEnv.yml: list of conda packages to be installed in conda.
- config.r: R config file with data paths and thresholds used in analysis.
- _ptest_ R library is currently not available to install through CRAN, so the source file needs to be downloaded and made available from the config directory.