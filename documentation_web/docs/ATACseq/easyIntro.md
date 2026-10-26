---
sidebar_position: 1
displayed_sidebar: documentationSidebar

---

# PIPPA overview

Pipeline for Integrative Peak calling in Population-scale ATAC-seq (**PIPPA**) provides end-to-end quality control and processing of ATAC-seq datasets and curated of reproducible peak sets. It is designed to run quality-control (QC) of samples in 3 stages, followed by a final stage of peak set curation. 

The quality-control part of the pipeline is structured in 3 main stages:
1. ATAC-sequencing data: confirms quality of sequencing data and ATAC-seq specific properties.
2. Sample identity: confirms samples belong to the assigned individual.
3. Cell type identity: confirms sample is representative of the assigned cell type.

After QC has been completed for all samples, these are then used to generated curated peak sets for each cell type. This process is structured in three steps:
1. Default peak calling to generate peak set.
2. Segmented analysis is performed on default peak set to find optimal peak calling threshold.
3. Peak calling is repeated using the chosen stringent threshold.

You can find a more detailed explanation of PIPPA [here](/img/pipelineOverview.pdf)