# Cell Atlas and Novel Signal Processing Strategy in Primate Insular Cortex
<img src="Abstract_n.jpg" alt="图片描述" width="宽度" height="高度">
Rui-Feng Liu*, Mengyao Huang*, Yuhui Shen*, Mingting Shao*, Junzhan Jing, Nana Xu, Lei Tang, Biaodi Liu, Jianming Shi,  Fanrui Chen, Zhao-Zhe Hao, Xiaolong Jiang and Sheng Liu

#### This repository contains the analysis code and the preprocessed data for the above manuscript.
-----------------------------------------------
### Repository Structure

The repository is organized as follows:

&emsp;&emsp;● **10X**: Contains Jupyter notebooks for the analysis of scRNA-seq data (10X Genomics) and the generation of related figures.

&emsp;&emsp;● **Metabolism**: Contains the normalized metabolic data table and codes used to generate related figures (Figure 8 and Extended Data Figure 10).

&emsp;&emsp;● **Modeling**: Contains hoc scripts for neuron models (Neuron 8.0) and MATLAB codes for generating related figures.

&emsp;&emsp;● **Ephys**: Contains Python and MATLAB codes for electrophysiological data analysis and figure generation.

&emsp;&emsp;● **Morphology**: Contains Python and MATLAB codes for morphological data analysis and figure generation.

&emsp;&emsp;● **patchSeq**: Contains codes for analyzing the transcriptomic data of Patch-seq cells.

&emsp;&emsp;● **patchSeqMappingTo10X** Mapping to 10X: Contains codes for mapping Patch-seq neurons transcriptomically to the scRNA-seq neuronal atlas.

--------------------------------------
### Preprocessed data and meta data

Because of the large file sizes, some processed data files are not included directly in this repository. However, all processed data files have been deposited at **ZENODO** and are organized in the same hierarchical folder structure as presented here (Github). They are available for download from https://doi.org/10.5281/zenodo.18814764. The processed data (h5ad format) from the AIC scRNA-seq and Patch-seq experiments, along with the corresponding raw data (FASTQ format), have also been deposited in the NCBI GEO database (see following).

------------------------------------------------
### Reproduce the analysis figures

After downloading the processed data from Zenodo, place it into the corresponding directories within the GitHub repository's data folders. Then, execute the provided Python notebooks or MATLAB scripts to regenerate the figures.

We constructed the computational models using NEURON (version 8.0). The corresponding hoc scripts are available in the <code>'modelling/'</code> folder, with specific experimental models organized in relevant sub-folders. MATLAB scripts for plotting the figures are located in <code> 'modeling/FigPlot/' </code>. Simulated neuronal traces data are also deposited in **ZENODO**.

-------------------------------------------------
### Raw data

All primary datasets from this study have been deposited in public repositories.

&emsp;&emsp;● **Transcriptomics**: Raw transcriptomic data in FASTQ format, meta data and relavent processed data (.h5ad), including both scRNA-seq (10x) and Patch-seq data, are available at NCBI GEO under accession codes <code>GSExxxx1</code>  and <code>GSExxxx2</code>, respectively. They can be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM319557 and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE319369.

&emsp;&emsp;● **Electrophysiology**: Raw electrophysiological data in NWB format, along with corresponding metadata, can be accessed from the DANDI Archive at: 

&emsp;&emsp;&emsp;https://dandiarchive.org/dandiset/001746 (electrophysiology for Patch-seq) 

&emsp;&emsp;&emsp;https://dandiarchive.org/dandiset/001750 (Sodium chanel current) 

&emsp;&emsp;&emsp;https://dandiarchive.org/dandiset/001751 (Simutaneous multi-channel Patch-clamp recording) 

&emsp;&emsp;&emsp;https://dandiarchive.org/dandiset/001752 (Post-synaptic current dataset)

&emsp;&emsp;&emsp; These datasets include electrophysiology traces of firing patterns, outside-out sodium channel currents, multi-patch configurations, and postsynaptic current traces.

&emsp;&emsp;● **Morphology**: Raw morphological reconstructions are deposited at NeuroMorpho.Org and are available via the persistent link: https://doi.org/10.13021/y4be-0p18.

&emsp;&emsp;● **Metabolomics**: Raw metabolic data are available at MetaboLights under accession code https://www.ebi.ac.uk/metabolights/editor/MTBLS13927.

------------------------------------------
### Level of Support

We will update this code as necessary.







