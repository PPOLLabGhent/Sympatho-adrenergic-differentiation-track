## Single cell RNA sequencing guide for sympatho-adrenergic-differentiation-track

Citation: Stéphane Van Haver, Yujie Fan, Sarah-Lee Bekaert, Celine Everaert, Wouter Van Loocke, Vittorio Zanzani, Joke Deschildre, Vanessa Vermeirssen, Katleen De Preter, Ting Zhou, Alex Kentsis, Lorenz Studer, Frank Speleman, Stephen S. Roberts (2023). Human iPSC modeling recapitulates in vivo sympathoadrenal development and reveals an aberrant developmental subpopulation in familial neuroblastoma. IScience, https://www.sciencedirect.com/science/article/pii/S2589004223021739

Link to publication
https://www.sciencedirect.com/science/article/pii/S2589004223021739

**Preparing input files**

Prior to single-cell RNA-Seq (scRNAseq) the cells of selected time-points were hashed and pooled together for scRNAseq analysis. The time-points D25-26, D27-28, D29-30 and D31-32 were hashed for the WT cells. For the ALK mutant cells time-points D25-27-29-32 were hashed. The scRNAseqq was performed on Chromium instrument (10X genomics) following the user guide manual (Reagent Kit 3’ v3.1). The indexed DNA libraries were double-size purified (0.6–0.8X) with SPRI beads and sequenced on Illumina NovaSeq S4 platform (R1 – 26 cycles, i7 – 8 cycles, R2 – 70 cycles or higher).

**Running Cell Ranger 10XGenomics**

The sequenced libraries were quantified with Cell Ranger software (version 6.0.0) (Zheng et al., 2017) with default parameters. The standard reference provided with Cell Ranger (version 2020-A) was used. This version was build on a GRCh38 human genome (bundle 3.0.0) and a filtered version of Gencode v32. In addition, HTO surface antibodies were quantified with the Cell Ranger count feature reference option. The feature reference files are provided on GitHub. Each sample has it own library and feature file.

Once the input files are ready the scripts *CellrangerCount_FASTQ.sh* and *CellrangerCount_Libraries.sh* can be executed for each sample separately on the bash shell prompt as follows:

$ bash CellrangerCount_FASTQ.sh

$ bash CellrangerCount_Libraries.sh

**Create SeuratObjects**

The count files were loaded in R and further analysis with Seurat (version 4.1.1) (Hao et al., 2021) was performed. Cells with less than 300 or more than 4000 detected genes or with more than 10% mitochondrial genes were omitted. Subsequently, demultiplexing of the mixes was executed by HTODemux function and only singlet were retained. The data was normalized and scaled with the standard Seurat functions.

The public data (GSE147821) of Adameyko (2021) was integrated with the own samples by making use of the Seurat integration features.
This was done for the ALK (*SeuratObject_ALK_Data.R*) and the noALK (*SeuratObject_NoALK_Data.R*) samples separately.

**Cells of interest**

The focus of this research is the differentiation track of the sympatho adrenergic cells, therefore cells of interest were subjected to further in-depth analysis (SCP, Bridging, CPC, SAP and proliferating SAP). They were extracted, re-embedded and re-clustered, followed by a second post-clustering quality control phase (*StéphaneOnlyData_NoALK.Rmd*).

* SCP: schwann cell precursors cells
* CPC: chromaffin-like connecting progenitor
* SAP: sympatho-adrenergic progenitors

Marker genes that defined clusters by differential expression were identified using the Seurat FindAllMarkers function. Clusters were annotated to cell types by comparison of marker genes for each cluster to canonical cell type markers from the literature.

**Differentiation between ALK and noALK data**

Using the ALK(R1275Q) mutant iPSCs, we performed serial scRNAseq at four time points between days 25 and 32 following the same differentiation protocol. Next, we integrated the scRNAseq data from these AL(KR1275Q) mutant iPSCs with the previously generated WT iPSC-data for a comparative UMAP analysis. This revealed a large overlap between the ALK-mutant and WT populations, suggesting that the majority of ALK-mutant cells resemble their WT counterparts and undergo normal differentiation. The script describes the differential gene expression analysis between populations and between cell types (*DifferentialGeneExpression_ALK_vs_noALK.Rmd*).

**Diffusion and Pseudo Time Point Analysis**

The destiny R package v.3.10.0 (Angerer et al., 2016) was used to calculate the diffusion map embeddings. Single-cell pseudotime trajectories were constructed using the slingshot R package v2.4.0 (Street et al., 2018) based on clusters identified by Seurat and the DC components from destiny. A node in cluster 12 (SCPs) was selected as the starting point for the trajectory (*Diffusion_PseudoTimePlot_Analysis.Rmd*).
