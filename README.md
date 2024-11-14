# How to navigate this repository: #
There are three main folders- OV90, HCT116, and HEK293. These correspond to the cell type the experiment was performed in.
Within these folders there are sub folders divided by experiment type; PRO-seq, RNA-seq, ChIP-seq or MS-proteomics. See the experimental design figure below for more details.

Additional folders:
- metadata: incredibly helpful metadata tables that layout what samples are what (metadata files), as well as which comparisons were made in this study (interaction tables).
- scripts: general "scripts" folder that contains unpolished miscelanous scripts and notebooks related to this study.
- summary: contains files where multiple data types were combined to make summary plots and analyses, as well as gene target lists and other misc data used in this study.

### Within OV90: ###
#### PRO-seq ####
Contains the data associated with **GSE228346**.
This data is PRO-seq data generated from OV90 cells treated +/- 50nM SY-5609 (CDK7 inhibitor) for 30min, then was heat shocked at 42C for 0min, 20min or 45min.

#### RNA-seq ####
There are two subsets of data within this folder:
1) HS dataset contains the data associated with **GSE228346**.
This data is RNA-seq data generated from OV90 cells treated +/- SY-5609 for 30min, then was heat shocked at 42C for 0min, 30min or 120min. The 30min dataset has a 60min recovery before the mRNA was harvestd via polyA selection.

2) Long exposure SY-5609 contains the data associated with **GSE262536**.
This data is a RNA-seq data set prepared separately from the HS data set. In this experimental design OV90 cells were treated with 0nM, 40nM or 400nM SY-5609 for 4hr before mRNA harvest.

#### Proteomics ####
This data was generated from OV90 cells treated +/- SY-5609 for 30min or 120min, then the peptides were isolated for MS analysis.
The associated file is prot_res, it contains the log2FC values and p-values associated with DMSO vs SY-5609 30min and DMSO vs SY-5609 120min.
The raw data is also uploaded in the form of an excel file.

### Within HCT116: ###
#### ChIP-seq ####
Contains the data associated with **GSE268531**.
This data is ChIP-seq data generated from HCT116 cells treated +/- SY-5609 for 30min, then was treated with IFNg for 45min.

#### PRO-seq ####
Contains the data associated with **GSE261575**.
This data is PRO-seq data generated from HCT116 cells treated +/- SY-5609 for 30min, then was treated with IFNg for 45min.
There is also data +/- Palbociclib (a CDK4/6 inhibitor) with a 30min treatment.

#### Proteomics ####
This data was generated from HCT116 cells treated +/- SY-5609 for 30min or 120min, then the peptides were isolated for MS analysis.
The associated file is prot_res, it contains the log2FC values and p-values associated with DMSO vs SY-5609 30min and DMSO vs SY-5609 120min.
The raw data is also uploaded in the form of an excel file.

![gitTREATMENTscheme](https://github.com/Dowell-Lab/CDK7_inhibition/assets/48491008/b2f2ba94-8fc8-4a1c-9e66-852c1c99e971)

### Within HEK293: ###
#### MNase-ChIP-seq ####
Contains the data associated with **GSE218269**.
This data is MNase-ChIP-seq data generated from CDK7as cells treated +/- 1-NM-PP1 for 30min.

# Main programs used: #
**General programs:**
- DEseq2 (v1.26.0, R v3.6; https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
 
Counts were performed with Rsubread featurecounts (v2.0.1; https://www.rdocumentation.org/packages/Rsubread/versions/1.22.2/topics/featureCounts).

Associated files: deseq2_res; contains ROI, baseMean, sample counts (batch corrected and normalized), multiple comparisons and their associated log2FC and adjusted p-values.

- GSEA (v4.3.2, c5.all.v7.5.symbols.gmt; https://www.gsea-msigdb.org/gsea/index.jsp)

Associated files: gsea_res; contains the PATHWAY, GO-database type from C5 all (biological processes, reactome, etc), pathway size, multiple comparisons and their associated normalized enrichment scores (NES) and FDR q-values.

**RNA-seq specific programs:**
- TFEA.ChIP (v1.6.0; https://www.bioconductor.org/packages/release/bioc/vignettes/TFEA.ChIP/inst/doc/TFEA.ChIP.html)

Associated files: tfea.chip; contains enrichment results from TFEA.ChIP, which consist of enrichment scores (ES) and adjusted p-values.

**PRO-seq specific programs:**
- TFEA (v1.1.4; https://github.com/Dowell-Lab/TFEA)

Associated files: tfea_res; contains the HOCOMOCO TF name,  multiple comparisons and their associated number of motif events, corrected enrichment scores and corrected adjusted p-values. The correction is the correction for GC bias in TF motifs. Note that a conversion table for HOCOMOCO TF name and geneIDs is available in summary files.

- TF Profiler (v1.0; https://github.com/Dowell-Lab/TF_profiler)

Associated files: tf_profile; contains the HOCOMOCO TF name, the number of hits used to calculate the MD-score (small hits/large hits; hits within 300bp/hits within 3000bp of mu) for both experimental data and the statistical model, the MD-score for both experimental data and the statistical model and the associated significance values.

In order to run TFEA and TF Profiler bidirectional calls are necessary. In order to annotate bidirectional regions we used Bidirectional-Flow to run tfit (v1.2;https://github.com/Dowell-Lab/Tfit; https://github.com/Dowell-Lab/Bidirectional-Flow) on each sample. To make a consensus set of bidirectionals we used mumerge (v1.1.0; https://github.com/Dowell-Lab/mumerge). This consensus sets (one per cell line) were used for all bidirectional analyses in this study and is uploaded to the OV90 and HCT116 PROseq folders.
