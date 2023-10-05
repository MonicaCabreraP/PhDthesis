# PCHi-C (Promoter Capture Hi-C)

PCHi-C combines the power of Hi-C with targeted enrichment of promoter regions, allowing for higher-resolution analysis interactions of each promoter in a cell type with other DNA elements of the genome, such as enhancers or other promoters.

In this thesis, to process the PCHi-C data I used the softwares [HiCUP](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5540598/pdf/pcbi.1005665.pdf]) for mapping, filter and calculate the capture efficiency, and [CHiCAGO](https://bioconductor.org/packages/release/bioc/html/Chicago.html) to perform the interaction calling.

A general workflow to process PCHi-C follows these steps:

1.	**Mapping and filtering**: [Here](https://stevenwingett.github.io/HiCUP/) you can find the documentation to run HiCUP. 

2.	**Capture Efficiency**: To perform this step there is a *Miscellaneous Functioality* inside HiCUP not shown in the documentation, here you will be able to find how to run it.
First, to use this script a bam file from HiCUP and a bed file with the coordinates of the captured restriction fragments are needed. For this purpose I generated a bed file with the coordinates of the captured fragments and also the annotation for each fragment, this extra information will be necessary for the other steps in further analysis (The bed file used in this thesis can be found in the data folder named: baits.bed). Then, execute:

```bash
perl HICUP/0.8.2/Misc/hicup_capture --baits baits.bed hicup.bam
```

3.	**Interaction calling**: In [this article](https://www.nature.com/articles/s41596-021-00567-5) you can find how to prepare the required files to perform this step and the commands. Once the files are prepare, run the runChicago.R script available in their [bitbucket repository](https://bitbucket.org/chicagoTeam/chicago/src). There are a lot of different arguments and you can personalize and make your our own script to run Chicago with your options.

4.	**Weight Recalibration of Interacting Calling**: Default parameters for interaction calling are calibrated on high-confidence calls from seven human Macrophage data sets (i.e. interactions that pass our p-value threshold in all seven samples). If the cell type of study is not too dissimilar to these calibration data, it should be fine to leave the parameters at their default settings. However, if your data set is from an unusual cell type, you may wish to recalibrate these parameters using data from cell types similar to yours. To do this, you should used the fitDistCurve.R script from ChicagoTools found in Chicago's bitbucket.

As in my thesis the datasets are from HCT116 cell line, a cancer cell line, I decided to recalibrate the parameters used to our control data as follows:

fitDistCurve.R --inputs 1stFile.Rda,2ndFile.Rda…

This procedure generates a file called ''cellType_summaryInput.Rda'' (you can find it in the data folder). This file was then used to re-run the CHiCAGO pipeline adding to the --settings-file the path where the cellType_summaryInput.Rda is.

5.	**Peak matrix**: Finally, CHiCAGO results from multiple experiments can be summarized in the form of a ‘peak matrix’ using makePeakMatrix.R from the chicago-Tools. The peakMatrix used in my thesis with all the PCHi-C samples can be found in the data folder.

Since PCHi-C specifically targets promoter regions, the analysis can be tailored to prioritize the identification and interpretation of promoter-enhancer interactions. The results of PCHi-C processing include two types of interactions:
  ●	**Promoter-promoter interactions**: interactions that involve physical contacts between promoters of different genes, which can provide insights into potential regulatory relationships, co-regulation, or functional coordination between genes.
  ●	**Promoter-other end interactions**: interactions between gene promoters and other genomic regions, such as enhancers, insulators, or other regulatory elements. These interactions reveal the spatial proximity and potential functional interactions between gene promoters and distal regulatory elements, playing a crucial role in gene regulation and transcriptional control.
