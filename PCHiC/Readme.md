# PCHi-C (Promoter Capture Hi-C)

PCHi-C combines the power of Hi-C with targeted enrichment of promoter regions, allowing for higher-resolution analysis interactions of each promoter in a cell type with other DNA elements of the genome, such as enhancers or other promoters.

In this thesis, to process the PCHi-C data I used the softwares [HiCUP](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5540598/pdf/pcbi.1005665.pdf]) for mapping, filter and calculate the capture efficiency, and [CHiCAGO](https://bioconductor.org/packages/release/bioc/html/Chicago.html) to perform the interaction calling.

A general workflow to process PCHi-C follows these steps:

1.	**Mapping and filtering**

3.	**Capture Efficiency**
4.	**Interaction calling**
  
The analysis of PCHi-C data follows similar principles to standard Hi-C data analysis with further steps (they will be mentioned in the methods section). Since PCHi-C specifically targets promoter regions, the analysis can be tailored to prioritize the identification and interpretation of promoter-enhancer interactions. The results can include two types of interactions:
●	**Promoter-promoter interactions**: interactions that involve physical contacts between promoters of different genes, which can provide insights into potential regulatory relationships, co-regulation, or functional coordination between genes.
●	**Promoter-other end interactions**: interactions between gene promoters and other genomic regions, such as enhancers, insulators, or other regulatory elements. These interactions reveal the spatial proximity and potential functional interactions between gene promoters and distal regulatory elements, playing a crucial role in gene regulation and transcriptional control.
