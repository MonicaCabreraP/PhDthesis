# Hi-C (High-throughput Chromosome Conformation Capture)

Hi-C provides valuable insights into how chromosomes are spatially arranged and interact with each other. 

In this thesis, to process the Hi-C data I used the [TADbit](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5540598/pdf/pcbi.1005665.pdf) pipeline with a specific version available at https://github.com/fransua/TADbit/tree/p53_javierre.

In the [TADbit tutorial website](https://3dgenomes.github.io/TADbit/tutorial.html#from-fastq-files-to-interaction-matrices) or in the methods section of my thesis you can find the commands to run the software, but a general workflow to process Hi-C follows these steps:

1.	**Quality control** 
2.	**Mapping** 
3.	**Read filtering** 
4.	**Hi-C data normalization**
5.	**Interaction matrix generation**: For this thesis interaction matrices were generated at two resolutions, 100 kb for A/B compartments and 50 kb for TADs.
  
    ●	*A/B compartment analysis*: A/B compartments were determined independently for each chromosome and time-point at a resolution of 100 kb. The [ICE-normalized](https://pubmed.ncbi.nlm.nih.gov/22941365/) interaction matrices were distance corrected and transformed into [Pearson correlation matrices](https://pubmed.ncbi.nlm.nih.gov/22941365/). The A compartment was assigned to genomic bins with a positive first principal component (PC1), while the B compartment was assigned to genomic bins with negative PC1. The resulting transformed eigenvectors were referred to as compartment scores.
  	
    ●	*TAD identification*: TADs were identified at a resolution of 50 kb t. TAD borders were assigned scores ranging from 1 to 10 based on their robustness in the TAD border detection. The strength of a TAD border was calculated by the number of times it appeared in the optimal pathway. TAD borders were also identified using the [Insulation score](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4498965/). 

After running all the steps of TADbit, you should have one folder per sample with these folders inside:

![image](https://github.com/MonicaCabreraP/PhDthesis/assets/31327141/0499c5c3-6094-4606-b06f-ffe1feb1fbfd)

The 06_segmentation folder will contain each chromosome's eigenvector values for each bin. Merge in a single file all the chromosomes eigenvector values in different ros for each sample in different columns to obtain a single file. In my case my file is named: ``aligned_Compartments_TADbit.tsv``

Then run the Compartments.R script to clean the file obtained and prepare it for visualization.

