
WaspsProject
========

This project contains the supplementary materials for the paper 
[Transcriptomics of an extended phenotype: parasite manipulation of wasp social behavior correlates with shifting expression of caste-related genes]()

Abstract
========
We investigated the relationship between the primitively eusocial paper wasp, Polistes dominula, and its obligate endoparasite, Xenos vesparum (Insecta: Strepsiptera) to understand the effects of a parasite on its host’s brain transcriptome. Previous research suggests that X. vesparum exerts control over physio-behavioral toolboxes involved in host caste determination in ways that benefit the parasite. Specifically, we hypothesized that X. vesparum-infested (stylopized) females that would typically be workers show a shift in their transcriptomic profiles towards being more like pre-overwintering queens (gynes), as the behavior of stylopized females resembles that of typical gynes. We used RNA-sequencing data to determine patterns of brain gene expression in stylopized females, and compared this to unstylopized workers and gynes. In support of our hypothesis, we found that stylopized females, despite sharing numerous physiological and life history characteristics with members of the worker caste, show similar brain expression patterns to those of unstylopized gynes. These data suggest the parasite affects its host by utilizing preexisting host neurogenetic machinery to shift naturally-occurring social behavior in a way that is beneficial to the parasite, rather than inducing completely novel host behavior. 

DE Analysis script using guide
==================

1. Use `git clone https://github.com/ruolin/WaspsProject/` to download or use the download link in the upper-right corner on the project page. 
2. Unzip the input files in the folder HTSeq-Count-BFAST-raw-reads-count
3. Before Running the R script, change the working dir to your local dir by uncommenting and change this line `# setwd("/path/to/your/current/dir")`
The R script also generates two excel spreadsheets, one for the condition-averaged raw reads counts for the 8484 expressed genes, the other for the edges R FDR values for each pairwise comparison.  
