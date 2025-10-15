# Diapause_temporal_modeling

A model-based temporal analysis tailored to time-resolved transcriptomics.

This approach enabled robust modelling of gene-expression trajectories across temporal samples, 
avoiding local inaccuracies that might arise in the absence of replicates. 

## Part 1: RNA-Seq modeling

Folder: rnaseq_analysis

This model uses a hierarchy of three nested parametric functions: 
1. constant (null),
2. sigmoid (a single sustained change)
3. “impulse” (a transient change), which was shown to successfully capture transcriptional dynamics (Chechik and Koller, 2009).

For each gene, a likelihood ratio test selects the best-fitting function of these three options,
identifying transcripts with significant temporal changes during diapause exit.

## Part 2: 2-condition RNA-Seq modeling

Folder: rnaseq_romney

This model expands the temporal modeling framework to compare gene expression between two conditions,
and identify condition-specific expression differences, arising either from distinct dynamics or levels.
