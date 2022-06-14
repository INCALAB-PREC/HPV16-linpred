# HPV16-linpred 

Human Papillomavirus (HPV) is the causal agent of 5% of cancers worldwide and the
main cause of cervical cancer and it is also associated with a significant percentage of
oropharyngeal and anogenital cancers. More than 60% of cervical cancers are caused by
HPV16 genotype, which has been classified into lineages (A, B, C, and D). Lineages are
related to the progression of cervical cancer and the current method to assess lineages
is by building a Maximum Likelihood Tree (MLT); which is slow, it cannot assess poor
sequenced samples, and annotation is done manually. In this study, we have developed
a new model to assess HPV16 lineage using machine learning tools. A total of 645
HPV16 genomes were analyzed using Genome-Wide Association Study (GWAS), which
identified 56 lineage-specific Single Nucleotide Polymorphisms (SNPs). From the SNPs
found, training-test models were constructed using different algorithms such as Random
Forest (RF), Support Vector Machine (SVM), and K-nearest neighbor (KNN). A distinct set
of HPV16 sequences (n = 1,028), whose lineage was previously determined by MLT, was
used for validation. The RF-based model allowed a precise assignment of HPV16 lineage,
showing an accuracy of 99.5% in the known lineage samples. Moreover, the RF model
could assess lineage to 273 samples that MLT could not determine. In terms of computer
consuming time, the RF-based model was almost 40 times faster than MLT. Having a fast
and efficient method for assigning HPV16 lineages, could facilitate the implementation of
lineage classification as a triage or prognostic marker in the clinical setting.

**REFERENCES**
Asensio-Puig L, Alemany L and Pavón MA 2022 A Straightforward HPV16 Lineage Classification Based on Machine Learning, Front. Artif. Intell. 5:851841. doi: 10.3389/frai.2022.851841
