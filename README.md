# deltabit

Welcome to Deltabit! Deltabit is tool designed for detecting clusters of horizontally transferred genes in genomes. It's still a work in progress, but check back for updates! Deltabit was originally designed for the purpose of detecting _Starships_, giant horizontally-transferring transposons in fungi, though it isn't limited to this use case. Deltabit is based off the existing metric for detecting horizontal gene transfer called Alien Index. Deltabit uses the genomics software suite Mycotools.

# Summary

A user provides a genome of interest and a focal gene (referred to as a captain) in that genome. For all genes in a specified distance from the captain, Deltabit BLASTs it against a designated set of "Ingroup" proteins from closely related species and "Outgroup" proteins from distantly related species. Deltabit looks at the best hit from each of these BLASTs and takes the difference in the bitscores (Best Ingroup Bitscore - Best Outgroup Bitscore) and assigns this value as the Deltabit score for that gene. Deltabit then plots each gene's log2Deltabit and distance from the captain (with some smoothing to make trends easier to see and pick out). Deltabit will also display the best Ingroup and Outgroup bitscores for each gene for reference. Clusters of genes with low Deltabit scores will appear as peaks in these plots and are indicative of potential horizontal transfer with the Outgroup.

Deltabit can also run a pipeline referred to as Deltabit Deep that involves the standard Deltabit analysis along with some additional downstream analyses. In Deltabit Deep, Deltabit scores are also calculated for BUSCO genes in the genome which are used as a null distribution for determining statisical significance of any peaks in the plot. Deltabit Deep will identify a stastically significant peak near the captain and use BLAST to search for similar clusters in the Ingroup and Outgroup genomes and produce an interactive synteny plot. 

For more detailed information on Deltabit and the Deltabit Deep pipeline, please refer to the documentation. 

# Authors

Deltabit has been developed by Conor Bendett at the University of Wisconsin-Madison with significant contributions from Dr. Samuel O'Donnell and Dr. Emile Gluck-Thaler. A publication detailing Deltabit is in the works. 
