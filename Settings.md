# The Settings Tab #

## General ##

  * K-Value Range: Repeat the clustering procedure using each K-Value found in this range.  The K-Value is passed to each clustering algorithm to indicate how many clusters it should "find" in the data.
  * Subsamples: Number of times to create a random dataset and cluster it.  ConsensusCluster attempts to avoid sample bias by randomly subsampling from your data, clustering each subsample, then forming a consensus from all of these clustering iterations.  Larger values will increase accuracy, with diminishing returns.
  * Fraction to Sample: The percentage of the genes/samples/both to subsample each clustering iteration.  See Subsamples, above.

## PCA ##

  * Log2: Log2 re-express each data value.
  * Sub Medians: Subtract the median of medians from each data value.
  * Center: Center the data by subtracting the mean from each data value.
  * Scale: Divide each data value by the root-mean-square.

  * PCA Fraction: Keep the top PCA Fraction percentage of Principle Components in the data.
  * Eigenvalue Weight: Within those principle component eigenvectors determined by PCA Fraction (see above), use only the top features with eigenvalues which occur with this weight by absolute value.  I.e., 0.25 would take the top 25% of features in this eigenvectors.

## Algorithm ##

  * K-Means: K-Means algorithm.  Recommended.
  * SOM: Self-Organising Map Algorithm.  Lengthy.  Recommended.
  * PAM: Partition Around Medoids Algorithm.  Not recommended unless your data would benefit from medoid cluster representation.
  * Hierarchical: Hierarchical clustering.  Will attempt to perform clustering with Single, Complete, and/or Average linkage depending on those settings.  Not recommended except for _exceptionally_ clean data.
  * Cluster Consensus Using: Once a consensus matrix is made, this is clustered to deliver the "final" clusters.  It is HIGHLY recommended that Hierarchical is used, because this will allow ConsensusCluster to produce a dendrogram of its clustering.

  * Distance Metric: Euclidean and Pearson Correlation metrics are provided.

## Misc ##

  * Set Variance to 1: Standardise the variance of each feature to be 1 over all samples each clustering iteration.  This is NOT recommended for microarray data, though other data sources with features that can dominate clustering due to differing scales may benefit from this (try "scale" first). Can be used for crude combination of multiple datasets, though virtually any other method is better.