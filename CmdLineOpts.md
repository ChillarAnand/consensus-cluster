# General Options #

**-f** _filename_
> Load _filename_ for clustering.  Default Parser: Normal
**-p** _parser_
> Parse _filename_ with _parser_.  Only valid with the **-f** option.  E.g. 'Normal'
**-d**
> > Don't init display, run from console. This happens automatically if there is no display or the required libraries are unavailable.
**-h**, **--help**

> Print help and exit.

# Data Normalisation #

**--log2**, **--nolog2**
> Perform log2 reexpression, or turn it off. Default is off.
**--submedians**, **--nosubmedians**
> Perform median centring, or turn it off. Default is off.
> NOTE: Turning this on will turn off mean centring.
**--center**, **--nocenter**
> Perform mean centring, or turn it off. Default is on.
> NOTE: Turning this on will turn off median centring.
**--scale**, **--noscale**
> Perform RMS scaling, or turn it off. Default is off.
**--normvar**, **--nonormvar**
> Normalise variance to 1 for each subsample, or turn it off. Default is off.

# PCA and Feature Selection #

**--nopca**
> Do not perform PCA at all. This precludes feature selection. Useful if your data is known to be singular.
**--pcafraction** _fraction_
> Select features from the top _fraction_ principle components. Default is 0.85
**--eigweight** _fraction_
> Select the top _fraction_ features by weight in each principle component. Default is 0.25
**--noselection**
> Do not perform feature selection. Simply sets pcafraction and eigweight to 1.

# Sample Selection #

**-c** _file1 file2 file3 .._
> Define samples (one on each line) in file1, etc as clusters.  Sample set will be reduced to these samples, and their labels will be shown in logs and PCA plot.
**--krange** _fst_ _snd_
> Repeat for each kvalue between _fst_ and _snd_ inclusive. Default is 2 to 6.
**--subsamples**  _number_
> Number of clustering iterations to perform. Default is 300.
**--subfraction** _fraction_
> Select a random _fraction_ of the samples each iteration. Default is 0.80

# Clustering Options #

**--kmeans**
> Run the K-Means algorithm
**--som**
> Run the Self-Organising Map algorithm
**--pam**
> Run the Partition Around Medoids algorithm
**--hier**
> Run the Hierarchical Clustering algorithm. Note that this option adds the Hierarchical algorithm to clustering iterations, rather than the 'final' consensus clustering.
**--euclid**
> Cluster using the Euclidean distance metric. Default.
**--corr**
> Cluster using the Pearson Correlation distance metric
**--coring**
> Turns on EXPERIMENTAL coring support. Additional logfiles and images are generated which detail suggested 'core' clusters. Take its advice at your own risk!

# Example #

> python common.py -f mydata.txt -d --kmeans --log2 --submedians --noselection -c clusterdefs/`*`

> Opens mydata.txt, log2 reexpresses and median centres the data, performs no feature selection, and begin k-means clustering using the cluster definitions in the clusterdefs folder without using the GUI.