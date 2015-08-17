# How do I figure out how my clusters changed from run to run? #

It turns out ConsensusCluster has a mechanism for dealing with this problem, and that is to allow the user to "name" samples in clusters by adding them to a text file.

For example, say in one clustering run we have the following clusters after one run:

| **Cluster 1** | **Cluster 2** |
|:--------------|:--------------|
| twas          | did           |
| brillig       | gyre          |
| and           | gimble        |
| the           | in            |
| slithy        | wabe          |
| toves         | mimsy         |

and we want to discover their whereabouts in the next clustering run. To do this, we do the following:
  * Create two text files, cluster1.txt and cluster2.txt
  * Put in each file the contents of the cluster by sample id, one sample id per line
  * Click define clusters and **select BOTH files** (hold ctrl to select multiple files) or use the -c option on the command line, e.g. -c cluster`*`.txt
  * Begin clustering

### IMPORTANT NOTE ###
Each time you define clusters within the GUI, the previous ones are overwritten! Selecting multiple files in a single window is mandatory for using multiple cluster definitions!

ConsensusCluster will automatically label your samples in the PCA plot:

![http://consensus-cluster.googlecode.com/files/PCAexample.png](http://consensus-cluster.googlecode.com/files/PCAexample.png)

And tell you where each sample came from in the logfiles:

```
Cluster 0 (blue):

	toves		cluster1.txt

Cluster 1 (green):

	brillig		cluster1.txt
	gyre		cluster2.txt
	in		cluster2.txt
	did		cluster2.txt
	wabe		cluster2.txt
	slithy		cluster1.txt
	and		cluster1.txt
	mimsy		cluster2.txt
	the		cluster1.txt
	gimble		cluster2.txt
```

### IMPORTANT NOTE ###

When you define clusters, ConsensusCluster will ONLY use these sample ids you have defined in the clustering! This is also a useful way to cluster only a subset of your data.