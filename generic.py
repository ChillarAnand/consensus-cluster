"""

(c) Michael Seiler 2008

Rutgers University
miseiler@gmail.com


"""

import numpy, sys
import pca, parsers, cluster, display

from mpi_compat import *

try:
    import psyco
    psyco.full()
except:
    pass

try:
    import euclidean
    EUC_C_EXT_ENABLED = 1
except:
    EUC_C_EXT_ENABLED = 0


#Methods

def makeplot(sdata, V, label):
    """

    Use matplotlib and display.py's Plot function to draw the samples along the first two Principle Components

    Usage: makeplot(sdata, V, label)
        
        sdata   - parsers.Parse object containing sample points and data in numpy.array form
        V       - The eigenvectors of the covariance matrix as determined by SVD
        label   - The filename will be of the form "label - timestamp.png"

    If matplotlib isn't installed, this function will simply do nothing.

    WARNING:    Depending on how the matrix is decomposed you may find different, but also correct, values of V
                This will manifest itself as the same plot, but reflected in one or both directions
    
    """

    M = numpy.array([ x.data for x in sdata.samples ])
    N = numpy.dot(V[:2], numpy.transpose(M))

    plots = [N]
    legend = ['Generic']

    display.Plot(plots, legend = legend, fig_label = label)

def run_pca(sdata, pca_fraction=0.85, eigenvector_weight=0.25):
    """

    Create a binary matrix via gen_matrix, normalise it, and then run PCA to reduce dimensionality.

    Usage: run_pca(sdata, pca_fraction, eigenvector_weight

        sdata               - parsers.Parse object with sample data as raw sequences
        pca_fraction        - The top pca_fraction fraction of principle components to keep
        eigenvector_weight  - The top fraction of SNPs to keep which occur with high weights in those principle components

    Returns: modified parsers.Parse object

    This function runs makeplot once the data in sdata has been converted to binary and then normalised.
    It calls console to log its results to screen and to logfile.

    """

    console = display.ConsoleDisplay(logname = 'PCA results')
    
    M = numpy.array([ x.data for x in sdata.samples ])

    console.log("Normalising %sx%s matrix" % (len(sdata.samples), len(sdata.samples[0].data)))

    M = pca.normalise(M, log2=False, sub_medians=False, center=True, scale=False)   #Only center the data

    #Unrolling pca.select_genes_by_pca...
    V = pca.pca(M, pca_fraction)    #From SVD
    SNP_indices = pca.select_genes(V, eigenvector_weight)

    console.log("Found %s principle components in the top %s fraction" % (len(V), pca_fraction)) #166
    console.log("Found %s reliable SNPs occurring with high weight (top %s by absolute value)" % (len(SNP_indices), eigenvector_weight)) #410

    #Don't reduce dimensionality right away, we need to take a picture
    for i in xrange(len(sdata.samples)):
        sdata.samples[i].data = M[i]
    
    makeplot(sdata, V, 'PCA results - All samples')

    #Reduce dimensions
    for i in xrange(len(sdata.samples)):
        sdata.samples[i].data = M[i].take(SNP_indices)

    return sdata

def clust_init(scale = 50):
    """

    Initialise clustering procedure and tell the user what's going on.

    Usage: clust_init(scale = 50)

        scale   - Number to randomly sample from each clade, using scale_by_clade

    Returns: modified sdata object

    Don't worry if MPI fails.  It's supposed to if you aren't using it.

    """

    console = display.ConsoleDisplay(log=False)
    
    console.write("Parsing data...")
    sdata = parsers.ParseRaw('testcase')
    console.success()
    
    console.write("Running PCA...")
    sdata = run_pca(sdata)
    console.success()
    
    console.write("Using MPI?")

    if MPI_ENABLED:
        console.success()
    else:
        console.fail()
    
    return sdata

def run_cluster(sdata, num_clusters = 6, subsamples = 500, subsample_fraction = 0.8):
    """

    Run the clustering routines, generate a heatmap of the consensus matrix, and fill the logs with cluster information.

    Each time this is run it will create a logfile with the number of clusters and subsamples in its name.  This contains
    information on which samples from which clade and haplogroup where clustered together for that particular K value.

    This function assumes clust_init has already been run.
    
    Usage: run_cluster(sdata, num_clusters, subsamples, subsample_fraction)

        sdata               - parsers.Parse object
        num_clusters        - K value, or the number of clusters for the clustering functions to find for each subsample.
        subsamples          - The number of subsampling iterations to run.  In each subsample, the SNPs, samples, or both could
                            be randomly selected for clustering.  This helps to ensure robust clustering.  More subsamples, more
                            robust clusters.
        subsample_fraction  - The fraction of SNPs, samples, or both to take each subsample.  0.8 is a good default.

    This function is designed to be run in iterations, so that one can call it for K = 2 to 20 unattended, for example.

    """
    
    console = display.ConsoleDisplay(logname = '%s clusters - %s subsamples' % (num_clusters, subsamples))
    
    console.log("\nSamples: %s" % len(sdata.samples))

    console.write("\nClustering data...")

    if EUC_C_EXT_ENABLED:
        euc_method = euclidean.euclidean
    else:
        euc_method = cluster.euclidean

    #Actual work
    clust_data = cluster.ConsensusCluster(sdata, distance_metric=euc_method, num_clusters=num_clusters, subsamples=subsamples, subsample_fraction=subsample_fraction)
    
    @only_once
    def report():

        clusters = dict()
    
        for clust_obj in [ (clust_data.datapoints[x].sample_id, clust_data.datapoints[x].cluster_id) for x in clust_data.reorder_indices ]:
            if clust_obj[1] not in clusters:
                clusters[clust_obj[1]] = [clust_obj[0]]
            else:
                clusters[clust_obj[1]].append(clust_obj[0])

        clust_count = 1
    
        for cluster in clusters:
            console.log("\nCluster %s:\n" % clust_count)
    
            for sample_id in clusters[cluster]:
                console.log("\t%s" % sample_id)
    
            clust_count += 1

    console.log("\nClustering results")
    console.log("---------------------")

    console.log("\nNumber of clusters: %s\nNumber of subsamples clustered: %s\nFraction of samples/genes used in subsample: %s" % (num_clusters, subsamples, subsample_fraction))
    console.log("\n---------------------")
    console.log("\nClusters")
    
    report()

    console.write("\n\nBuilding heatmap...")

    @only_once
    def save_hmap():
        filename = lambda s: "%s - %s clusters - %s subsamples" % (s, num_clusters, subsamples)

        display.Clustmap(clust_data).save(filename('Consensus Matrix'))

    save_hmap()

    if display.HMAP_ENABLED:
        console.success()
    else:
        console.fail()

    #If we're repeating, this is a good idea
    clust_data._reset_clusters()


if __name__ == '__main__':
    
    sdata = clust_init(scale = 20)

    for i in range(2, 7):
        run_cluster(sdata, num_clusters=i, subsamples = 300)
