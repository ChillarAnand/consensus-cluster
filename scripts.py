#Random utilities

import parsers, numpy, os, log_analyse, pca


def union(list1, list2):
    """
    Remove the indices which make up the union between two lists in log n time
    
    Returns a tuple, where tuple[0] is the list of indices in list1 which is in common with list2, and tuple[1]
    is the same list for list2

    """

    swapped = False
    if len(list1) > len(list2):         #Make list2 the longer one
        list1, list2 = list2, list1
        swapped = True

    indices_list1 = numpy.argsort(list1)
    indices_list2 = numpy.argsort(list2)

    union_indices_list1 = []
    union_indices_list2 = []
    
    breakpoint = 0

    for i in indices_list1:    
        for j in range(len(indices_list2))[breakpoint:]:    #Ugly, but reduces complexity
            idx = indices_list2[j]

            if list1[i] == list2[idx]:
                union_indices_list1.append(i)
                union_indices_list2.append(idx)
                breakpoint = j
                break

    if not swapped:
        return union_indices_list1, union_indices_list2

    return union_indices_list2, union_indices_list1

def scale_to_set(sdata, filenames):
    """

    scale_to_set(filename)

        Removes all but those sample_ids you specifiy.

        filenames    - array of filenames, each file containing list of sample ids to use

    Returns: modified sdata object, dict of cluster->indices

    """

    defined_clusters = {}

    for filename in filenames:
        name = os.path.basename(filename)
        defined_clusters[name] = parsers.get_list_from_file(filename)

    samples_to_keep = sum([ defined_clusters[x] for x in defined_clusters ], [])

    sample_id_list = [ x.sample_id for x in sdata.samples ]
    
    sample_indices = union(sample_id_list, samples_to_keep)[0]
    sdata.samples = [ sdata.samples[x] for x in sample_indices ] #Adjustment

    sample_id_list = [ x.sample_id for x in sdata.samples ] #This is different!
    
    for name in defined_clusters: #If samples aren't in the main, ignore them
        sample_list = defined_clusters[name]
        def_indices = union(sample_list, sample_id_list)[0]
        defined_clusters[name] = [ sample_list[x] for x in def_indices ]

    return sdata, defined_clusters

def scale_probes(sdata, filename):
    """

    scale_probes(sdata, filename)
        
        Removes all gene probes except those you specify

        filename    - File containing a list of probes, one on each line

    Returns: modified sdata object

    NOTE: Currently there is no way to call this from common.py!  This will change in the future.
          For now, stick this in your _preprocess(self) subclass.

    """

    probes_to_keep = union(sdata.gene_names, parsers.get_list_from_file(filename))[0]

    sdata.gene_names = sdata.gene_names.take(tuple(probes_to_keep))

    for i in xrange(len(sdata.samples)):
        sdata.samples[i].data = sdata.samples[i].data.take(tuple(probes_to_keep))

    return sdata

def crop_sample_list(filen, remfile):
    #filen: sample list to cropped
    #remfile: samples you'd like to remove from filen

    rem = parsers.get_list_from_file(remfile)
    sams = numpy.array(parsers.get_list_from_file(filen)) #Take is amazing.

    ind_to_remove = numpy.array(union(sams, rem)[0])
    ind_to_keep = numpy.lib.arraysetops.setdiff1d(numpy.arange(len(sams)), ind_to_remove)

    sams = sams.take(tuple(ind_to_keep))

    f = open(filen, 'w')
    for sam in sams:
        f.write(sam)
        f.write("\n")
    f.close()

def new_defined_clusters(sdata, conv):
    #sdata: sample data obj
    #conv: conversion dict, keys sample ids values new cluster assignments
    #Stick this in your preprocess function (see common.py for subclassing help)

    new_clusts = {}
    s_ids = [x.sample_id for x in sdata.samples]

    for s_id in s_ids:
        if s_id in conv:
            new_clusts.setdefault(conv[s_id], []).append(s_id)
        else:
            new_clusts.setdefault('Unknown', []).append(s_id)

    return new_clusts

def write_normal(sdata, filename):
    #Takes an sdata obj and writes out a tab-delimited datafile, suitable for ParseNormal
    #useful to convert between data formats
    
    if not len(sdata.gene_names) > 0: 
        raise ValueError, "No gene names found! Unsuitable for this data format."

    sids = [x.sample_id for x in sdata.samples]

    f = open(filename, 'w')

    #Sample line, first row
    f.write("\t".join(['SAMPLE ID'] + sids))
    f.write("\n")

    #Data
    for i in xrange(len(sdata.gene_names)):
        f.write("\t".join([sdata.gene_names[i]] + [ str(x.data[i]) for x in sdata.samples ]))
        f.write("\n")

    f.close()

def make_def_clusters_from_log(logfile):
    """Takes a logfile and writes cluster definition files"""

    logdict = log_analyse.gen_cluster_dict(logfile)

    for clustname in logdict:
        name = clustname.split() #Most currently: # (colour)
        filen = 'cluster_' + str(name[0]) #cluster_0, etc

        f = open(filen, 'w')

        for sample in logdict[clustname]:
            f.write(sample)
            f.write("\n")

        f.close()

def remove_pc(sdata, num=1):
    """Remove the first num principle components from the data"""

    M = numpy.array([x.data for x in sdata.samples]) #Sample in rows, genes in columns
    M = pca.normalise(M, log2=False, center=True, scale=False, sub_medians=False)

    u, s, V = numpy.linalg.svd(M, 0)        #Decompose
    S = numpy.identity(s.shape[0]) * s

    for i in xrange(num):
        S[i][i] = 0.        #Sets the offending eigenvalue to 0

    M = numpy.dot(numpy.dot(u, S), V)       #Recompose

    for i in xrange(len(sdata.samples)):
        sdata.samples[i].data = M[i]

    return sdata

def write_table(ndict, filename):
    """Write a tab delimited flat file, one key per line"""

    ls = ndict.keys()
    ls.sort()

    f = open(filename, 'w')
    
    for key in ls:
        f.write("\t".join([key, ndict[key]]))
        f.write("\n")

    f.close()

def write_ratio(s, clust1, clust2, filename, threshold=1.):
    """Write SNR ratios given sdata obj, clust1 filename, clust2 filename, file to write to, threshold for SNR"""

    M = numpy.array([x.data for x in s.samples])

    f = open(filename, 'w')
    
    c1ind = get_indices(s, clust1)
    c2ind = get_indices(s, clust2)
    
    ratios = pca.snr(M, c1ind, c2ind, threshold)
    
    f.write("%s vs %s:\n" % (clust1, clust2))
    f.write("--------------------\n")
    f.write("Gene ID\t\t%s Avg\t%s Avg\tSNR Ratio\n" % (clust1, clust2))
    
    for result in ratios:
        f.write("\n%10s\t\t%f\t\t%f\t\t%f" % (s.gene_names[result[1]], result[2], result[3], result[0]))

    f.close()
                        
def write_classifier(s, filename, clust1, clust2=None, threshold=None):
    """
    Writes a bayesian binary classifier
    
    See pca.binary_classifier for details
    
    s - SampleData obj
    filename - What to write the classification information to
    clust1 - a file with one sample name on each line which composes the cluster you're trying to define
    clust2 - an optional second cluster.  Otherwise it's every other sample.
    threshold - If you want to classify using only genes over a certain SNR threshold, use this.  None uses all.

    """

    M = numpy.array([x.data for x in s.samples])

    f = open(filename, 'w')
    
    c1ind = get_indices(s, clust1)
    
    if clust2 is not None:
        c2ind = get_indices(s, clust2)
    else:
        c2ind = numpy.lib.arraysetops.setdiff1d(numpy.arange(len(s.samples)), numpy.array(c1ind))

    rlist, w0 = pca.binary_classifier(M, c1ind, c2ind, threshold)

    f.write("%s vs %s:\n" % (clust1, clust2))
    f.write("--------------------\n\n")

    #Returns (a, b), where a is w in (wi, i) pairs and b is w0
    f.write("w0 is %s\n" % w0)
    f.write("\nGene ID\t\tMultiplier\n\n")

    rlist.sort() #FIXME: abs?

    for result in rlist:
        f.write("%10s\t%f\n" % (s.gene_names[result[1]], result[0]))

    f.close()

def get_indices(s, filename):
    """
    Return the indices of the samples in filename in the sdata object
    
    """

    sams = parsers.get_list_from_file(filename)
    return union([ x.sample_id for x in s.samples ], sams)[0]
