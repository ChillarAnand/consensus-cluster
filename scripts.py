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

    Returns: modified sdata object

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

    f = open(remfile, 'r')
    rem = [x.strip() for x in f]
    f.close()

    f = open(filen, 'r')
    sams = [x.strip() for x in f if x.strip() not in rem]
    f.close()

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
