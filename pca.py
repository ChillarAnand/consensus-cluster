"""

PCA and other normalisation methods.


Copyright 2008 Michael Seiler
Rutgers University
miseiler@gmail.com

This file is part of ConsensusCluster.

ConsensusCluster is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ConsensusCluster is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ConsensusCluster.  If not, see <http://www.gnu.org/licenses/>.


"""

import numpy

def select_genes_by_pca(M, pca_fraction=0.85, eigenvector_weight=0.15):
    """
    Expects an array with samples on the rows and genes on the columns

    pca_fraction        - Fraction of eigenvalues which explain pca_fraction of the variance to accept
    eigenvector_weight  - The top eigenvector_weight (by absolute value) fraction of those genes which occur with high weights
                          in those eigenvectors which correspond to the eigenvalues explained by pca_fraction
    
    Returns a tuple of the indices of those genes which comprise the top weight% in each eigenvector

    """

    V = pca(M, pca_fraction)

    return select_genes(V, eigenvector_weight)

def pca(M, frac):
    """Takes a matrix M and returns those eigenvectors which explain frac of the variance"""

    u, s, v = numpy.linalg.svd(M, 0) #Will run out of memory from U otherwise
    
    variances = s**2/M.shape[1]
    total_variances = numpy.sum(variances, 0)

    variance_fractions = numpy.divide(variances, total_variances)

    for i in range(1, len(variance_fractions) + 1):
        if numpy.sum(variance_fractions[:i], 0) >= frac:
            break

    if i < 2:
        i = 2   #Minimum 2

    #return numpy.transpose(numpy.dot(v[:i], numpy.transpose(M))) #The transformed data
    return v[:i]
    
def select_genes(v, weight):
    """Returns a tuple of the indices of those genes which comprise the top weight% in each eigenvector"""

    genes = [0] * len(v[0])
    gene_indices = []

    for vec in v:
        min_value = (1 - weight) * numpy.min(vec)
        max_value = (1 - weight) * numpy.max(vec)
 
        for i in xrange(len(vec)):
            if vec[i] <= min_value or vec[i] >= max_value:
                genes[i] = 1

    for i in xrange(len(genes)):
        if genes[i]:
            gene_indices.append(i)

    if len(gene_indices) < 2:
        raise TypeError, "Not enough genes at %s%% weight in eigenvectors" % (weight)

    return tuple(gene_indices)

def normalise(M, log2=True, sub_medians=True, center=True, scale=True):
    """
    Perform a number of normalisation routines on M
    
    log2                - log2 transform the data (Yes if this is raw gene data)
    sub_medians         - subtract the median of medians from each data value
    center              - subtract the data by its average, making the overall mean 0
    scale               - subtract the root-mean-square of the data AFTER centering

    """

    if log2:
        M = log2_transform(M)

    if sub_medians:
        M = subtract_medians(M)

    if scale or center:
        M = center_matrix(M)

    if scale:
        M = scale_matrix(M)

    return M

def center_matrix(M):
    """Subtract the mean from matrix M, resulting in a matrix with mean 0"""

    return (M - numpy.average(M, 0))

def scale_matrix(M):
    """Subtract the root-mean-square from each data member"""

    numpy.seterr(all='raise')   #Catch divide-by-zero, otherwise SVD won't converge

    T = numpy.transpose(M)
    for i in range(M.shape[1]):
        if T[i].any():
            T[i] = T[i] / numpy.sqrt(numpy.sum(T[i]**2) / (M.shape[0] - 1))

    return M

def get_median(lst, n):
    """Gets the median of a list or array"""
    
    if isinstance(lst, numpy.ndarray):
        c = lst.copy()
    else:
        c = lst[:]

    c.sort()

    if n & 1:
        return c[n // 2]

    return (c[n // 2] + c[n // 2 - 1]) / 2.

def log2_transform(M):
    """Take the log2 of each value in M"""

    return numpy.log(M)/numpy.log(2)

def subtract_medians(M):
    """Subtract each value in M by the median of medians"""

    medians = []

    for sample_row in M:
        medians.append(get_median(sample_row, len(sample_row)))

    return (M - get_median(medians, len(medians)))
    
def get_columns(M, list1, list2):
    """

    Return a list of (col1, col2) pairs

        list1 - list of row indices for the first group
        list2 - list of row indices for the second group

    This function assumes samples are in rows and genes are in columns

    """

    cols = []

    N1 = M.take(tuple(list1), 0)
    N2 = M.take(tuple(list2), 0)

    for i in xrange(len(M[0])):
        col1 = N1[:,i]
        col2 = N2[:,i]
    
        cols.append((col1, col2))

    return cols

def snr(M, list1, list2, threshold = None):
    """

    Performs a signal-to-noise ratio test on M, assuming samples are in rows and genes are in columns

        list1   - List of row indices for first group
        list2   - List of row indices for second group
        threshold - Minimum SNR ratio to report

    Returns a reverse-ordered list of (ratio, index, mean1, mean2) pairs, where index is the column index of the gene,
    and mean1 and mean2 correspond to the mean for that particular gene in list1 and list2, respectively

    """

    ratios = []
    cols = get_columns(M, list1, list2)

    for i in xrange(len(cols)):

        col1, col2 = cols[i][0], cols[i][1]
        mean1, mean2 = col1.mean(), col2.mean()

        ratios.append( (abs((mean1 - mean2) / (col1.std() + col2.std())), i, mean1, mean2) )

    ratios.sort(reverse=True)

    if threshold is not None:
        for i in xrange(len(ratios)):
            if ratios[i][0] < threshold:
                break
    else:
        i = len(ratios)

    return ratios[:i]

def binary_classifier(M, list1, list2, threshold = None, prior_prob = None):
    """

    Create a bayesian linear-discrimination function based on two clusters being list1 and list2

        list1   - List of row indices for first group
        list2   - List of row indices for second group
        threshold - Minimum SNR ratio to report

    WARNING: This function makes some assumptions of which you should be aware:
        
        Genes are independent events (likely false)
        Your prior probabilities are indeed indicative of a true population
        Your data is indicative of a true population
        Your clusters are definitively known, and perfect
        Up/down status determined by gene average IN YOUR DATA SET truly represents average and upgregulation/downregulation (definitely false!)

    NOTE: Prior probabilities unimplemented

    """
    
    w = []
    w0 = 0.0

    ran = numpy.random.random

    cols = get_columns(M, list1, list2)

    snr_genes = []
    for ratio in snr(M, list1, list2, threshold):
        snr_genes.append(ratio[1])

    for i in snr_genes:
        p_array = cols[i][0]
        q_array = cols[i][1]

        p = (p_array > 0.0).sum() / float(len(p_array))  #array([True, False, False, True... etc
        q = (q_array > 0.0).sum() / float(len(q_array))

        #hack
        if p > 0.999:
            p = 0.999
        if q > 0.999:
            q = 0.999
        if p < 0.001:
            p = 0.001
        if q < 0.001:
            q = 0.001

        w.append( (numpy.log( (p * (1 - q)) / (q * (1 - p)) ), i) ) #(wi, i) pairs
        w0 += numpy.log( (1 - p) / (1 - q) )

    return w, w0 #g(x) = w*xi for x in w + w0
