This software package allows the user to perform various consensus clustering operations on any numerical data using PCA to select informative features.

### News ###

**NEW** - The power of consensus clustering with the speed of the GPU: The first alpha of [CUDAConsensusCluster](http://code.google.com/p/cuda-consensus-cluster) has been released!

Defining clusters explained! Find out how to track your clusters between runs in DefiningClusters.

Clustering can now be fully automated! See command line options in the CmdLineOpts Wiki.

### Current Version - 0.6 ###

You can find the most current version for your operating system in the "Featured" box on the left of this page.

### Installation ###

If using the Win32 binary, simply extract the zipped folder.  The executable is called 'common.exe'

If using the source package, please have the following packages installed:

_Required, all OS_
  * Python 2.6 (3.0 NOT currently supported)
  * Python-numpy >= 1.1.0
  * Python-scipy

_Highly Recommended, all OS_
  * Python-matplotlib >= 0.98.0
  * Python-pil >= 1.1.6
  * Python-pygtk

_Linux Optional_
  * psyco
  * GCC
  * GNU Make
  * pyMPI

Please see the README file for operating system instructions for installation of the C-Source modules, which are HIGHLY RECOMMENDED for purposes of increased speed.  There are also instructions for the use of multiprocessing using MPI on linux systems.

### Quick Start Guide ###

Create a tab-delimited file for your data.  Put the 'samples' to be clustered in the columns, and the 'features' that define them in the rows.  Make sure the first row and the first column are labels.

_Example_

|SAMPLE\_ID|Sample 1|Sample 2|...|
|:---------|:-------|:-------|:--|
|Gene 1    |3.4234  |1.4244  |
|...       |

To use the GUI, run common.exe or type python common.py from a command line if you downloaded the Win32 binary or the python source, respectively.

Choose File->Open and select your datafile.  Make sure Normal is the selected parser (this is the parser which understands the file format described in the example above.)

Click Begin Clustering.  Clustering progress will be tracked in the text window and output in the form of PNG and LOG files.

Configuration option descriptions are found in the Settings Wiki: http://code.google.com/p/consensus-cluster/wiki/Settings

### Citing ConsensusCluster ###

Please use the following citation if you find ConsensusCluster useful in your work:

Michael Seiler, C. Chris Huang, Sandor Szalma, Gyan Bhanot. _OMICS: A Journal of Integrative Biology._ February 2010, 14(1): 109-113.