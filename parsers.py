"""

Text parsers for creating SampleData objects


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

import sys
from numpy import array

class SampleData(object):
    """

    BaseCluster

        Usage:
            
            BaseParserObj.samples.append(SampleData(sample_id, sample_num, sample_class, data, index))
        
            sample_id       - Label for the sample.  Need not be unique. (Optional, but highly recommended)
            sample_num      - Another label (Optional)
            sample_class    - Yet another label, usually for known subclasses (Optional)
            data            - Data for this particular sample.  Will be converted into a numpy array by
                              BaseParser, regardless of what format data takes when assigned.  Required.
            index           - Yet another label (Optional)

        Properties

            cluster_id      - Id of cluster to which this sample belongs.  Generally assigned by clustering algorithms.


    """

    def __init__(self, sample_id=None, sample_num=None, sample_class=None, data=None, index=None):

        if data is None:
            self.data = []
        else:
            self.data = data

        self.cluster_id = None
        self.sample_class = sample_class
        self.sample_id = sample_id
        self.sample_num = sample_num
        self.index = index

    def __getitem__(self, x):

        return self.data[x]


class BaseParser(object):
    """

    BaseParser

        Common methods and actions for all parsers
        Designed to be subclassed.  All parsers must do the following:
        
        1: Have an __init__ method which takes data_file as an argument

        2: Have a _parse_data_file method which adds SampleData objects to self.samples, using
           the data from data_file.  Optionally, this can return a list of gene_ids which
           have incomplete information.  All data from this gene_id will be stripped from the
           final dataset.

        3: If the data vector entries have labels (e.g., gene probe names), they should be assigned
           to self.gene_names as a list, in the same length and order as the data vectors in each sample
    
        4: Be called "ParseSOMETHING" where SOMETHING is any name
           The gtk front end looks for ParseX names and puts them in a drop-down list for convenience

        Properties

            samples     - A list of SampleData objects
            gene_names  - A list of data vector entry labels

            raw_data    - Boolean var.  If the data needs additional preprocessing and you don't want
                          it expressed as a numpy array now, set this to True.  The default behaviour
                          (False) is PROBABLY what you want.  Unless you're reading sequence data, which
                          numpy handles poorly.  In that case, make sure you make them numpy arrays in
                          the _preprocess method of your main class following appropriate numerical
                          conversion.

    """

    def __init__(self, data_file, raw_data=False):

        self.samples = []  #List of SampleData instances
        self.gene_names = []
        
        data_handle = open(data_file, 'r')

        incomplete = self._parse_data_file(data_handle)

        data_handle.close()

        #Convert to NumPy arrays:
        if not raw_data:
            for sample in self.samples:
                sample.data = array(sample.data)

        self.gene_names = array(self.gene_names)

        if incomplete and len(self.gene_names):
            incomplete = dict.fromkeys(incomplete) #kill duplicates

            if len(incomplete) >= len(self.gene_names):
                raise ValueError, 'No complete genes or incorrectly reported incomplete data!'

            kept_indices = tuple([ x for x in xrange(len(self.gene_names)) if self.gene_names[x] not in incomplete ])

            self.gene_names = self.gene_names.take(kept_indices)
            
            for sample in self.samples:
                sample.data = sample.data.take(kept_indices)

    def _parse_data_file(self, data_handle):
        """Parse datafile into sample name<->number pairs and load data"""
        pass

    def __delitem__(self, x):

        del self.samples[x]

    def __getitem__(self, x):

        return self.samples[x]


class ParseRaw(BaseParser):
    """

    ParseRaw
        
        Parses GSE data I received.  I'm not certain if this is the default file format.

    """

    def __init__(self, data_file):

        self.column_headers = None

        BaseParser.__init__(self, data_file)

    def _parse_data_file(self, data_handle):
        """Parse datafile into sample name<->number pairs and load probe data"""

        for line in data_handle:
            a = line.split()

            if a:

                if a[0][0] == '#':
            
                    if len(a) == 4:  #Name<->Sample number map
                        sample_id = a[0][1:]
                        sample_num = a[-1]

                        self.samples.append(SampleData(sample_id=sample_id, sample_num=int(sample_num)))

                else:
                    if not self.column_headers:
                        self.column_headers = a[1:]

                    else:
                        self.gene_names.append(a[0])

                        for i in range(len(self.column_headers)):    #FIXME: If this isn't conveniently in the same order, it's broke
                            self.samples[i].data.append(float(a[i + 1]))


class ParseWang(BaseParser):
    """

    ParseWang

        Strange data parser

    """

    def __init__(self, data_file):

        BaseParser.__init__(self, data_file)

    def _parse_data_file(self, data_handle):
        """Parse datafile into sample name<->number pairs and load probe data"""

        for line in data_handle:
            a = line.split()

            if a:

                if a[0] == 'Original':
            
                    for label in a[5:]:
                        sample_class, sample_id = label.split('_')
                        
                        self.samples.append(SampleData(sample_id=sample_id, sample_class=sample_class))

                elif a[0][:5] == 'label':
                    pass

                else:
                    self.gene_names.append(a[1])

                    for i in range(len(self.samples)):
                        self.samples[i].data.append(float(a[i + 3]))


class ParseClean(BaseParser):
    """

    ParseClean

        Yet another parser for strange data

    """

    def __init__(self, data_file):

        BaseParser.__init__(self, data_file)

    def _parse_data_file(self, data_handle):
        """Parse datafile into sample name<->number pairs and load probe data"""

        sample_classes = []
        sample_ids = []

        for line in data_handle:
            a = line.split("\t")

            if not a[0] and a[6]:
                if not a[1][:4] == 'Gene':
                    for s_class in a[6:]:
                        if s_class:
                            if s_class[-1] == '\n':
                                s_class = s_class[:-1]

                            sample_classes.append(s_class)
                else:
                    for s_id in a[6:]:
                        if s_id:
                            if s_id[-1] == '\n':
                                s_id = s_id[:-1]

                            sample_ids.append(s_id)

                    if len(sample_ids) != len(sample_classes):
                        print "Error! Sample ids and sample classes have different lengths!"
                        print "Sample_ids: %s, Sample_classes: %s" % (len(sample_ids), len(sample_classes))
                        sys.exit(1)
                    
                    for i in range(len(sample_classes)):
                        self.samples.append(SampleData(sample_id=sample_ids[i], sample_class=sample_classes[i]))

            elif not a[0]:
                pass

            else:
                self.gene_names.append(a[1])

                values = []

                for value in a[6:]:
                    if value:
                        values.append(value)

                if len(values) != len(self.samples):
                    print "Error! Data values and stored samples have different lengths!"
                    print "Values: %s, self.samples: %s" % (len(values), len(self.samples))
                    print "Lines of data so far: %s" % len(self.samples[0].data)
                    print "Line in question: %s" % a
                    sys.exit(1)

                for i in range(len(self.samples)):
                    self.samples[i].data.append(float(values[i]))


class ParseNormal(BaseParser):
    """

    ParseNormal

        This is the one you use when the data is just a table with no special characteristics
        Probes in rows, samples in columns.  Sample ids are in row 0, and gene ids
        are in column 0.
        
    """

    def __init__(self, data_file):

        BaseParser.__init__(self, data_file)

    def _parse_data_file(self, data_handle):
        """Parse datafile into sample name<->number pairs and load probe data"""

        incomplete = []

        for line in data_handle:
            a = line.split("\t")

            if a:

                if not self.samples:
                    for sample in a[1:]:
                        sample = sample.strip()

                        self.samples.append(SampleData(sample_id=sample))

                else:
                    self.gene_names.append(a[0])
                    
                    for i in range(len(self.samples)):
                        try:
                            self.samples[i].data.append(float(a[i+1]))
                        except:
                            print "Incomplete data for sample", self.samples[i].sample_id, "gene", a[0]
                            self.samples[i].data.append(0.0)
                            incomplete.append(a[0])

        return incomplete


class ParsePartek(BaseParser):
    """

    ParsePartek

        This parser is used when there's a table which has the probes in columns and the samples on rows
        The first row contains probe names, and there aren't any sample names, such as in a partek output file

    """

    def __init__(self, data_file):

        BaseParser.__init__(self, data_file)

    def _parse_data_file(self, data_handle):

        count = 0

        for line in data_handle:
            a = line.split("\t")
            if a:

                if not self.gene_names:
                    for name in a:
                        name = name.strip()

                        self.gene_names.append(name)

                else:
                    self.samples.append(SampleData(sample_id = 'Sample ' + str(count), data = [ float(x) for x in a ]))
                    count += 1


class ParseSTI(BaseParser):
    """

    ParseSTI
        
        Parse the table from Supplemental Table STI formatted as a text file using tab delimits

        Usage:
            
            sdata = parsers.ParseSTI(data_file)

            Where data_file is the filename containing the STI table.

    """

    def __init__(self, data_file):
        
        self.refs = []

        BaseParser.__init__(self, data_file, True)

    def _parse_data_file(self, data_handle):
        """Parse datafile into sample name<->number pairs and load probe data"""

        header_flag = False

        for line in data_handle:
            a = line.split("\t")

            if a:
                if not header_flag:
                    header_flag = True

                else:
                    seq = a[4].strip()
                    
                    if not a[0] or a[0] == '0':
                        self.refs.append(SampleData(sample_id=a[1], sample_class=a[3], data=seq))
                    else:
                        self.samples.append(SampleData(sample_id=a[1], sample_class=a[3], data=seq))


class ParseKidney(BaseParser):

    def __init__(self, data_file):

        BaseParser.__init__(self, data_file)

    def _parse_data_file(self, data_handle):
        """Parse datafile into sample name<->number pairs and load probe data"""

        for line in data_handle:
            a = line.split("\t")

            if a:

                if not self.samples:
                    for sample in a[5:]:
                        sample = sample.strip()

                        self.samples.append(SampleData(sample_id=sample))

                elif a[4] == 'labels':
                    i = 0
                    for label in a[5:]:
                        self.samples[i].sample_class = a[5 + i]
                        i += 1

                else:
                    self.gene_names.append(a[1])
                    
                    for i in range(len(self.samples)):
                        try:
                            self.samples[i].data.append(float(a[i+5]))
                        except:
                            print "Incomplete data for sample", self.samples[i].sample_id, "gene", a[1], "...using 0.0"
                            self.samples[i].data.append(0.0)


def read_table(filename):
    """

    Read a simple conversion table between one name and another, useful for converting between probe names and gene symbols.

    This really should parse the entire CSV in the future, but memory concerns have held me back for now.
    Maybe an SQLite database?

    Returns a dict of first-column: second-column associations

    """

    conv = dict()

    handle = open(filename, 'r')

    for line in handle:
        a = line.split("\t")

        if a:

            conv[a[0]] = " - ".join(a[1:]).strip()

    return conv
