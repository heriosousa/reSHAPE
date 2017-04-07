#!/software/bin/python

import sys
import subprocess
#from scipy import stats
import os.path
import re
import gzip

class Bsub:
    """Class for submitting bsub jobs on the farm"""
    def __init__(self, out, error, name, queue, mem, cmd, start=0, end=0, more="", depend=""):
        self.q = queue
        #self.m = str(mem)
        self.m = mem
        self.extra = more
        self.c = cmd

        # is this a job array?
        if 0 < start:
            self.o = out + '.%I'
            self.e = error + '.%I'
            self.n = '"' + name + '[' + str(start) + '-' + str(end) + ']"'

            # translate every appearance of 'INDEX' in command to \$LSB_JOBINDEX
            self.c = self.c.replace('INDEX', '\$LSB_JOBINDEX')
        else:
            self.o = out
            self.e = error
            self.n = name

        self.deps = ""

        if depend != "":
            self.add_dependency([depend])

    def __str__(self):
        # work out the memory in megs/gigs
        mb = int(self.m * 1000)
        gb = int(self.m * 1000000)

        #syscmd = "bsub  -E 'test -e /nfs/users/nfs_m/mh12/ && test -e /lustre/scratch103/sanger/mh12/' " + self.extra + " " + self.deps + " -q " + self.q + ' -R "select[type==X86_64 && mem > ' + str(mb) + '] rusage[mem=' + str(mb) + ']" -M' + str(gb)
        syscmd = "bsub -E 'test -e /nfs/users/nfs_m/mh12/'" + self.extra + " " + self.deps + " -q " + self.q + ' -R "select[type==X86_64 && mem > ' + str(mb) + '] rusage[mem=' + str(mb) + ']" -M' + str(gb)
        syscmd += ' -o ' + self.o + ' -e ' + self.e + ' -J ' + self.n + ' ' + self.c
        return syscmd

    def add_dependency(self, dep):
        self.deps = " -w'"

        for d in dep:
            self.deps += "done(\"" + d + "\") && "

        self.deps = self.deps[0:-3] + "' "


    def run(self):
        retcode = subprocess.call(str(self), shell=True)

        if retcode != 0:
            sys.exit("Error in bsub.run(). Command called:\n" + str(self))
#_____________________________________________________________________________

# takes a new-style artemis plot file and averages the values over a given
# window width, writing this new plot to a file
def artplot_new2windows(infile, outfile, window):
    # open files
    try:
        f_in = open(infile)
    except IOError:
        sys.exit("Error opening infile " + infile)

    try:
        f_out = open(outfile, 'w')
    except IOError:
        sys.exit("Error opening output file " + outfile)

    # print the 1st line, which is just a header
    line = f_in.readline()
    if line.startswith('#'):
        f_out.write(line)
        line = f_in.readline()
    else:
        print >> f_out, '#blah'

    pos = 1
    vals = []

    while line:
        # data[0] = position, data[1] = value of plot
        data = [int(x) for x in line.rstrip().split()]

        # if the position of the current line lies later than the
        # end of the current window being calculated
        while data[0] >= pos + window:
            f_out.write(str(pos) + "\t0\n")
            pos += window

        while data and line and data[0] < pos + window:
            vals.append(data[1])
            line = f_in.readline()
            if line:
                data = [int(x) for x in line.rstrip().split()]

        # work out mean plot value for the window, print it
        mean = sum(data) / window
        f_out.write(str(pos) + "\t" + str(round(1.0 * sum(vals) / window, 1)) + "\n")

        pos += window

        line = f_in.readline()
        if line:
            vals = [data[1]]

    f_in.close()
    f_out.close()
#_____________________________________________________________________________


# takes a cigar string and returns the number of bases which were hard clipped
# from the ends.  Works for ssaha sam output
#
# usage: cigarstring2HSsum(cigar_string)
# returns: int, number of bases
def cigarstring2HSsum(s_in):
    hard_re = re.compile('\d+[HS]')
    hits = hard_re.findall(s_in)

    if len(hits) == 0:
        return 0
    else:
        total = 0

        for x in hits:
            total += int(x[:-1])

        return total
#_____________________________________________________________________________


# takes a cigar string and returns the length of the hit
def cigarstring2hitlength(s_in):
    length = 0
    p = re.compile('[0-9]+[A-Z]')

    for hit in p.findall(s_in):
        if hit[-1] in ['M','D']:
            length += int(hit[:-1])

    return length
#_____________________________________________________________________________


# fills a dictionary with sequence id => length, using an fai file from samtools
#
# usage: fai2dic(filename, dictionary)
# returns: nothing
def fai2dic(filename, d):
    f = open_file_read(filename)

    for line in f:
        # first field=id, second field=length
        info = line.split()
        d[info[0]] = int(info[1])

    f.close()
#_____________________________________________________________________________


# splits a file, bsub job array-friendly named
#
# usage: file_splitter(filename, outprefix, #lines per split file)
# returns: number of split files made
def file_splitter(filename, outprefix, lines):
    # check if file exists
    if not os.path.exists(filename):
        sys.exit(filename + " not found in mh12_utils.file_splitter! aborting")

    # split the file
    cmd  = "split -l " + str(lines) + " -d -a 5 " + filename + " " + outprefix + "-tmp"
    retcode = subprocess.call(cmd, shell=True)

    if retcode != 0:
        sys.exit("Error in split command in mh12_utils.file_splitter. Command called:\n" + cmd)


    # put the split files into a list
    files = subprocess.Popen(["ls"], stdout=subprocess.PIPE).communicate()[0].split()

    # rename the split files
    i = 0

    for x in files:
        if x.startswith(outprefix + "-tmp"):
            i += 1
            cmd = "mv " + x + " " + " " + x[:-9] + "." + str(i)
            #print cmd

            retcode = subprocess.call(cmd, shell=True)
            if retcode != 0:
                sys.exit("Error in renaming files in mh12_utils.file_splitter. Command called:\n" + cmd)

    return i
#_____________________________________________________________________________


# takes a string (expected to be a filename) and replaces every . with an
# underscore, except the last dot.  LaTeX can only handle image filenames
# with one dot (before the extension). Bah.
def filename2tex_friendly(s_in):
    tmp = s_in.split(".")

    return "_".join(tmp[0:-1]) + "." + tmp[-1]
#_____________________________________________________________________________



# plots a list in R
# usage: pist2Rplot(list, R file, plot file, x axis label, y axis label)
#  - optional:  type = string [barplot]
#               clean = booli [True] (if true, deletes the Rout file)
# returns: nothing
# NOTES:  plot type noyet supported, current;y only blot bar chart
def list2Rplot(l_in, r_file, plot_out, xlabel, ylabel, type="barplot", clean=True):
    try:
        f_r = open(r_file, 'w')
    except IOError:
        sys.exit("Error opening R file " + r_file + " in mh12_utils.list2Rplot\n")

    # in R, a = vector of y coordinates
    s = "a=c("

    for x in l_in:
        s += str(x) + ", "

    s = s[:-2] + ")\n" # remove extra comma at the end and close brackets

    s += "x=c(0:" + str(len(l_in) - 1) + ")\n" # x = vector of x coordinates
    s += 'pdf("' + plot_out + "\")\n"
    s += "  barplot(a, names.arg=x, ylab=\"" + ylabel + "\", xlab=\"" + xlabel + "\")\n"
    s += "dev.off()\n"

    f_r.write(s)
    f_r.close()

    # run the R script
    cmd = "R CMD BATCH " + r_file
    retcode = subprocess.call(str(cmd), shell=True)

    if retcode != 0:
        sys.exit("Error running R script in mh12_utils.list2Rplot. Command called:\n" + cmd)

    if clean:
       os.remove(r_file + "out")
#_____________________________________________________________________________


# open file for reading and return file handle.
#If ilfename ends with '.gz', assumes file is gzipped.
def open_file_read(filename, binary=False):
    if filename == '-':
        return sys.stdin
    if filename.endswith('.gz'):
        try:
            f = gzip.open(filename, 'rb')
        except IOError:
            sys.exit('Error opening for reading gzipped file ' + filename)
    else:
        try:
            if binary:
                f = open(filename, 'rb')
            else:
                f = open(filename)
        except IOError:
            sys.exit('Error opening for reading file ' + filename)

    return f
#_____________________________________________________________________________


# open file for writing and return file handle.
# If ilfename ends with '.gz', makes gzipped file.
def open_file_write(filename, binary=False):
    if filename.endswith('.gz'):
        try:
            f = gzip.open(filename, 'wb')
        except IOError:
            sys.exit('Error opening for writing gzipped file ' + filename)
    else:
        try:
            if binary:
                f = open(filename, 'wb')
            else:
                f = open(filename, 'w')
        except IOError:
            sys.exit('Error opening fo writing file ' + filename)

    return f
#_____________________________________________________________________________


def dna_code_converter(base_in):
    b = base_in.upper()

    d = {'R':'A/G',
         'Y':'C/T',
         'K':'G/T',
         'M':'A/C',
         'S':'C/G',
         'W':'A/T',
         'B':'C/G/T',
         'D':'A/G/T',
         'H':'A/C/T',
         'V':'A/C/G',
         'N':'A/C/G/T',
         'X':'X',
         'A':'A',
         'C':'C',
         'G':'G',
         'T':'T',
         'U':'U'
        }

    if b in d:
        return d[b]
    else:
        sys.exit("Unrecognised base " + base_in + " in mh12_utils.dna_code_converter. Aborting\n")
#_____________________________________________________________________________

def plot_test(out_prefix):
    try:
        figure()
        plot([1,2,3,4,5])
        xlabel('Test')
        ylabel('Test')
        savefig(out_prefix + '-test.pdf')
        savefig(out_prefix + '-test.png')
        close()
    except:
        print >> sys.stderr, 'Error making test plot file.  X11 connection bad?'
        sys.exit(1)



# given a sam flag, returns a string of summary of the flag
def samflag_translate(f_in):
    samflag_out = subprocess.Popen(["~mh12/bin/samflag", f_in], stdout=subprocess.PIPE).communicate()[0].split("\n")

    output = f_in + "|"

    for x in samflag_out:
        if x.startswith("[1]"):
            output += x[11:] + "|"

    return output[:-1]
#_____________________________________________________________________________


# given a sam flag, returns array of breakdown of the flag
def samflag2set(flag):
    if not isinstance(flag, int):
        sys.exit('Error in mh12_utils.samflag2set(). Argument ' + flag
              + ' not an int.  Type:' + type(flag).__name__)

    flags = [0x0001, 0x0002, 0x0004, 0x0008, 0x0010, 0x0020, 0x0040,
             0x0080, 0x0100, 0x0200, 0x0400]

    return set([n & flag for n in flags if n & flag > 0])
#_____________________________________________________________________________


class Samflag:
    # creates samflag from either an int, or anything which can be
    # converted to a set
    def __init__(self, input):
        if isinstance(input, int):
            self.number = input
        else:
            try:
                input = set(input)
            except TypeError:
                sys.exit('Cannot create Samflag from ' + input)

            self.number = sum(input)

        if self.number < 0:
            sys.exit('Error!  Cannot create negative Samflag!')

        self._update_breakdown()

    def __add__(self, other):
        return Samflag(self.number + other)

    def __radd__(self, other):
        return Samflag(self.number + other)

    def __str__(self):
        return str(self.number)

    def __repr__(self):
        return str(self.number) + ' ' + str(self.breakdown)

    def __len__(self):
        return len(slf.breakdown)

    def _update_breakdown(self):
        flags = [0x0001, 0x0002, 0x0004, 0x0008, 0x0010, 0x0020, 0x0040,
                 0x0080, 0x0100, 0x0200, 0x0400]

        self.breakdown = set([n & self.number for n in flags if n & self.number > 0])

    def __contains__(self, flag):
        return flag in self.breakdown

    # returns '+' or '-' for the strand
    def strand_char(self):
        if 0x0010 in self.breakdown:
            return '-'
        else:
            return '+'
#_____________________________________________________________________________


# takes a line from sam file and returns the flag as a Samflag object
def samline2samflag(s):
    if s:
        return Samflag(int(s.split()[1]))
    else:
        sys.exit('Error in samline2samflag.  Argument passed: ' + s)


# given a filename and a string, writes that string to a file
def string2file(filename, string):
    # open file
    try:
        f = open(filename, 'w')
    except IOError:
        sys.exit("Error opening file " + filename + " in method mh12_utils.string2file")

    # write string to file and close file
    f.write(string)
    f.close()
#_____________________________________________________________________________

# don't use this! Use sam.Sam insead!
class Samline:
    def __init__(self, s):
        data = s.rstrip().split('\t')

        self.qname = data[0]
        self.flag = Samflag(int(data[1]))
        self.rname = data[2]
        self.pos = int(data[3])
        self.mapq = int(data[4])
        self.cigar = data[5]
        self.mrname = data[6]
        self.mpos = int(data[7])
        self.isize = int(data[8])
        self.seq = data[9]
        self.qual = data[10]
        self.opt = {}

        if len(data) > 11:
            extra = ' '.join(data[11:]).split()

            for x in extra:
                tmp = x.split(':')

                if tmp[1] == 'i':
                    self.opt[tmp[0]] = (tmp[1], int(tmp[2]))
                else:
                    self.opt[tmp[0]] = (tmp[1], tmp[2])

    def __str__(self):
        opt_fields = ''
        tmp = []

        for x in self.opt:
            tmp.append(x + ':' + self.opt[x][0] + ':' + str(self.opt[x][1]))

        opt_fields = '\t'.join(tmp)

        s = '\t'.join([self.qname, str(self.flag), self.rname, \
                          str(self.pos), str(self.mapq), self.cigar, self.mrname, \
                          str(self.mpos), str(self.isize), self.seq, self.qual])

        if opt_fields:
            return s + '\t' + opt_fields
        else:
            return s

    # returns end coordinate in reference of the hit
    def get_end_pos(self):
        return self.pos + cigar2ref_length(self.cigar) - 1

    # returns the length of a hit in the reference defined by the cigar string
    def get_ref_hit_length(self):
        return cigar2ref_length(self.cigar)


    def strand_char(self):
        return self.flag.strand_char()

#_____________________________________________________________________________
# returns the length of a hit in the reference defined by the cigar string
def cigar2ref_length(cigar):
    pos = 0

    numbers = [int(x) for x in re.split('[A-Z]', cigar)[:-1]]
    letters = re.split('[0-9]*', cigar)[1:]

    for i in range(len(numbers)):
        if letters[i] in ['M', 'D', 'N']:
            pos += numbers[i]

    return pos
#_____________________________________________________________________________

# takes a hash and uses R to make a barplot (keys = x axis, values
# = height of bar).
# minimum/maximum options will deal with outliers, depending on chop:
#  - chop=true: delete the outliers
#  - chop=False: set the value of each outlier to be = min or max
#
# Type of output is determined by the extensio of outfile (either pdf or png)
def hist2Rplot(h, outfile, minimum=None, maximum=None, chop=False, x_label=None, y_label=None, title=None, clean=True):
    # chop down the values to min/max?
    if minimum != None:
        if minimum not in h:
            h[minimum] = 0

        to_delete = set()

        for x in h:
            if x < minimum:
                h[minimum] += h[x]
                to_delete.add(x)

        for x in to_delete:
            del h[x]

        if chop:
            del h[minimum]

    if maximum != None:
        if maximum not in h:
            h[maximum] = 0

        to_delete = set()

        for x in h:
            if x > maximum:
                h[maximum] += h[x]
                to_delete.add(x)

        for x in to_delete:
            del h[x]

        if chop:
            del h[maximum]



    if maximum != None:
        for x in h:
            if h[x] > maximum: h[x] = maximum

    # make data for R to plot
    f_out = open_file_write(outfile + '.R')

    # make vector for x axis bars from keys of hash and y axis heights
    k = h.keys()
    k.sort()

    x_data = 'x=c(' + ', '.join([str(z) for z in k]) + ')'
    y_data = 'y=c(' + ', '.join([str(h[z]) for z in k]) + ')'

    print >> f_out, x_data + '\n' + y_data

    #print >> f_out, 'x=c(' + ', '.join([str(z) for z in range(min(h), max(h) + 1)]) + ')'
    #s = 'y=c('

    #for i in range(min(h), max(h) + 1):
    #    if i in h:
    #        s += str(h[i]) + ', '
    #    else:
    #        s += '0, '

    #s = s[:-2] + ')'
    #print >> f_out, s

    # determine out file type
    if outfile.endswith('pdf'):
        ext = 'pdf'
    elif outfile.endswith('png'):
        ext = 'png'
    else:
        sys.exit('Aborting.  Must have .png or .pdf filename in mh12_utils.hist2Rplot()')

    # start the R plot
    print >> f_out, ext + '("' + outfile + '")'

    # Sort out labels
    barplot = '    barplot(y, names.arg = x'

    if x_label:
        barplot += ', xlab="' + x_label + '"'
    if y_label:
        barplot += ', ylab="' + y_label + '"'
    if title:
        barplot += ', main="' + title + '"'

    barplot += ')'

    print >> f_out, barplot
    print >> f_out, 'dev.off()'
    f_out.close()


    syscmd = 'R CMD BATCH ' + outfile + '.R'
    retcode = subprocess.call(syscmd, shell=True)

    if retcode != 0:
        sys.exit("Error running R.  Command called:\n" + syscmd)

    if clean:
        os.remove(outfile + '.R')
        os.remove(outfile + '.Rout')


