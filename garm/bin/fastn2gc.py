#!/software/bin/python

# Script to make gc plot and text file of gc content of a
# fastq file.

import sys
from optparse import OptionParser
import mh12_utils
import fastn
from math import floor

parser = OptionParser('[options] <fasta or fastq file> <text outfile> <plot outfile.{png, pdf}>')
parser.disable_interspersed_args()
parser.set_defaults(noclean = False)

parser.add_option('--window', type='int', help='Output GC content in windows of width WINDOW.  If not used, GC content of each sequence is calculated [%defualt]')
parser.add_option('--noclean', action='store_true',  help="Don't delete temporary files [%default]")

(options, args) = parser.parse_args()


# die if not enough mandatory arguments given
if len(parser.rargs) != 3:
    parser.print_help()
    sys.exit(1)


options.infile = parser.rargs[0]
options.txtout = parser.rargs[1]
options.plotout = parser.rargs[2]


filetype = fastn.fasta_or_fastq(options.infile)

if filetype == 'unknown':
    sys.exit('File ' + infile + ' not recognised as a fasta/q')


gc_hist = dict(zip(range(101), [0] * 101))

f_in = mh12_utils.open_file_read(options.infile)
f_out = mh12_utils.open_file_write(options.txtout)


while 1:
    seq = fastn.get_next_seq_from_file(f_in, filetype)

    if not seq:
        break

    if options.window:
        i = 0

        while i < len(seq):
            tmp = fastn.Fasta(seq.id, seq.seq[i:i + options.window])
            gc = tmp.gc()
            gc_hist[floor(gc)] += 1
            print >> f_out, seq.id, str(i+1), gc
            i += options.window
    else:
        gc = seq.gc()
        gc_hist[floor(gc)] += 1
        print >> f_out, seq.id, gc

f_in.close()
f_out.close()

c = not options.noclean

mh12_utils.hist2Rplot(gc_hist, options.plotout, x_label='%GC content', y_label='Frequency', clean=not options.noclean)


