#!/usr/bin/env python

import sys
from optparse import OptionParser
from fastn import fastn_split

parser = OptionParser('[options] <fasta or fastq> <out prefix> <#bases per split file>')
parser.disable_interspersed_args()
parser.set_defaults()

(options, args) = parser.parse_args()

if len(parser.rargs) != 3:
    parser.print_help()
    sys.exit(1)

options.infile = parser.rargs[0]
options.preout = parser.rargs[1]
options.bases = int(parser.rargs[2])

fastn_split(options.infile, options.preout, options.bases)
