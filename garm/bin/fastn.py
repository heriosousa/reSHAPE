#!/usr/bin/env python

# Python module for dealing with fasta and/or fastq sequences

# to do:
#  - change gc() to have a window width option
#  - add 454/solexa recognition to seq_id2.... functions

import sys
import re
import os
import string
import mh12_utils
import random


#_____________________________________________________________________________
# Fasta class: for doing things with fasta sequences
#_____________________________________________________________________________
# TO DO:
# find wherever I used length() method and replace with overloaded len
# same as above, with subseq
class Fasta:
    def __init__(self,id_in,seq_in):
        self.id = id_in
        self.seq = seq_in

    def __str__(self):
        return ">" + self.id + "\n" + self.seq

    def __getitem__(self, key):
        try:
            return Fasta(self.id, self.seq[key])
        except IndexError:
            sys.stderr.write('Index out of range in Fasta slice\n')
            raise

    def __len__(self):
        return len(self.seq)


    # returns length of the sequence
    def length(self):
        return len(self.seq)

    # returns gc content of the sequence in [0,100]
    def gc(self):
        c_count = self.seq.upper().count('C')
        g_count = self.seq.upper().count('G')
        no_n_length = len(self.seq) - self.seq.upper().count('N')

        if (no_n_length == 0):
            return 0
        else:
            return 100.0 * (c_count + g_count) / no_n_length

    # returns count of GC, as opposed to % GC
    def gc_count(self):
        c_count = self.seq.upper().count('C')
        g_count = self.seq.upper().count('G')
        no_n_length = len(self.seq) - self.seq.upper().count('N')

        if (no_n_length == 0):
            return 0
        else:
            return c_count + g_count


    # returns the name of the reads mate
    def mate(self):
        return seq_id2mate(self.id)

    # returns the first part of the read (i.e. part before it tells you F or R)
    def id_header(self):
        return seq_id2firstpart(self.id)

    # add some info to the id
    def add_id_info(self, type, sep = '_'):
        if type == "gc":
            self.id += sep + str(round(self.gc(), 1))
        elif type == "length":
            self.id += sep + str(self.length())
        else:
            sys.exit("Error in fastn.Fasta.id_add_info. " + type + " not recognised")


    # reverse the read
    def rev(self):
        self.seq = self.seq[::-1]


    # reverse complement the read
    def revcomp(self):
        self.seq = self.seq.translate(string.maketrans("ATCGatcg", "TAGCtagc"))[::-1]

    # same as __str__, but adds line breaks for multi-line printing
    def multi_line_str(self, line_length):
        s = '>' + self.id + '\n'
        offset = 0

        while offset < len(self.seq):
            s += self.seq[offset:(offset + line_length)] + '\n'
            offset += line_length

        return s[:-1]


    # returns a Fasta which is a subsequence of the original
    def subseq(self, start, end):
        return  Fasta(self.id, self.seq[(start - 1):end])


#_____________________________________________________________________________
# Fastq class: subclass of Fasta
#_____________________________________________________________________________
class Fastq(Fasta):
    def __init__(self,id_in,seq_in,qual_in):
        self.id = id_in
        self.seq = seq_in
        self.qual = qual_in

    def __str__(self):
        return "@" + self.id + "\n" + self.seq + "\n+\n" + self.qual

    def __getitem__(self, key):
        try:
            return Fastq(self.id, self.seq[key], self.qual[key])
        except IndexError:
            sys.stderr.write('Index out of range in Fasta slice\n')
            raise

    # reverse the read
    def rev(self):
        self.seq = self.seq[::-1]
        self.qual = self.qual[::-1]

    # reverse complement the read
    def revcomp(self):
        self.seq = self.seq.translate(string.maketrans("ATCGatcg", "TAGCtagc"))[::-1]
        self.qual = self.qual[::-1]

    # returns a Fastq which is a subsequence of the original
    def subseq(self, start, end):
        return  Fastq(self.id, self.seq[(start - 1):end], self.qual[(start - 1):end])



# takes an fai file and dictionary.  Adds lengths of sequences from
# fai file into dictionary
def fai2length_hash(filename, d):
    f = mh12_utils.open_file_read(filename)

    for line in f:
        tmp = line.split()
        d[tmp[0]] = int(tmp[1])

    f.close()

#_____________________________________________________________________________
# determines whether a file is fasta, fastq or other
#
# usage: fasta_or_fastq(filename)
# returns: string = fasta|fastq|unknown
def fasta_or_fastq(filename):
    f = mh12_utils.open_file_read(filename)
    x = f.readline()
    f.close()

    if x.startswith('>'):
        return 'fasta'
    elif x.startswith('@'):
        return 'fastq'
    else:
        return 'unknown'
#_____________________________________________________________________________



def fasta2singleLine(infile, outfile):
    overwrite = False

    if infile == outfile:
        overwrite = True
        outfile = outfile + '-tmp-' +  str(random.randint(0,1000000))
        if infile.endswith('.gz'):
            outfile = outfile + '.gz'

        if os.path.exists(outfile):
            sys.exit('Yowzer! Unlikely to happen, but ' + outfile + ' already exists.  Aborting')

    regex = re.compile('^\d')

    fasta = True

    # determine if it's fasta or a fasta.qual file
    f_in = mh12_utils.open_file_read(infile)

    f_in.readline()
    tmp  = f_in.readline()

    if regex.match(tmp):
        fasta = False

    f_in.close()

    f_in = mh12_utils.open_file_read(infile)
    f_out = mh12_utils.open_file_write(outfile)

    first = True

    for line in f_in:
        if line.startswith("\n"):
            continue
        elif line.startswith(">"):
            if not first:
                f_out.write("\n")
            else:
                first = False

            f_out.write(line)
        else:
            f_out.write(line.rstrip())

            if not fasta:
                f_out.write(' ')

    f_out.write("\n")

    f_in.close()
    f_out.close()

    if overwrite:
        os.rename(outfile, infile)


# takes a fasta fiel and writes a new one with line breaks every
# line_length characters of each sequence (except for the last
# line of each sequence, obviously)
def fasta2multiline(infile, outfile, line_length):
    if infile == outfile:
        sys.exit('infile = outfile in fastn.fasta2multiline. Aborting')

    f_in = mh12_utils.open_file_read(infile)
    f_out = mh12_utils.open_file_write(outfile)

    while 1:
        seq = get_next_seq_from_file(f_in, 'fasta')

        if not seq:
            break

        print >> f_out, seq.multi_line_str(line_length)

    f_in.close()
    f_out.close()



#_____________________________________________________________________________
# checks for unneccessary linebreaks in a fasta file
# returns True if ok, False if bad
def fasta_singleline_ok(filename):
    line_count = 0
    seq_count = 0

    f = mh12_utils.open_file_read(filename)

    for line in f:
        line_count += 1

        if line.startswith('>'):
            seq_count += 1

            if 2 * (seq_count - 1) != line_count - 1:
                return False

    f.close()

    return True
#_____________________________________________________________________________

#_____________________________________________________________________________
# reads a fasta/q file and fills a dictionary with id => Fasta object
#
# usage: fastn2dictionary(filename, dictionary name)
# returns: nothing, but fills hash

def fastn2dictionary(fname,d, first_only=False):
    filetype = fasta_or_fastq(fname)
    f_in = mh12_utils.open_file_read(fname)

    while 1:
        seq = get_next_seq_from_file(f_in, filetype)

        if not seq:
            break

        if first_only:
            seq.id = seq.id.split()[0]

        d[seq.id] = seq

    f_in.close()
#_____________________________________________________________________________

#_____________________________________________________________________________
# fills a dictionary with the length of each seqeunce in the file
def fastn2lengthdic(fname, d, min_length = 1, ignoreN = False, first_only=False):
    filetype = fasta_or_fastq(fname)
    if filetype == "unknown":
        sys.exit("Unknown file format of " + fname + " in method mh12_utils.fastn2lengthdic")

    f = mh12_utils.open_file_read(fname)

    while 1:
        seq = get_next_seq_from_file(f, filetype)

        if not seq:
            break

        if len(seq) < min_length:
            continue

        if ignoreN:
           seq.seq = seq.seq.replace('N', '')
           seq.seq = seq.seq.replace('n', '')

        if first_only:
            seq.id = seq.id.split()[0]

        d[seq.id] = len(seq)

    f.close()
#_____________________________________________________________________________


#_____________________________________________________________________________
# splits a fastn file into separate files, one file per sequence
# Output files are called outprefx + sequence name + '.fasta' (or ....fastq)
# Spaces in fasta header converted to underscores for file naming
# returns list of files which were written
def fastn_splitter(fname, outprefix):
    file_list = []
    filetype = fasta_or_fastq(fname)
    if filetype == 'unknown':
        sys.exit('Unknown file format of ' + fname + ' in method fastn.fastn_splitter')

    f_in = mh12_utils.open_file_read(fname)


    # loop through file, writing each sequence to new file
    while 1:
        seq = get_next_seq_from_file(f_in, filetype)

        if not seq:
            break

        outname = outprefix + seq.id.replace(' ', '_')[1:] + '.' + filetype
        f_out = mh12_utils.open_file_write(outname)
        print >> f_out, seq
        f_out.close()
        file_list.append(outname)


    f_in.close()

    return file_list
#
#_____________________________________________________________________________


# takes a subset of a fasta or fastq file given by the ids in a set
# If complement = True, take the reads which are not in the set.
# If start_only=False, use whole name.
#    start_only='', then take everything up to first whitespace
#    start_only=string (!=''), take everything up to 1st occurence of string.
# Writes subset of reads to outfile
def fastn2subset(infile, ids, outfile, complement = False, start_only=False):
    filetype = fasta_or_fastq(infile)

    if filetype == 'unknown':
        sys.exit('File ' + infile + ' not recognised as a fasta/q')

    f_in = mh12_utils.open_file_read(infile)
    f_out = mh12_utils.open_file_write(outfile)

    while 1:
        seq = get_next_seq_from_file(f_in, filetype)

        if not seq:
            break

        if start_only == '':
            seq.id = seq.id.split()[0]
        elif start_only != False:
            seq.id = seq.id.split(start_only)[0]

        if (not complement and seq.id in ids) or (complement and seq.id not in ids):
            print >> f_out, seq

    f_in.close()
    f_out.close()
#_____________________________________________________________________________


# takes two fasta/q files and shuffles them
def file_shuffler(infwd, inrev, outfile):
    filetype = fasta_or_fastq(infwd)
    if filetype == 'unknown':
        sys.exit('File type not recognised' + infwd)

    lineskip = {'fasta': 2, 'fastq': 4}[filetype]

    f_infwd = mh12_utils.open_file_read(infwd)
    f_inrev = mh12_utils.open_file_read(inrev)
    f_out = mh12_utils.open_file_write(outfile)

    line_fwd = f_infwd.readline()
    line_rev = f_inrev.readline()

    while line_fwd:
        for i in range(lineskip):
            f_out.write(line_fwd)
            line_fwd = f_infwd.readline()

        for i in range(lineskip):
            f_out.write(line_rev)
            line_rev = f_inrev.readline()

    f_infwd.close()
    f_inrev.close()
    f_out.close()
#_____________________________________________________________________________


# takes fasta or q file and return hash of summary stats
def get_stats(infile, min_length=1, ignore_N=False):
    d = {}

    # get sorted list of sequence lengths from biggest to smallest
    lengths = {}
    fastn2lengthdic(infile, lengths, min_length = min_length, ignoreN = ignore_N)

    lengths = sorted(lengths.values())
    lengths.reverse()

    d['total_length'] = sum(lengths)
    d['number'] = len(lengths)
    d['mean_length'] = 1.0 * d['total_length'] / d['number']
    d['longest'] = max(lengths)

    # work out the N50
    cumulative_length = 0

    for i in range(d['number']):
        if cumulative_length >= d['total_length'] / 2.0:
            if i == 1:
                d['n50'] = lengths[i-1]
            else:
                d['n50'] = lengths[i-1]

            d['n50_n'] = i
            break

        cumulative_length += lengths[i]


    if 'n50' not in d:
        d['n50'] = cumulative_length
        d['n50_n'] = d['number']

    return d
#_____________________________________________________________________________


#_____________________________________________________________________________
# reads next sequence and returns a fasta/fastq object.
# returns None iff there isn't another sequence in the file
def get_next_seq_from_file(filehandle, seqtype):
    id = filehandle.readline().rstrip()[1:]

    if id == '':
        return None

    # easy if it's fastq: don't need to worry about linebreaks in the sequence
    if seqtype == 'fastq':
        seq = filehandle.readline().rstrip()
        qual = filehandle.readline()
        qual = filehandle.readline().rstrip()

        if (len(seq) != len(qual)):
            sys.stderr.write('Mismatch in sequence/quality length! I got this from file:\n')
            sys.stderr.write('id:' + id + '\nseq:  ' + seq + '\nqual: ' + qual + '\n')
            sys.exit(1);

        return Fastq(id, seq, qual)

    # if it's fasta, need to keep reading until hit the next header,
    # and then return the position in the file back to right place
    seq = ''
    pos = filehandle.tell()

    while 1:
        line = filehandle.readline()

        if line.startswith('>') or not line:
            filehandle.seek(pos)
            break

        seq += line.rstrip()
        pos = filehandle.tell()

    return Fasta(id, seq)
#_____________________________________________________________________________
# returns name of next sequence in file, leaves the posiion in the
# unchanged.  Returns None if there isn't a next sequence in the file
def peek_next_seq(filehandle):
    pos = filehandle.tell()
    id = filehandle.readline().rstrip()[1:]
    filehandle.seek(pos)
    return id
#_____________________________________________________________________________
# moves position in file to start of next sequence.  Like
# get_next_seq_from file, but doesn't put anything into memory.
# returns True if skipped a sequence, False if didn't (i.e. was at end
# of file to begin with)
#def skip_next_seq(filehandle, seqtype):
#    line = filehandle.readline()  # line with id in it
#    if not line:
#       return False
#
#    if seqtype == 'fastq':
#        for i in [1,2,3]:
#            line = filehandle.readline()
#            if not line:
#                print >> sys.stderr, 'Error in fastn.skip_next_seq().  Hit EOF in the middle of skipping a fastq sequence'
#                sys.exit(1)
#    elif seqtype == 'fasta':
#        pos = filehandle.tell()
#
#        while 1:
#            line = filehandle.readline()
#
#            if line.startswith('>') or not line:
#                filehandle.seek(pos)
#                return True
#
#            pos = filehandle.tell()
#    else:
#        print >> sts.stderr, 'Error in fastn.skip_next_seq().  Filetype must be in {fasta, fastq}'
#        sys.exit(1)
#_____________________________________________________________________________
# removes duplicated reads from infile, writes uniq file to outfile
def fastn2uniq(infile, outfile):
    filetype = fasta_or_fastq(infile)
    reads = set()

    f_in = mh12_utils.open_file_read(infile)
    f_out = mh12_utils.open_file_write(outfile)

    while 1:
        seq = get_next_seq_from_file(f_in, filetype)

        if not seq:
            break

        if seq.id not in reads:
            reads.add(seq.id)
            print >> f_out, seq

    f_in.close()
    f_out.close()
#_____________________________________________________________________________




#_____________________________________________________________________________
# splits a fastn file into separate files, with each file having up to a
# given numer of bases.  Files are named prefix.counter.fasta, starting
# from 1.
# Spaces in fasta header converted to underscores for file naming
# returns number of files written
def fastn_split(fname, outprefix, no_of_bases):
    # writes array of sequences to a file
    def write_file(fname, a):
        f = mh12_utils.open_file_write(fname)

        for seq in a:
            print >> f, seq

        f.close()


    filetype = fasta_or_fastq(fname)
    if filetype == 'unknown':
        sys.exit('Unknown file format of ' + fname + ' in function fastn.fastn_split')

    f_in = mh12_utils.open_file_read(fname)

    file_counter = 1
    seqs = [get_next_seq_from_file(f_in, filetype)]
    base_counter = len(seqs[0])

    # loop through file, writing out files when enough data gathered
    while 1:
        next_seq = get_next_seq_from_file(f_in, filetype)

        if not next_seq:
            write_file(outprefix + str(file_counter) + '.' + filetype, seqs)
            break

        # either write out sequences, or add the next one to array
        if base_counter + len(next_seq) >= no_of_bases:
            write_file(outprefix + str(file_counter) + '.' + filetype, seqs)
            seqs = [next_seq]
            base_counter = next_seq.length()
            file_counter += 1
        else:
            seqs.append(next_seq)
            base_counter +=  next_seq.length()

    f_in.close()

    return file_counter
#_____________________________________________________________________________






# converts a phred score to fastq character score
#_____________________________________________________________________________
def phred_score2fastq(x):
    return chr(x + 33)
#_____________________________________________________________________________

# converts a fastq character quality score to a phred score (integer)
#_____________________________________________________________________________
def fastq_score2phred(x):
    return ord(x) - 33
#_____________________________________________________________________________


# converts a string of fastq scores (chars) to phred scores (ints)
#_____________________________________________________________________________
def fastq_score_string2phred(s):
    return ' '.join([str(fastq_score2phred(x)) for x in s])
#_____________________________________________________________________________


# converts a string of phred scores (ints) to fastq scores (chars)
#_____________________________________________________________________________
def phred_score_string2fastq(s):
    return ''.join([phred_score2fastq(int(x)) for x in s.split()])
#_____________________________________________________________________________


#_____________________________________________________________________________
# takes a sequence id and returns its mate name
def seq_id2mate(seq):
    regex_abi_p = re.compile('\.p1k')
    regex_abi_q = re.compile('\.q1k')

    if regex_abi_p.search(seq):
        return seq.replace(".p1k", ".q1k")
    elif regex_abi_q.search(seq):
        return seq.replace(".q1k", ".p1k")
    else:
        sys.exit("Error in fastn.seq_is2mate: " + seq + " not recognised")


#_____________________________________________________________________________
# takes a sequence id and returns the head part (i.e. everything except the
# end part which defines whether it's forward or reverse)
def seq_id2firstpart(seq):
    regex_abi = re.compile('\.[pq]1k')

    # if it's abi, the read header should be everything before the first dot
    if regex_abi.search(seq):
        a = seq.rsplit('.',1)
    else:
        sys.exit("Error in fastn.seq_id2firstpart: " + seq + " not recognised")

    return a[0]

