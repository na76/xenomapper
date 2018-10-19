#!/usr/bin/env python3
# encoding: utf-8
"""
xenomapper.py

A script for parsing pairs of sam files and returning sam files containing only reads where no better mapping exist in other files.
Used for filtering reads where multiple species may contribute (eg human tissue xenografted into mouse).

Created by Matthew Wakefield.
Copyright (c) 2011-2016  Matthew Wakefield, The Walter and Eliza Hall Institute and The University of Melbourne. All rights reserved.


   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


"""
import sys
import os
import argparse, textwrap
import subprocess
import re
from collections import Counter
from copy import copy

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2011-2016 Matthew Wakefield, The Walter and Eliza Hall Institute and The University of Melbourne"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Production/Stable"


def get_bam_header(bamfile): #pragma: no cover #not tested due to need for samtools
    p = subprocess.Popen('samtools view -H -',stdin=bamfile,stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
    header = []
    for line in p.stdout:
        header.append(line.decode('ascii'))
    bamfile.seek(0) #reset to start of file for next samtools call
    return [x.strip('\n') for x in header]


def bam_lines(f): #pragma: no cover #not tested due to need for samtools
    """Use samtools in a subprocess to yield lines of sam from a bam file
        Arguments: a file or file like object in bam format
        Yields:    ascii sam file lines
    """
    p = subprocess.Popen('samtools view -',stdin=f,stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
    header_lines = []
    for line in p.stdout:
        yield line.decode('ascii')


def getBamReadPairs(bamfile1, skip_repeated_reads=False): #pragma: no cover #not tested due to need for samtools
    """Process two bamfiles to yield the equivalent line from each file
        Arguments:
        bamfile1, bamfile2  - file or file like objects in binary bam format
                              containing the same reads in the same order
                              mapped in two different species
        Yields:    a tuple of lists of sam fields split on white space
    """
    bam1 = bam_lines(bamfile1)
    line1= next(bam1).strip('\n').split() #split on white space. Results in 11 fields of mandatory SAM + variable number of additional tags.
    while line1 and line1 !=['']:
        yield line1
        previous_read1 = line1[0]
        line1= next(bam1).strip('\n').split()
    pass


def add_pg_tag(sam_header_list,comment=None):
    new_header = copy(sam_header_list)
    if not [x[0] for x in new_header] == ['@',]*len(new_header):
        raise ValueError('Incorrect SAM header format :\n{0}'.format('\n'.join(new_header))) #pragma: no cover
    if new_header[-1][:3] == '@PG':
        PP = 'PP'+[x for x in new_header[-1].split() if x[:2] == 'ID'][0][2:]+'\t'
    else:
        PP = ''
    new_header.append('@PG\tID:Xenomapper\tPN:Xenomapper\t'+PP+'VN:{0}'.format(__version__))
    if comment:
        new_header.append('@CO\t'+comment)
    return new_header


def process_headers(file1, primary_specific=sys.stdout, primary_multi=None, unassigned=None):
    """Process headers from two sam or bam files and write appropriate
    header information to the correct output files
        Arguments:
        file1, file2  - file or file like objects in binary bam format
                        or ascii sam format
        primary_specific, secondary_specific, primary_multi,
        secondary_multi, unassigned, unresolved
                      - ascii file or file like objects for outputs
        bam           - Boolean flag indicating file1 & file2 are
                        in binary bam format.  Default = False
    """
    samheader1 = get_bam_header(file1)
    print('\n'.join(add_pg_tag(samheader1,
                    comment='species specific reads'
                    )), file=primary_specific)
    if primary_multi:
        print('\n'.join(add_pg_tag(samheader1,
                    comment='species specific multimapping reads'
                    )), file=primary_multi)
    if unassigned:
        print('\n'.join(add_pg_tag(samheader1,
                    comment='reads that could not be assigned'
                    )), file=unassigned)
    pass


def get_tag(sam_line,tag='AS'):
    """Return the value of a SAM tag field
    Arguments:
        sam_line  - list of elements from a SAM file line
        tag       - the name of the optional tag to be returned
                    Only suitable for numeric tags eg AS or XS
    Returns
        tag_value - the value of the SAM tag converted to a float
                    or -inf if tag is not present.
    """
    tag_list = [x for x in sam_line[11:] if tag in x]
    if not tag_list:
        return float('-inf') #this will always be worse than any bowtie score
    if len(tag_list) > 1:
        raise ValueError('SAM line has multiple values of {0}: {1}'.format(tag,sam_line)) #pragma: no cover
    return float(tag_list[0].split(':')[-1])


def get_mapping_state(AS1,XS1, min_score=float('-inf')):
    """Determine the mapping state based on scores in each species.
    Scores can be negative but better matches must have higher scores
    Arguments:
        AS1  - float or interger score of best match in primary species
        XS1  - float or interger score of other match in primary species
        AS2  - float or interger score of best match in secondary species
        XS2  - float or interger score of other match in secondary species
        min_score - the score that matches must exceed in order to be
                    considered valid matches. Note scores equalling this
                    value will also be considered not to match.
                    Default = -inf
    Returns
        state - a string of 'primary_specific', 'secondary_specific',
                'primary_multi', 'secondary_multi',
                'unresolved', or 'unassigned' indicating match state.
    """
    if AS1 <= min_score:  #low quality mapping
        return 'unassigned'
    elif AS1 > min_score: #maps in primary better than secondary
        if not XS1 or AS1 > XS1:       #maps uniquely in primary
            return 'primary_specific'
        else:            #multimaps in primary better than secondary
            return 'primary_multi'
    else: raise RuntimeError('Error in processing logic with values {0} '.format((AS1,XS1,AS2,XS2))) # pragma: no cover


def main_single_end(readpairs,
                    primary_specific=sys.stdout,
                    primary_multi=None,
                    unassigned=None,
                    min_score=float('-inf'),
                    tag_func=get_tag):
    """Main loop for processing single end read files
    Arguments:
        readpairs - an iterable of tuples of lists of sam fields
        primary_specific, secondary_specific, primary_multi,
        secondary_multi, unassigned, unresolved
                  - ascii file or file like objects for outputs
        min_score - the score that matches must exceed in order to be
                    considered valid matches. Note scores equalling this
                    value will also be considered not to match.
                    Default = -inf
        tag_func  - a function that takes a list of sam fields and a
                    tag identifier (at least 'AS' and 'XS')
                    returns a numeric value for that tag
    Returns:
        category_counts - a dictionary keyed by category containing
                    occurance counts
    """
    #assume that reads occur only once and are in the same order in both files

    category_counts = Counter()

    for line1 in readpairs:
        AS1 = tag_func(line1, tag='AS')
        XS1 = tag_func(line1, tag='XS')

        state = get_mapping_state(AS1,XS1,min_score)

        category_counts[state] += 1

        if state == 'primary_specific':
            if primary_specific:
                print('\t'.join(line1),file=primary_specific)
        elif state == 'primary_multi':
            if primary_multi:
                print('\t'.join(line1),file=primary_multi)
        elif state == 'unassigned':
            if unassigned:
                print('\t'.join(line1),file=unassigned)
        else: raise RuntimeError('Unexpected state {0} '.format(state)) # pragma: no cover
    return category_counts


def output_summary(category_counts, outfile=sys.stderr):
    print('-'*80, file=outfile)
    print('Read Count Category Summary\n', file=outfile)
    print('|       {0:45s}|     {1:10s}  |'.format('Category','Count'), file=outfile)
    print('|:','-'*50,':|:','-'*15,':|',sep='', file=outfile)
    for category in sorted(category_counts):
        print('|  {0:50s}|{1:15d}  |'.format(str(category),category_counts[category]), file=outfile)
    print(file=outfile)
    pass


def command_line_interface(*args,**kw): #pragma: no cover
    parser = argparse.ArgumentParser(prog = "xenomapper",
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    description=textwrap.dedent("""\
                    A script for parsing pairs of sam files and returning sam files
                    containing only reads where no better mapping exist in other files.
                    Used for filtering reads where multiple species may contribute
                    (eg human tissue xenografted into mouse, pathogen growing on plant).

                    Files should contain an AS and XS score and better matches must have
                    a higher alignment score (but can be negative).
                    Reads must be in the same order in both species.

                    In practice this is best acchieved by using Bowtie2 in --local mode.
                    If the -p option is used you must also use --reorder.

                    Limited support is provided for aligners that do not produce AS and XS
                    score tags via the --cigar_score option.

                    All input files must be seekable
                    (ie not a FIFO, process substitution or pipe)'
                    """),
                    epilog = textwrap.dedent("""\
                    To output bam files in a bash shell use process subtitution:
                        xenomapper --primary_specific >(samtools view -bS - > outfilename.bam)

                    This program is distributed in the hope that it will be useful,
                    but WITHOUT ANY WARRANTY; without even the implied warranty of
                    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n

                    """),
                    )
    parser.add_argument('--primary_bam',
                        type=argparse.FileType('rb'),
                        default=None,
                        help='a BAM format Bowtie2 mapping output file corresponding to the primary species of interest')
    parser.add_argument('--primary_specific',
                        type=argparse.FileType('wt'),
                        default=sys.stdout,
                        help='name for SAM format output file for reads mapping to a specific location in the primary species')
    parser.add_argument('--primary_multi',
                        type=argparse.FileType('wt'),
                        default=None,
                        help='name for SAM format output file for reads multi mapping in the primary species')
    parser.add_argument('--unassigned',
                        type=argparse.FileType('wt'),
                        default=None,
                        help='name for SAM format output file for unassigned (non-mapping) reads')
    args = parser.parse_args()
    return args


def main(): #pragma: no cover
    args = command_line_interface()

    process_headers(args.primary_bam,
                            primary_specific=args.primary_specific,
                            primary_multi=args.primary_multi,
                            unassigned=args.unassigned)

    readpairs = getBamReadPairs(args.primary_bam, skip_repeated_reads=False)
#
    category_counts = main_single_end(readpairs,
                        primary_specific=args.primary_specific,
                        primary_multi=args.primary_multi,
                        unassigned=args.unassigned,
                        tag_func=get_tag)

    output_summary(category_counts=category_counts)
    pass

if __name__ == '__main__':
    main()
