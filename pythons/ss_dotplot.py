#!/usr/bin/env python

"""Calculate the mfes and the partition function of the sequences for both temperatures."""

# from __future__ import print_function
import sys
import RNA
import matplotlib
import Bio
import re # for regular expressions
import matplotlib.pyplot as plt
from math import exp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from matplotlib.backends.backend_pdf import PdfPages

md_low = RNA.md()
md_high = RNA.md()

md_low.temperature = 30
md_high.temperature = 37
RRS = 'AAGGAG'
spacer = 6
start = 'AUG'
k = 1.9872041e-3

def str2filename ( string ):
    # replace every non-alphanumeric character (except _ and -) with "-"
    # turn string into 'raw string'
    # strip all '-' from start and end
    rawstring = repr(string)
    newstring = re.sub(r'[\W]', '-',rawstring)
    newstring = newstring.strip('-')
    # cut at 31 characters, remove remaining '-' or '_' at the end
    newstring = newstring[:31]
    newstring = newstring.strip('-_')
    return newstring

def plot_2bppms ( bppm_low, bppm_high, name ):
    # call str2filename to make sure the filename is okay
    goodname = str2filename(name)
    # initialize the pdf file
    pp = PdfPages('{}_bppm.pdf'.format(goodname))
    # plot base pair probability matrix, include plottitle, write plot to pdf
    plt.matshow(bppm_low, fignum='low', cmap=plt.get_cmap('BuPu'))
    plt.title(md_low.temperature)
    pp.savefig()
    plt.matshow(bppm_high, fignum='high', cmap=plt.get_cmap('BuPu'))
    plt.title(md_high.temperature)
    pp.savefig()
    pp.close()
    plt.close()
    return

def seqconstraints ( sequence, str_1, str_2, spacer ):
    sequence = sequence.upper() 
    str_1 = str_1.upper()  
    str_2 = str_2.upper()  
    constr_1 = re.sub(
        '{:s}(?=.{{{:d}}}{:s})'.format(str_1, spacer, str_2),
        'x'*len(str_1),
        sequence
    )
    constr_1 = re.sub(
            '[ACGTU]',
            '.',
            constr_1
        )
    constr_2 = re.sub(
        '(?<={:s}.{{{:d}}}){:s}'.format(str_1, spacer, str_2),
        'x'*len(str_2),
        sequence
    )
    constr_2 = re.sub(
            '[ACGTU]',
            '.',
            constr_2
        )
    print constr_1, constr_2
    return constr_1, constr_2

def accessibility ( sequence, md, constr, pf_noconstr ):
    k = 1.9872041e-3
    fc_constr = RNA.fold_compound(sequence, md)
    fc_constr.constraints_add(constr, RNA.CONSTRAINT_DB_DEFAULT)
    pf_constr_struct, pf_constr = fc_constr.pf()
    if re.search('x', constr):
        acces = exp((pf_noconstr - pf_constr)/(k * (md.temperature + 273.15)))
    else:
        acces = 0
    return acces

def versions_used():
    return "\n_____________________________\n\
Biopython {:s}\n\
VRP {:s}\n\
matplotlib {:s}\n\
Python {:s}".format(Bio.__version__, RNA.__version__, matplotlib.__version__, sys.version)


print 'name mfe_low mfe_high pf_low pf_high RRS_acces_low RRS_acces_high AUG_acces_low AUG_acces_high'

for seq_file in SeqIO.parse(sys.stdin, "fasta"):
    # obtain the sequence string
    sequ = str(seq_file.seq)
    # generate the fold compounds
    fc_low = RNA.fold_compound(sequ, md_low)
    fc_high = RNA.fold_compound(sequ, md_high)
    # calculate mfe structure and mfe
    struct_low, mfe_low = fc_low.mfe()
    struct_high, mfe_high = fc_high.mfe()
    # calculate partition function
    pfstruct_low, pf_low = fc_low.pf()
    pfstruct_high, pf_high = fc_high.pf()
    # calculate base pair probability matrix
    bppm_low = fc_low.bpp()
    bppm_high = fc_high.bpp()
    # plot both bppms to one pdf
    plot_2bppms(bppm_low, bppm_high, seq_file.id)
    # plot mfe secondary structure
    RNA.PS_rna_plot(sequ, struct_low, '{:s}_low_ss.ps'.format(str2filename(seq_file.id)))    
    RNA.PS_rna_plot(sequ, struct_high, '{:s}_high_ss.ps'.format(str2filename(seq_file.id)))

    
    # calculate sequence constraints for RRS and AUG
    constr1, constr2 = seqconstraints(sequ,RRS,start,spacer)
    # calculate accessibilities with function
    RRS_acces_low = accessibility(sequ,md_low,constr1,pf_low)
    RRS_acces_high = accessibility(sequ,md_high,constr1,pf_high)
    AUG_acces_low = accessibility(sequ,md_low,constr2,pf_low)
    AUG_acces_high = accessibility(sequ,md_high,constr2,pf_high)
    # print sequence name, mfes and partition functions for both temperatures
    print seq_file.id, mfe_low, mfe_high, pf_low, pf_high, RRS_acces_low, RRS_acces_high, AUG_acces_low, AUG_acces_high
    print versions_used()

