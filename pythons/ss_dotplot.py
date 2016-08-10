#!/usr/bin/env python

"""Calculate the mfes and the partition function of the sequences for both temperatures."""

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

MODEL_LOW_TEMPERATURE = RNA.md()
MODEL_HIGH_TEMPERATURE = RNA.md()

MODEL_LOW_TEMPERATURE.temperature = 30
MODEL_HIGH_TEMPERATURE.temperature = 37
RRS = 'AAGGAG'
SPACER = 6
START = 'AUG'
BOLTZMANN_K = 1.9872041e-3

def str2filename ( string ):
    """ turn string into suitable filename """
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
    plt.title(MODEL_LOW_TEMPERATURE.temperature)
    pp.savefig()
    plt.matshow(bppm_high, fignum='high', cmap=plt.get_cmap('BuPu'))
    plt.title(MODEL_HIGH_TEMPERATURE.temperature)
    pp.savefig()
    pp.close()
    plt.close()
    return

def seqconstraints ( sequence, str_1, str_2, SPACER ):
    sequence = sequence.upper() 
    str_1 = str_1.upper()  
    str_2 = str_2.upper()  
    constr_1 = re.sub(
        '{:s}(?=.{{{:d}}}{:s})'.format(str_1, SPACER, str_2),
        'x'*len(str_1),
        sequence
    )
    constr_1 = re.sub(
            '[ACGTU]',
            '.',
            constr_1
        )
    constr_2 = re.sub(
        '(?<={:s}.{{{:d}}}){:s}'.format(str_1, SPACER, str_2),
        'x'*len(str_2),
        sequence
    )
    constr_2 = re.sub(
            '[ACGTU]',
            '.',
            constr_2
        )
    return constr_1, constr_2

def accessibility ( sequence, md, constr, pf_noconstr ):
    fc_constr = RNA.fold_compound(sequence, md)
    fc_constr.constraints_add(constr, RNA.CONSTRAINT_DB_DEFAULT)
    pf_constr_struct, pf_constr = fc_constr.pf()
    if re.search('x', constr):
        acces = exp((pf_noconstr - pf_constr)/(BOLTZMANN_K * (md.temperature + 273.15)))
    else:
        acces = 0
    return acces

def versions_used():
    return "\n_____________________________\n\
Biopython {:s}\n\
VRP {:s}\n\
matplotlib {:s}\n\
Python {:s}".format(Bio.__version__, RNA.__version__, matplotlib.__version__, sys.version)




def main():
    """ for sequence string, calculate mfe structure, mfe, pf, base pair probability matrix, plot structures and bppms, calculate accessibilities"""
    print 'name mfe_low mfe_high pf_low pf_high RRS_acces_low RRS_acces_high AUG_acces_low AUG_acces_high'
    for seq_file in SeqIO.parse(sys.stdin, 'fasta'):
         sequ = str(seq_file.seq)
         fc_low = RNA.fold_compound(sequ, MODEL_LOW_TEMPERATURE)
         fc_high = RNA.fold_compound(sequ, MODEL_HIGH_TEMPERATURE)
         struct_low, mfe_low = fc_low.mfe()
         struct_high, mfe_high = fc_high.mfe()
         pfstruct_low, pf_low = fc_low.pf()
         pfstruct_high, pf_high = fc_high.pf()
         bppm_low = fc_low.bpp()
         bppm_high = fc_high.bpp()
         
         plot_2bppms(bppm_low, bppm_high, seq_file.id)
         RNA.PS_rna_plot(sequ, struct_low, '{:s}_low_ss.ps'.format(str2filename(seq_file.id)))    
         RNA.PS_rna_plot(sequ, struct_high, '{:s}_high_ss.ps'.format(str2filename(seq_file.id)))
    
         constr1, constr2 = seqconstraints(sequ,RRS,START,SPACER)
         RRS_acces_low = accessibility(sequ,MODEL_LOW_TEMPERATURE,constr1,pf_low)
         RRS_acces_high = accessibility(sequ,MODEL_HIGH_TEMPERATURE,constr1,pf_high)
         AUG_acces_low = accessibility(sequ,MODEL_LOW_TEMPERATURE,constr2,pf_low)
         AUG_acces_high = accessibility(sequ,MODEL_HIGH_TEMPERATURE,constr2,pf_high)

         print seq_file.id, mfe_low, mfe_high, pf_low, pf_high, RRS_acces_low, RRS_acces_high, AUG_acces_low, AUG_acces_high
    print versions_used()

if __name__ == '__main__':
    sys.exit(main())
