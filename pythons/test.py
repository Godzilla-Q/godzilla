#!/usr/bin/env python

"""Calculate accessibilities of RRS and start codon for both temperatures."""

# from __future__ import print_function
import sys
import RNA
from math import exp

import Bio
import re # for regular expressions

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

    
md_low = RNA.md()
md_high = RNA.md()

md_low.temperature = 30
md_high.temperature = 37
RRS = 'AAggAG'
spacer = 6
start = 'AUG'
k = 1.9872041e-3

def seqconstraints ( sequence, str_1, str_2, spacer ):
    sequence = sequence.upper() 
    str_1 = str_1.upper()  
    str_2 = str_2.upper()  
    constr_1 = re.sub(
        '%s(?=.{%s}%s)' % (str_1, spacer, str_2),
        'x'*len(str_1),
        sequence
    )
    constr_1 = re.sub(
            '[ACGTU]',
            '.',
            constr_1
        )
    constr_2 = re.sub(
        '(?<=%s.{%s})%s' % (str_1, spacer, str_2),
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
    k = 1.9872041e-3
    fc_constr = RNA.fold_compound(sequence, md)
    fc_constr.constraints_add(constr, RNA.CONSTRAINT_DB_DEFAULT)
    pf_constr_struct, pf_constr = fc_constr.pf()
    if re.search('x', constr):
        acces = exp((pf_noconstr - pf_constr)/(k * (md.temperature + 273.15)))
    else:
        acces = 0
    return acces

for seq_file in SeqIO.parse(sys.stdin, "fasta"):
    # obtain the sequence string
    sequ = str(seq_file.seq)
    # calculate sequence constraints for RRS and AUG
    constr1, constr2 = seqconstraints(sequ,RRS,start,spacer)
    # generate the fold compounds for both temperatures and each constraint
    fc_low = RNA.fold_compound(sequ, md_low)
    fc_high = RNA.fold_compound(sequ, md_high)
    # calculate mfe structure and mfe (no constraints)
    struct_low, mfe_low = fc_low.mfe()
    struct_high, mfe_high = fc_high.mfe()
    # calculate partition function for both temperatures and all constraints
    pfstruct_low, pf_low = fc_low.pf()
    pfstruct_high, pf_high = fc_high.pf()
    # calculate accessibilities with function
    RRS_acces_low = accessibility(sequ,md_low,constr1,pf_low)
    RRS_acces_high = accessibility(sequ,md_high,constr1,pf_high)
    AUG_acces_low = accessibility(sequ,md_low,constr2,pf_low)
    AUG_acces_high = accessibility(sequ,md_high,constr2,pf_high)
    if RRS_acces_low or RRS_acces_high or AUG_acces_low or AUG_acces_high:
        print '\n%s' % (seq_file.id)
    if RRS_acces_low:
        print RRS_acces_low
    if RRS_acces_high:
        print RRS_acces_high
    if AUG_acces_low:
        print AUG_acces_low
    if AUG_acces_high:
        print AUG_acces_high
    
print '\n_____________________________\nBiopython version %s' % (Bio.__version__)
print 'ViennaRNA Package version %s' % (RNA.__version__)



