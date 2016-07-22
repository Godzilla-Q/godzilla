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

for seq_file in SeqIO.parse(sys.stdin, "fasta"):
    # obtain the sequence string
    sequ = str(seq_file.seq)
    # calculate sequence constraints for RRS and AUG
    constr1, constr2 = seqconstraints(sequ,RRS,start,spacer)
    # generate the fold compounds for both temperatures and each constraint
    fc_low = RNA.fold_compound(sequ, md_low)
    fc_high = RNA.fold_compound(sequ, md_high)
    RRS_fc_low = RNA.fold_compound(sequ, md_low)
    RRS_fc_low.constraints_add(constr1, RNA.CONSTRAINT_DB_DEFAULT)
    AUG_fc_low = RNA.fold_compound(sequ, md_low)
    AUG_fc_low.constraints_add(constr2, RNA.CONSTRAINT_DB_DEFAULT)    
    RRS_fc_high = RNA.fold_compound(sequ, md_high)
    RRS_fc_high.constraints_add(constr1, RNA.CONSTRAINT_DB_DEFAULT)
    AUG_fc_high = RNA.fold_compound(sequ, md_high)
    AUG_fc_high.constraints_add(constr2, RNA.CONSTRAINT_DB_DEFAULT)    
    # calculate mfe structure and mfe (no constraints)
    struct_low, mfe_low = fc_low.mfe()
    struct_high, mfe_high = fc_high.mfe()
    # calculate partition function for both temperatures and all constraints
    pfstruct_low, pf_low = fc_low.pf()
    pfstruct_high, pf_high = fc_high.pf()
    RRS_pfstruct_low, RRS_pf_low = RRS_fc_low.pf()
    RRS_pfstruct_high, RRS_pf_high = RRS_fc_high.pf()
    AUG_pfstruct_low, AUG_pf_low = AUG_fc_low.pf()
    AUG_pfstruct_high, AUG_pf_high = AUG_fc_high.pf()
    # calculate accessibilities
    RRS_acces_low = exp((pf_low - RRS_pf_low)/(k * (md_low.temperature + 273.15)))
    RRS_acces_high = exp((pf_high - RRS_pf_high)/(k * (md_high.temperature + 273.15)))
    AUG_acces_low = exp((pf_low - AUG_pf_low)/(k * (md_low.temperature + 273.15)))
    AUG_acces_high = exp((pf_high - AUG_pf_high)/(k * (md_high.temperature + 273.15)))

    # print sequence name, mfes and partition functions for both temperatures
    print seq_file.id
    print pf_low
    if pf_low != RRS_pf_low and pf_low != AUG_pf_low:
        print RRS_pf_low
        print AUG_pf_low
    print pf_high
    if pf_high != RRS_pf_high and pf_high != AUG_pf_high:
        print RRS_pf_high
        print AUG_pf_high
        print repr(RRS_acces_low)
        print repr(RRS_acces_high)
        print repr(AUG_acces_low)
        print repr(AUG_acces_high)
    print pfstruct_low
    print pfstruct_high
    print mfe_low, mfe_high
    print constr1
    print constr2
    print k
    print (md_low.temperature + 273.15)

    
print '\n_____________________________\nBiopython version %s' % (Bio.__version__)
print 'ViennaRNA Package version %s' % (RNA.__version__)



