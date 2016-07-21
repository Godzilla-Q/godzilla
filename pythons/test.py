#!/usr/bin/env python

"""Calculate accessibilities of RRS and start codon for both temperatures."""

# from __future__ import print_function
import sys
import RNA

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
#    fc_low.constraints_add(constr1, RNA.CONSTRAINT_DB_DEFAULT)
#    RRS_fc_low = fc_low.constraints_add(constr1, RNA.CONSTRAINT_DB_DEFAULT)
#    AUG_fc_low = fc_low.constraints_add(constr2, RNA.CONSTRAINT_DB_DEFAULT)    
#    RRS_fc_high = fc_high.constraints_add(constr1, RNA.CONSTRAINT_DB_DEFAULT)
#    AUG_fc_high = fc_high.constraints_add(constr2, RNA.CONSTRAINT_DB_DEFAULT)    
    # calculate mfe structure and mfe (no constraints)
    struct_low, mfe_low = fc_low.mfe()
    struct_high, mfe_high = fc_high.mfe()
    # calculate partition function for both temperatures and all constraints
    pfstruct_low, pf_low = fc_low.pf()
    pfstruct_high, pf_high = fc_high.pf()
    """    RRS_pfstruct_low, RRS_pf_low = RRS_fc_low.pf()
    RRS_pfstruct_high, RRS_pf_high = RRS_fc_high.pf()
    AUG_pfstruct_low, AUG_pf_low = AUG_fc_low.pf()
    AUG_pfstruct_high, AUG_pf_high = AUG_fc_high.pf()"""
    # print sequence name, mfes and partition functions for both temperatures
#    print RRS_fc_low
    print seq_file.id
    print mfe_low, mfe_high
    print pf_low, pf_high
    print constr1
    print constr2
    
print '\n_____________________________\nBiopython version %s' % (Bio.__version__)
print 'ViennaRNA Package version %s' % (RNA.__version__)



