#!/usr/bin/env python

"""Calculate accessibilities of RRS and start codon for both temperatures."""

import sys
import RNA
from math import exp

import Bio
import re # for regular expressions

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from ss_dotplot import versions_used
    
MODEL_LOW_TEMPERATURE = RNA.md()
MODEL_HIGH_TEMPERATURE = RNA.md()

MODEL_LOW_TEMPERATURE.temperature = 30
MODEL_HIGH_TEMPERATURE.temperature = 37
RRS = 'AAggAG'
SPACER = 6
START = 'AUG'
BOLTZMANN_K = 1.9872041e-3

def seqconstraints ( sequence, str_1, str_2, SPACER ):
    sequence = sequence.upper() 
    str_1 = str_1.upper()  
    str_2 = str_2.upper()  
    constr_1 = re.sub(
        '%s(?=.{%s}%s)' % (str_1, SPACER, str_2),
        'x'*len(str_1),
        sequence
    )
    constr_1 = re.sub(
            '[ACGTU]',
            '.',
            constr_1
        )
    constr_2 = re.sub(
        '(?<=%s.{%s})%s' % (str_1, SPACER, str_2),
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

def main():
    for seq_file in SeqIO.parse(sys.stdin, "fasta"):
        sequ = str(seq_file.seq)
        constr1, constr2 = seqconstraints(sequ,RRS,START,SPACER)
        fc_low = RNA.fold_compound(sequ, MODEL_LOW_TEMPERATURE)
        fc_high = RNA.fold_compound(sequ, MODEL_HIGH_TEMPERATURE)
        pfstruct_low, pf_low = fc_low.pf()
        pfstruct_high, pf_high = fc_high.pf()
        RRS_acces_low = accessibility(sequ,MODEL_LOW_TEMPERATURE,constr1,pf_low)
        RRS_acces_high = accessibility(sequ,MODEL_HIGH_TEMPERATURE,constr1,pf_high)
        AUG_acces_low = accessibility(sequ,MODEL_LOW_TEMPERATURE,constr2,pf_low)
        AUG_acces_high = accessibility(sequ,MODEL_HIGH_TEMPERATURE,constr2,pf_high)
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
    print versions_used()

if __name__ == '__main__':
    sys.exit(main())


