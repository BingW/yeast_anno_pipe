#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"
# dependents: 

from Bio import SeqIO
scaffold = "/Users/bingwang/zen/yeast_anno_pipe/SGD_data/scer.fsa"

def find_orfs_with_trans(seq, trans_table=11, min_protein_length=100):
    answer = []
    seq_len = len(seq)
    for strand, nuc in [("+", seq), ("-", seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_begin = 0 #for irritate
            aa_start = 0 #for translate
            aa_end = 0
            while aa_begin < trans_len:
                aa_end = trans.find("*",aa_begin)
                aa_start = trans.find("M",aa_begin,aa_end)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_start > -1 and aa_end-aa_start >= min_protein_length:
                    if strand == "+":
                        start = frame+aa_start*3
                        end = min(seq_len,frame+aa_end*3+3)
                    else:
                        start = seq_len-frame-aa_end*3-3
                        end = seq_len-frame-aa_start*3                        
                    answer.append((start, end, strand,
                                   trans[aa_start:aa_end]))
                aa_begin = aa_end+1
    answer.sort()
    return answer

for record in SeqIO.parse(scaffold,"fasta"):
    orf_list = find_orfs_with_trans(record.seq)
    for start, end, strand, pro in orf_list:
        print "%s...%s - length %i, strand %s, %i:%i" \
              % (pro[:30], pro[-3:], len(pro), strand, start, end)

