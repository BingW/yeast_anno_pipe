#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"

Out_Path = "/Users/bingwang/zen/yeast_anno_pipe/output/YGAPvsDevin/"

from Bio import SeqIO
import os

input_contig_file = "/Users/bingwang/zen/yeast_anno_pipe/suva.mergedscaf.fsa"
output_contig_file = "/Users/bingwang/zen/yeast_anno_pipe/SuvaCBS7001_3000.fasta"
minimal_contig_len = 3000

seq_list = []
for record in SeqIO.parse(input_contig_file,"fasta"):
    if len(record.seq) > minimal_contig_len:
        seq_list.append(record)

print len(seq_list)
SeqIO.write(seq_list,open(output_contig_file,"w"), "fasta")




