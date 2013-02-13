#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"

from Bio import SeqIO

alignment_path = "/Users/bingwang/zen/yeast_anno_pipe/SeubVSSuva/alignment4/"
fsa_file = alignment_path + "SeubFM1318_5000_130121.fasta"
tab_file = alignment_path + "/SeubFM1318_5000_130121_contigs.tab"
output_fsa = alignment_path + "SeubFM1318_5000_130121_linked.fasta"

seq_dict = {}
for record in SeqIO.parse(fsa_file,"fasta"):
    seq_dict[record.id] = record

sequence = ""
f = open(tab_file)
f.readline()
f.readline()
for line in f:
    if line.strip():
        elements = line.split("\t")
        contig_id = elements[1]
        strand = elements[3]
        start = int(elements[4])
        end = int(elements[5])
        if strand == "complement":
            sequence += str(seq_dict[contig_id].seq.reverse_complement())
        else:
            sequence += str(seq_dict[contig_id].seq)
        assert end == len(sequence)
    else:
        break

g = open(output_fsa,"w")
g.write(">SeubFM1318_5000_130121\n")
g.write("%s"%sequence)

for record in SeqIO.parse(output_fsa,"fasta"):
    pass
SeqIO.write(record,output_fsa,"fasta")
