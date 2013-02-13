#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"

from Bio import SeqIO

alignment_path = "/Users/bingwang/zen/yeast_anno_pipe/SeubVSSuva/alignment4/"
fsa_file_1 = alignment_path + "Suvarum_sequence.fsa"
fsa_file_2 = alignment_path + "SeubFM1318_5000_130121.fasta"
aln_file = alignment_path + "alignment4"
temp = alignment_path + "temp"
open(temp,"w").write(open(aln_file).read().replace("=",""))
for record in SeqIO.parse(temp,"fasta"):
    pos = record.id
    (start,end) = pos.split(":")[1].split("-")
    assert (int(end)-int(start)+1) == len(record.seq)-record.seq.count("-")
    print start,end,"|",len(record.seq)



