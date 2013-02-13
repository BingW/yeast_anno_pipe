#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"

from Bio import SeqIO
import os

input_fsa = "/Users/bingwang/zen/yeast_anno_pipe/SGD_data/scer.fsa"
out_put_dir = "/Users/bingwang/zen/yeast_anno_pipe/SGD_data/single_fsa/"
os.system("mkdir %s"%(out_put_dir))
for record in SeqIO.parse(input_fsa,"fasta"):
    out_file = out_put_dir + record.id + ".fsa"
    SeqIO.write(record,open(out_file,"w"),"fasta")

