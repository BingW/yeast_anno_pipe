#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import random
import os

species = "Skudriavzevii"
work_path = "/Users/bingwang/zen/get_up_out/"
n = 1
sp_homo_files = [f for f in os.listdir(work_path) if species in f and "homo_1000_up.fsa" in f]
assert len(sp_homo_files) == 1
sp_homo_file = work_path+sp_homo_files[0]
sp_freq_files = [f for f in os.listdir(work_path) if species in f and "1000_up_freq.txt" in f]
assert len(sp_freq_files) == 1
sp_freq_file = work_path+sp_freq_files[0]

out_path = "%s%s_1000random/"%(work_path,species)
cmd = "mkdir -p %s"%out_path
os.system(cmd)
seq_dict = {}
for record in SeqIO.parse(sp_homo_file,"fasta"): 
    seq_dict[record.id] = record.seq

for i in xrange(n,1000):
    records = []
    file_name = out_path + "%s_homo_random.%d.fsa"%(species,i+1)
    out_file = out_path + "%s_homo_meme.%d.out"%(species,i+1)
    for gene_name in random.sample(seq_dict.keys(),20):
        records.append(SeqRecord(seq_dict[gene_name],gene_name))
    SeqIO.write(records, file_name, "fasta")
    cmd = "meme %s -bfile %s -dna -nmotifs 3 -mod zoops -revcomp -minw 8 -maxw 13 -text | cat > %s" % (file_name,sp_freq_file,out_file)
    os.system(cmd)

