#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"
# dependents: 

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import random
import os

#data_path = "/Users/kmeihua1017/Documents/G875data/data/"
#out_path = "/Users/kmeihua1017/Documents/G875data/upseq_95/"
data_path = "/Users/bingwang/zen/yeast_anno_pipe/YGOB_dataset/"
out_path = "/Users/bingwang/zen/get_up_out/"

os.system("mkdir -p %s"%(out_path+"1000random"))
fsa_file = (out_path + "Klactis_62_homo_1000_up.fsa")
seq_dict = {}
for record in SeqIO.parse(fsa_file,"fasta"):
    seq_dict[record.id] = record.seq

for i in xrange(1000):
    records = []
    file_name = out_path+"1000random/Klactis_homo_random.%d.fsa"%(i+1)
    for name in random.sample(seq_dict.keys(),20):
        records.append(SeqRecord(seq_dict[name],name))
    SeqIO.write(records, file_name, "fasta")


input_path = "/Users/bingwang/zen/get_up_out/1000random/"
freq_file = "/Users/bingwang/zen/get_up_out/Klactis_5288_genes_1000_up.freq"
for file_name in os.listdir(input_path):
    if file_name.endswith(".fsa") and "random." in file_name:
        input_file = input_path+file_name
        out_file = input_path+file_name.replace("random","meme")
        cmd = "meme %s -bfile %s -dna -mod zoops -revcomp -minw 8 -maxw 13 -text | cat > %s" % (input_file,freq_file,out_file)
        os.system(cmd)



