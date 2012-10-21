#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"

from Bio import SeqIO
from Bio import AlignIO
from Bio import Motif
import random
import os

def hmm_aln(fsafile,n=20,fixrange=(0,0),out_file="seq_aln.sth"):
    seq_dict = {}
    for record in SeqIO.parse(fsafile,"fasta"):
        seq_dict[record.id] = record
    max_intron = max([len(seq_dict[id].seq) for id in seq_dict])
    f = open(out_file,"w")
    for id in seq_dict:
        seq = list(seq_dict[id].seq)
        while len(seq) < max_intron:
            seq.insert(random.randint(fixrange[0],len(seq)-fixrange[1]),"-")
        f.write(">"+seq_dict[id].description+"\n")
        f.write("".join(seq)+"\n")
    f.close()
    alignment = AlignIO.read(out_file,"fasta")
    AlignIO.write([alignment],out_file,"stockholm")
    hmmfile = out_file + ".hmm"
    A_file = out_file
    B_file = out_file + ".B.sth"
    for i in range(n/2):
        os.system("hmmbuild "+hmmfile+" "+A_file)
        os.system("hmmalign -o "+B_file+" "+hmmfile+" "+A_file)
        os.system("hmmbuild "+hmmfile+" "+B_file)
        os.system("hmmalign -o "+A_file+" "+hmmfile+" "+B_file)
    alignment = AlignIO.read(out_file,"stockholm")
    AlignIO.write([alignment],fsafile+".aln", "fasta")
    os.system("rm "+out_file+" "+hmmfile+" "+B_file)

def dialign(fsafile,conf = None):
    if conf == None:
        conf = "/Users/bingwang/zen/bing/dialign/conf/"
    dialign = "dialign-tx "
    outfile = fsafile+".aln"
    cmd = dialign+" "+conf+" "+fsafile+" "+outfile
    os.system(cmd)

def separate_fsa(in_file,block,out_file_list=None):
    assert len(block[-1]) == 1
    if out_file_list == None:
        out_file_list = []
        for seq_range in block:
            try:
                out_file_list.append(in_file[:in_file.rfind(".")+1]+"%s_%s.fsa" % seq_range)
            except:
                out_file_list.append(in_file[:in_file.rfind(".")+1]+"%s.fsa" % seq_range)
    else:
        assert len(out_file_list) == len(block)
    seq = {}
    for record in SeqIO.parse(in_file,"fasta"):
        seq[record.id] = record
    for i,seq_range in enumerate(block[:-1]):
        f = open(out_file_list[i],"w")
        for id in seq:
            f.write(">"+seq[id].description+"\n")
            f.write(str(seq[id].seq[seq_range[0]-1:seq_range[1]])+"\n") #starts with 1
        f.close()
    f = open(out_file_list[-1],"w")
    for id in seq:
        f.write(">"+seq[id].description+"\n")
        f.write(str(seq[id].seq[block[-1][0]-1:])+"\n")
    f.close()
    return out_file_list

def merge_fsa(in_file_list,out_file):
    seq = {}
    for record in SeqIO.parse(in_file_list[0], "fasta"):
        seq[record.id] = record
    for in_file in in_file_list[1:]:
        for record in SeqIO.parse(in_file, "fasta"):
            seq[record.id].seq += record.seq
    f = open(out_file,"w")
    for id in seq:
        f.write(">"+seq[id].description+"\n")
        f.write(str(seq[id].seq)+"\n")
    f.close()
    return 1


intron_fsa = "/Users/bingwang/zen/yeast_anno_pipe/output/intron.50f.fsa"
fsa_file_list = separate_fsa(intron_fsa,[(1,50),(51,56),(57,-54),(-53,-50),(-49,)])

for in_file in [fsa_file_list[1],fsa_file_list[3]]:
    os.system("mv "+in_file+" "+in_file+".aln")
for in_file in [fsa_file_list[0],fsa_file_list[4]]:
    dialign(in_file)
for in_file in [fsa_file_list[2]:
    hmm_aln(in_file)

merge_fsa([fsa_file for fsa_file in fsa_file_list],intron_fsa+".aln")
for fsa_file in fsa_file_list:
    os.system("rm "+fsa_file+"*")
