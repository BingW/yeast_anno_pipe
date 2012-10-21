#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"
# dependents: scer.gff,scer.fsa file in SGD_dataset/ 
# biopython
# bing

from Bio import SeqIO
from bing import gff_parse

flanking = 0

#indicate the scer.* files 
scergfffile=open("/Users/bingwang/zen/yeast_anno_pipe/SGD_data/scer.gff")
scerfsafile="/Users/bingwang/zen/yeast_anno_pipe/SGD_data/scer.fsa"

#initial datasets
scergff = []
for record in gff_parse.gffIterator(scergfffile):
    scergff.append(record)

scerfsa = {}
for record in SeqIO.parse(scerfsafile,"fasta"):
    scerfsa[record.id] = record

for record in scergff:
    if record.strand == "+":
        record.seq = scerfsa[record.seqid].seq[\
                int(record.start)-1-flanking:int(record.end)+flanking]
    elif record.strand == "-":
        record.seq =  scerfsa[record.seqid].seq[\
                int(record.start)-1-flanking:int(record.end)+flanking].reverse_complement()
    else:
        if record.type in ["chromosome","ARS"]: continue
        raise ValueError("lack strand info")

f = open("/Users/bingwang/zen/yeast_anno_pipe/output/intron."+str(flanking)+"f.fsa","w")
count = 0
name_sp = []
for record in scergff:
    if record.type == "intron" and record.attributes["Name"][0].startswith("Y"):
        while record.attributes["Name"][0] in name_sp:
            i = 1
            record.attributes["Name"][0] = record.attributes["Name"][0] + "." + str(i)
            print record.attributes["Name"][0]
            i += 1
        name_sp.append(record.attributes["Name"][0])
        f.write(">"+record.attributes["Name"][0]+"\t"+record.seqid+"\t"+\
                record.start+"\t"+record.end+"\t"+record.strand+"\n")
        f.write(str(record.seq)+"\n")
        count += 1
print count
