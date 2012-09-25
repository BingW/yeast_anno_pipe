#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"
# dependents: Scer.* file in SSS_dataset/ 

from Bio import SeqIO
from Bio.Emboss.Applications import NeedleCommandline
from bing import gff_parse
from bing import SGDfeature_parse

#indicate the scer.* files 
scergfffile=open("/Users/bingwang/zen/yeast_anno_pipe/SSS_dataset/Scer.gff")
SGDfeaturefile=open("/Users/bingwang/zen/yeast_anno_pipe/SSS_dataset/SGD_features.tab")
scerfsafile="/Users/bingwang/zen/yeast_anno_pipe/SSS_dataset/Scer.fsa"
scerscaffile="/Users/bingwang/zen/yeast_anno_pipe/SSS_dataset/Scer.mergedscaf"

#initial datasets
scergff = {}
for record in gff_parse.gffIterator(scergfffile):
    scergff[record.attributes['Gene'][0]] = record

SGDfeature = {}
for record in SGDfeature_parse.featureIterator(SGDfeaturefile):
    SGDfeature[record.SGDID] = record

scerfsa = {}
for record in SeqIO.parse(scerfsafile,"fasta"):
    scerfsa[record.id] = record

scerscaf = {}
for record in SeqIO.parse(scerscaffile,"fasta"):
    scerscaf[record.id] = record

#get seq and originseq
for id in scergff:
    #TODO scergff.phase
    scergff[id].seq = scerfsa[id].seq
    if scergff[id].strand == "+":
        scergff[id].originseq = scerscaf["Scer_"+scergff[id].seqid].seq[\
                int(scergff[id].start)-1:int(scergff[id].end)]
    elif scergff[id].strand == "-":
        scergff[id].originseq =  scerscaf["Scer_"+scergff[id].seqid].seq[\
                int(scergff[id].start)-1:int(scergff[id].end)].reverse_complement()
    else:
        raise ValueError("lack strand info")

def compareHead(origin,seq):
    for i,char in enumerate(seq):
        if char != origin[i]:
            if origin.endswith(seq[i+1:]):
                return i
            else:
                print "two introns or error"
                return 1
    return 0



for id in scergff:
    intron = compareHead(scergff[id].originseq[:20],scergff[id].seq[:20])
    if intron:
        print id,"introns:",scergff[id].attributes["introns"][0],"phase:",scergff[id].phase
        print intron
        print scergff[id].originseq,"\n\n"
        print scergff[id].seq


