#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"
# dependents: 

from Bio import SeqIO
import os

def construct_homo():
    pillarfile="/Users/bingwang/zen/yeast_anno_pipe/YGOB_dataset/Pillars.tab"
    aafile="/Users/bingwang/zen/yeast_anno_pipe/YGOB_dataset/AA.fsa"
    ncrnafile = "/Users/bingwang/zen/yeast_anno_pipe/YGOB_dataset/ncRNA.fsa"
    rrnafile = "/Users/bingwang/zen/yeast_anno_pipe/YGOB_dataset/rRNA.fsa"
    snrnafile = "/Users/bingwang/zen/yeast_anno_pipe/YGOB_dataset/snRNA.fsa"
    snoranfile = "/Users/bingwang/zen/yeast_anno_pipe/YGOB_dataset/snoRNA.fsa"
    trnafile = "/Users/bingwang/zen/yeast_anno_pipe/YGOB_dataset/tRNA.fsa"

    aafsa = {}
    for record in SeqIO.parse(aafile,"fasta"):
            aafsa[record.id] = record

    otherfsa = {}
    for fsafile in [ncrnafile,rrnafile,snrnafile,snoranfile,trnafile]:
        for record in SeqIO.parse(fsafile,"fasta"):
            otherfsa[record.id] = record
    print "initial finished"

    f = open(pillarfile)
    for i,line in enumerate(f):
        genes = line.replace("\n","").split("\t")
        for gene in genes:
            if gene in aafsa:
                g = open("/Users/bingwang/zen/yeast_anno_pipe/output/homo/aa_"+\
                        str(i)+".fsa","w")
                break
            elif gene in otherfsa:
                g = open("/Users/bingwang/zen/yeast_anno_pipe/output/homo/other_"+\
                        str(i)+".fsa","w")
                break
            else:
                continue
        for gene in genes:
            if gene in aafsa:
                g.write(">"+gene+"\n")
                g.write(str(aafsa[gene].seq)+"\n")
            elif gene in otherfsa:
                g.write(">"+gene+"\n")
                g.write(str(otherfsa[gene].seq)+"\n")
            else:
                continue

homo_path = "/Users/bingwang/zen/yeast_anno_pipe/output/homo/"
better_homo = "/Users/bingwang/zen/yeast_anno_pipe/output/goodhomo/"
count_good = 0
for fsa_file in os.listdir(homo_path):
    fsa_count = 0
    for line in open(homo_path+fsa_file):
        if line.startswith(">"):
            fsa_count += 1
    if fsa_count > 14:
        count_good += 1
        #os.system("cp "+homo_path+fsa_file+" "+better_homo+fsa_file)
print count_good
