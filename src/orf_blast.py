#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"

import os
from Bio import SeqIO

orf_fsa = "/Users/bingwang/zen/yeast_anno_pipe/output/orf.fsa"
orf_out = "/Users/bingwang/zen/yeast_anno_pipe/output/orf_VS_ygob.out"
#scer_db = "/Users/bingwang/zen/yeast_anno_pipe/ncbi_ygob_db/Scer"
ygob_db = "/Users/bingwang/zen/yeast_anno_pipe/ncbi_ygob_db/AA"
possible_orf = {}
for record in SeqIO.parse(orf_fsa,"fasta"):
    possible_orf[record.id] = record

def blastp(query,out,db,evalue=0.001,max_targets=1,threads=4,outfmt=None, testmod=False):
    if outfmt == None:
        outfmt="6 qseqid sseqid pident length mismatch gapopen "+\
                "qstart qend sstart send evalue bitscore"
    cmd="blastp -db %s -query %s -out %s " % (db,query,out) +\
        "-evalue %f -max_target_seqs %d " % (evalue,max_targets)+\
        "-num_threads %d -outfmt \"%s\"" % (threads,outfmt)
    if testmod:
        print cmd
    else:
        os.system(cmd)

blastp(orf_fsa,orf_out,ygob_db,testmod=False)
