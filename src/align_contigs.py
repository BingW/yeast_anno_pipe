#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"

from Bio import SeqIO
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
import os

class Record():
    pass
def make_db(fsa,db):
    os.system("makeblastdb -in %s -dbtype nucl -out %s"%(fsa,db))

def run_blast(db,fsa,out_file):
    os.system("blastn -evalue 0.01 -max_target_seqs 5 "+\
        "-db %s -query %s -out %s "%(db,fsa,out_file) +\
        "-outfmt \"7 qseqid sseqid pident length mismatch gapopen qstart " +\
        "qend sstart send\"")

def blast_out_parse(handle):
    while 1:
        line = handle.readline()
        if not line:
            break
        r = Record()
        (r.qseqid,r.sseqid,r.pident,r.length,r.mismatch,r.gapopen,r.qstart,\
        r.qend,r.sstart,r.send,r.qseq,r.sseq) = line.split("\t")
        yield r

Ref_fsa = "/Users/bingwang/zen/yeast_anno_pipe/YGOB_dataset/Suvarum_sequence.fsa"
contig_fsa = "/Users/bingwang/zen/yeast_anno_pipe/SeubFM1318_5000_130121.fasta"
blast_output = "/Users/bingwang/zen/yeast_anno_pipe/SeubFM1318_5000_blast.out"
Ref_db = "/Users/bingwang/zen/yeast_anno_pipe/YGOB_dataset/Suvarum_db"

make_db(Ref_fsa,Ref_db)
run_blast(Ref_db,contig_fsa,blast_output)

ref_seq_dict = {}
for record in SeqIO.parse(Ref_fsa,"fasta"):
    ref_seq_dict[record.id] = record

contig_seq_dict = {}
for record in SeqIO.parse(contig_fsa,"fasta"):
    contig_seq_dict[record.id] = record

blast_result = {}
handle = open(blast_output)
for record in blast_out_parse(handle):
    if record.sseqid not in blast_result:
        blast_result[record.sseqid] = [record]
    else:
        blast_result[record.sseqid].append(record)

total = 0
for key in blast_result:
    total += len(blast_result[key])
    #print key,len(blast_result[key])

#align = Alignment(Gapped(IUPAC.unambiguous_dna,"-"))
#align.add_sequence("Seq_1","ACTCGATC")
#align.add_sequence("Seq_2","ACTCGATC")

for key in blast_result:
    total = 0
    for record in blast_result[key]:
        coverd_length = abs(int(record.send) - int(record.sstart))
        print key,"%s..%s\tlength:%s"%(record.sstart,record.send,record.length),
        print "mapped by:%s"%record.qseqid,
        print "coverage:%.3f"%(coverd_length*1.0/len(ref_seq_dict[key].seq)),
        print "contig_usage:%.3f"%(abs(int(record.qstart)-int(record.qend)*1.0)/\
                len(contig_seq_dict[record.qseqid].seq))
        total += int(coverd_length)
    print key,"total_coverage:%.5f"%(int(total)*1.0/len(ref_seq_dict[key].seq))


