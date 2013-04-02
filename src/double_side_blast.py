#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"

from bing import gff_parse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

genome_file_1 = "/Users/bingwang/zen/yeast_anno_pipe/suva.genome.txt"
anno_file_1 = "/Users/bingwang/zen/yeast_anno_pipe/suva.annotation_final.txt"
genome_file_2 = "/Users/bingwang/zen/yeast_anno_pipe/seub3000.genome.txt"
anno_file_2 = "/Users/bingwang/zen/yeast_anno_pipe/seub3000.annotation_final.txt"

nt_file_1 = "/Users/bingwang/zen/yeast_anno_pipe/suva.nt.fsa"
nt_file_2 = "/Users/bingwang/zen/yeast_anno_pipe/seub.nt.fsa"

final_out = "/Users/bingwang/zen/yeast_anno_pipe/suva_seub_pair.tab"
def YGAPIterator(handle):
    #for record in YGAPIterator(handle):
    #    record.name = gene name
    #    record.strand = strand
    #    record.start = lowest coordinate
    #    record.end = highest coordinate
    #    record.ygob = YGOB browser status (all ON)
    #    record.seqid = chromosome/scaffold num
    #    record.shortname = gene short name
    #    record.coord = gene coordinates (including introns)
    #    record.orth = Ancestral Scer1|Scer2
    #    record.type = gene type(empty,PROTEIN,TKP/TY)
    #    record.pillar = YGOB pillar num
    #    record.tag = annotation tag(Cor/Fra/Man/NNN)
    #        Cor: successful frameshift correction
    #        Fra: should have frameshift but cannot correct
    #        Man: need to be manual attention(may untranslatable, overlap, strange length)
    #        NNN: N at begins or ends of the coding region
    #    record.anno = Annotation route(simple/multi/GETORF/DOGS)
    #        Simple: easily assigned to YGOB pillar
    #        Multi: gene has multigen family and its pillar may be unreliable
    #        GETORF: generate by GETORF
    #        DOGS: generate by DOGS
    class record():
        def __init__(self):
            pass
    while 1:
        line = handle.readline()
        if line.count("..") > 1:
            print "intron @:%s"%line
        if line.startswith("#"): continue
        if not line: return
        r = record()
        (r.name,r.strand,r.start,r.end,r.ygob,r.seqid,r.shortname,\
            r.coord,r.orth,r.type,r.pillar,r.tag,r.anno) = line.split("\t")
        r.anno = r.anno.strip()
        r.strand = "+" if r.strand == "1" else "-"
        #r.orth = [r.orth.split(" ")[0]] + r.orth.split(" ")[1].split("|")
        yield r

def write_nt_file(genome_file,anno_file,nt_file):
    sp2sequence = {}
    sp_chr2seq = {}
    sp_records = []
    for record in SeqIO.parse(genome_file,"fasta"):
        chr_id = record.id.split("_")[-1]
        sp_chr2seq[chr_id] = record.seq
    for record in YGAPIterator(open(anno_file)):
        start = int(record.start)
        end = int(record.end)
        if record.strand == "+":
            sp_records.append(SeqRecord(sp_chr2seq[record.seqid][start-1:end],\
                    id=record.name, description=""))
        else:
            sp_records.append(SeqRecord(sp_chr2seq[record.seqid][start-1:end].reverse_complement(),\
                    id=record.name, description=""))
    SeqIO.write(sp_records,nt_file,"fasta")
        
#1. generate two nucle tide sequence file
write_nt_file(genome_file_1,anno_file_1,nt_file_1)
write_nt_file(genome_file_2,anno_file_2,nt_file_2)

#2. run double side balst
def double_blast(nt_file_a,nt_file_b):
    #first make 2 db
    #use blastn do the blast
    #generate 2 blast results
    def make_db(nt_file):
        os.system("makeblastdb -in %s -dbtype nucl -out %s"%(nt_file, nt_file+".db"))
    def run_blast(db,nt_file,out_file):
        os.system("blastn -evalue 0.00001 -max_target_seqs 1 "+\
            "-db %s -query %s -out %s "%(db,nt_file,out_file) +\
            "-outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart " +\
            "qend sstart send \"")

    make_db(nt_file_a)
    make_db(nt_file_b)
    run_blast(nt_file_a+".db",nt_file_b,nt_file_b+".out")
    run_blast(nt_file_b+".db",nt_file_a,nt_file_a+".out")
    os.system("rm %s"%(nt_file_a+".db*"))
    os.system("rm %s"%(nt_file_b+".db*"))

double_blast(nt_file_1,nt_file_2)

#3. get 1:1 homolog relation
sp12sp2 = {}
for line in open(nt_file_1+".out"):
    sp1,sp2 = line.split("\t")[:2]
    sp12sp2[sp1] = sp2

sp22sp1 = {}
for line in open(nt_file_2+".out"):
    sp2,sp1 = line.split("\t")[:2]
    sp22sp1[sp2] = sp1

f = open(final_out,"w")
for sp1 in sp12sp2.keys():
    sp2 = sp12sp2[sp1]
    if sp2 in sp22sp1 and sp22sp1[sp2] == sp1:
        f.write("%s\t%s\n"%(sp1,sp2))

