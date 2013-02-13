#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"

from Bio import SeqIO
from bing import gff_parse

class Record():
    pass

ref_file = "/Users/bingwang/zen/yeast_anno_pipe/SGD_data/scer.gff"
contig_file = "/Users/bingwang/zen/yeast_anno_pipe/SeubFM1318_5000_130121.fasta"
annotation_file = "/Users/bingwang/zen/yeast_anno_pipe/seub.annotation_final.txt"
correspondances_file = "/Users/bingwang/zen/yeast_anno_pipe/correspondances.txt"
out_put_file = "/Users/bingwang/zen/yeast_anno_pipe/seubONscer.txt"

numid2scaf = {}
f = open(correspondances_file)
for line in f:
    scaffold,numid,alpheid = line.strip().split("\t")
    numid2scaf[numid] = scaffold

scer2seub = {}
seub2record = {}
for line in open(annotation_file):
    annos = line.rstrip().split("\t")
    seub,start,end,strand,scaf_num,scer = \
            annos[0],annos[2],annos[3],annos[1],annos[5],annos[8].split(" ")[1].split("|")
    record = Record()
    record.id = seub
    record.start = int(start)
    record.end = int(end)
    record.seqid = numid2scaf[scaf_num]
    record.strand = "+" if strand == "1" else "-"
    seub2record[seub] = record
    for gene in scer:
        if gene in scer2seub:
            scer2seub[gene].append(seub)
        else:
            scer2seub[gene] = [seub]

scer_genes = []
for record in gff_parse.gffIterator(open(ref_file)):
    if record.type == "gene" and record.attributes["orf_classification"][0] == "Dubious":
        #dubious.append(record.attributes["ID"][0])
        continue
    elif record.type == "gene" and not record.attributes["ID"][0].startswith("Q"):
        seqid = record.seqid
        start = int(record.start)
        end = int(record.end)
        record.id = record.attributes["ID"][0]
        scer_genes.append(record)

f = open(out_put_file,"w")
for scer in scer_genes:
    seub_ids = scer2seub[scer.id] if scer.id in scer2seub else ""
    line = "%s\t"*9+"%s\n"
    if seub_ids:
        if len(seub_ids) == 1:
            seub = seub2record[seub_ids[0]]
            main_seqid = seub.seqid
            f.write(line%(scer.id,scer.seqid,scer.start,scer.end,scer.strand,\
                    seub.id,seub.seqid,seub.start,seub.end,seub.strand))
        else:
            for seub_id in seub_ids:
                seub = seub2record[seub_id]
                if seub.seqid == main_seqid:
                    f.write(line%(scer.id,scer.seqid,scer.start,scer.end,scer.strand,\
                            seub.id,seub.seqid,seub.start,seub.end,seub.strand))
                else:
                    f.write(("---"+line)%(scer.id,scer.seqid,scer.start,scer.end,scer.strand,\
                            seub.id,seub.seqid,seub.start,seub.end,seub.strand))
    else:
        f.write(line%(scer.id,scer.seqid,scer.start,scer.end,scer.strand,\
                "","","","",""))


