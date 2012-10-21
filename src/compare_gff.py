#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"
# dependents: 
# bing

from bing import gff_parse
deven_gff = open("/Users/bingwang/zen/yeast_anno_pipe/SSS_dataset/Skud.gff")
augustus_gff = open("/Users/bingwang/zen/yeast_anno_pipe/SSS_dataset/skud.augustus.gff")
ygap_anno = open("/Users/bingwang/zen/yeast_anno_pipe/SSS_dataset/skud.YGAP.txt")
#pc_gff = open("/Users/bingwang/zen/yeast_anno_pipe/SSS_dataset/Skud.gff")


def YGAPIterator(handle):
    class record():
        def __init__(self):
            pass
    while 1:
        line = handle.readline()
        if line.startswith("#"): continue
        if not line: return
        r= record()
        (r.name,r.strand,r.start,r.end,r.ygob,r.seqid,r.shortname,\
            r.coord,r.orth,r.type,r.pillar,r.tag,r.anno) = line.split("\t")
        r.strand = "+" if r.strand == "1" else "-"
        yield r
    
    
class list(list):
    def __init__(self):
        self.name = ""

deven = list()
deven.name = "deven"
for record in gff_parse.gffIterator(deven_gff):
    deven.append(record)

augustus = list() 
augustus.name = "augustus"
for record in gff_parse.gffIterator(augustus_gff):
    record.seqid = record.seqid[5:]
    augustus.append(record)

ygap = list()
ygap.name = "ygap"
for record in YGAPIterator(ygap_anno):
    ygap.append(record)

def compare_gff(template,*gffs):
    temp = [read.seqid+"_"+read.start+"_"+read.end+"_"+read.strand for \
            read in template if read.type == "CDS"]
    other = {}
    for gff in gffs:
        other[gff.name] = [read.seqid+"_"+read.start+"_"+read.end+"_"+read.strand for \
            read in gff if read.type == "gene" or read.type == "PROTEIN"]
    sn = {}
    sp = {}
    for gff in gffs:
        TP = len([pos for pos in other[gff.name] if\
                pos in temp])
        sn[gff.name] = TP*1.0/len(temp)
        sp[gff.name] = TP*1.0/len(other[gff.name])
        
    print "sn:",sn
    print "sp:",sp

compare_gff(deven,augustus,ygap)

