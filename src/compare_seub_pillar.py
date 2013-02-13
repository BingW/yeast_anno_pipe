#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"

Out_Path = "/Users/bingwang/zen/yeast_anno_pipe/output/YGAPvsDevin/"

from bing import gff_parse
from Bio import SeqIO
import os

pillar_file = "/Users/bingwang/zen/yeast_anno_pipe/YGOB_dataset/Pillars.tab"
seub_annotation = "/Users/bingwang/zen/yeast_anno_pipe/YGAP_predict/seub.annotation_final.txt"
seub_pillar = "/Users/bingwang/zen/yeast_anno_pipe/YGAP_predict/seub.pillar.tab"

pillar_column =["Vpolyspora_1","Tphaffii_1","Tblattae_1","Ndairenensis_1",\
        "Ncastellii_1","Knaganishii_1","Kafricana_1","Cglabrata_1","Suvarum_1",\
        "Skudriavzevii_1","Smikatae_1","Scerevisiae_1","Ancestor","Zrouxii",\
        "Tdelbrueckii","Klactis","Egossypii","Ecymbalariae","Lkluyveri",\
        "Lthermotolerans","Lwaltii","Scerevisiae_2","Smikatae_2",\
        "Skudriavzevii_2","Suvarum_2","Cglabrata_2","Kafricana_2","Knaganishii_2",\
        "Ncastellii_2","Ndairenensis_2","Tblattae_2","Tphaffii_2","Vpolyspora_2"]

scer2seub = {}
seub2anc = {}
seub_synteny = []
for line in open(seub_annotation):
    annos = line.rstrip().split("\t")
    seub,anc,homo = annos[0],annos[8].split(" ")[0],annos[8].split(" ")[1].split("|")
    if not seub.startswith("SEUB"):
        continue
    seub_synteny.append(seub)
    if anc:
        seub2anc[seub] = anc
    else:
        seub2anc[seub] = "---"
    for gene in homo:
        if gene: 
            scer2seub[gene] = seub

all_synteny = []
seub_has_homo = []
f = open(seub_pillar,"w")
for line in open(pillar_file):
    if not "trna" in line and not "snR" in line and not "Scer_" in line:
        genes = line.rstrip().split("\t")
        scer_1,scer_2,anc,suva_1,suva_2 = genes[11],genes[-12],genes[12],genes[8],genes[-7]
        if scer_1 in scer2seub:
            seub_1 = scer2seub[scer_1]
            if anc != "---":
                if not seub2anc[seub_1] == anc:
                    print "!conflict:%s is annotated has anc: %s and scer: %s while %s has anc homo %s"\
                            %(seub_1,seub2anc[seub_1],scer_1,scer_1,anc)
        else:
            seub_1 = "---"
        if scer_2 in scer2seub:
            seub_2 = scer2seub[scer_2]
            if anc != "---":
                if not seub2anc[seub_2] == anc:
                    print "!conflict:%s is annotated has anc: %s and scer: %s while %s has anc homo %s"\
                            %(seub_2,seub2anc[seub_2],scer_2,scer_2,anc)
        else:
            seub_2 = "---"
        all_synteny.append((scer_1,scer_2,anc,suva_1,suva_2,seub_1,seub_2))
        seub_has_homo += [seub_1,seub_2]

for line in all_synteny:
    if line != ("---","---","---","---","---","---","---"):
        print "\t".join(line)
        pass

seub_only = [gene for gene in seub_synteny if gene not in seub_has_homo]
assert len(seub_only) == len(set(seub_only))
for gene in seub_only:
    gene_up = seub_synteny[seub_synteny.index(gene)-1]
    gene_down = seub_synteny[seub_synteny.index(gene)+1]
    if gene_up in seub_has_homo and gene_down in seub_has_homo:
        seub_has_homo.append(gene)
seub_only = [gene for gene in seub_synteny if gene not in seub_has_homo]
print len(seub_only)


'''
scer_only = set([g for g in scer2seub if not scer2seub[g]])
anc_only = set([g for g in anc2seub if not anc2seub[g]])

print len(seub_type["+Anc+Scer"])
print len(seub_type["+Anc-Scer"])
print len(seub_type["-Anc-Scer"])
print len(seub_type["-Anc+Scer"])
print len(scer_only)
print len(anc_only)

for gene in scer_only:
    i =sdacer_synteny_1.index(gene)
    gene_up = None if i == 0 else scer_synteny[i-1]
    gene_down = None if i == len(scer_synteny)-1 else scer_synteny[i+1]
    if gene_up and scer2seub[gene_up]:
        if gene_down and scer2seub[gene_down]:
            scer_only[gene] = "has_both"
        else:
            scer_only[gene] = "has_up"
    else:
        if gene_down and scer2seub[gene_down]:
            scer_only[gene] = "has_down"
        else:
            scer_only[gene] = "has_none"

print len([g for g in scer_only if scer_only[g] == "has_both"])
print len([g for g in scer_only if scer_only[g] == "has_up"])
print len([g for g in scer_only if scer_only[g] == "has_down"])
print len([g for g in scer_only if scer_only[g] == "has_none"])
'''
