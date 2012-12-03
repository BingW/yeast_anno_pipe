#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"
# dependents: YGOB/
# biopython
# bing

from Bio import SeqIO
YGOB_path="/Users/bingwang/zen/yeast_anno_pipe/YGOB_dataset/"
out_path="/Users/bingwang/zen/yeast_anno_pipe/output/"

class Yeast_Gene:
    #self.length
    def __init__(self,gene_name,chr_id=None):
        self.name = gene_name
        self.chr_id = chr_id

    def get_pos(self,gene_pos):
        self.pos_str = gene_pos
        self.strand = "-" if "complement" in gene_pos else "+"
        self.intron = True if gene_pos.count("..") > 1 else False
        pos_string = gene_pos[gene_pos.find("(")+1:gene_pos.rfind(")")]
        if self.intron:
            self.pos = []
            pairs = [a for a in pos_string.split(",")]
            for pair in pairs:
                self.pos.append((int(pair.split("..")[0]),int(pair.split("..")[1])))
        else:
            a,b = pos_string.split("..")
            self.pos = [int(a),int(b)] 

def intron_pillar(): 
    gene_list = []
    for line in open(out_path+"YGOB_intron.fsa"):
        if line.startswith(">"):
            gene_list.append(line.split("\t")[0][1:])

    f = open(out_path+"intron_pillar.tab", "w")
    for line in open(YGOB_path+"Pillars.tab"):
        flag = False
        names = line.split("\t")
        for name in names:
            if name != "---" and name in gene_list:
                flag = True
                break
        if flag:
            record = []
            for name in names:
                if name == "---" or name in gene_list:
                    record.append(name)
                else:
                    record.append("+++")
            f.write("\t".join(record)+"\n")

name2sp = {}
seq_dict = {}
for sp_name in ["Vpolyspora","Tphaffii","Tblattae","Ndairenensis",\
            "Ncastellii","Knaganishii","Kafricana","Cglabrata","Suvarum",\
            "Skudriavzevii","Smikatae","Scerevisiae","Zrouxii",\
            "Tdelbrueckii","Klactis","Egossypii","Ecymbalariae","Lkluyveri",\
            "Lthermotolerans","Lwaltii"]:
    for record in SeqIO.parse(YGOB_path + sp_name + "_sequence.fsa","fasta"):
        seq_dict[sp_name+"_"+record.id[record.id.find("_"):].replace("_","")] = record.seq
    for line in open(YGOB_path + sp_name + "_genome.tab"):
        gene_name = line.split("\t")[0]
        name2sp[gene_name] = sp_name

Gene_dict = {}
for record in SeqIO.parse(YGOB_path+"AA.fsa","fasta"):
    gene_name,chr_name,gene_pos,length = record.description.split(" ")[:4]
    if gene_pos.count("..") > 1:
        sp = name2sp[gene_name]
        Gene_dict[gene_name] = Yeast_Gene(gene_name,chr_name)
        Gene_dict[gene_name].get_pos(gene_pos)
        Gene_dict[gene_name].aa = record.seq
        Gene_dict[gene_name].sp = sp
        try:
            pos_start = Gene_dict[gene_name].pos[0][0]
        except:
            print gene_name
            assert 1 == 0
        pos_end = Gene_dict[gene_name].pos[-1][1]
        if Gene_dict[gene_name].strand == "+":
            Gene_dict[gene_name].full_nt = seq_dict[sp+"_"+chr_name][pos_start-1:pos_end]
        else:
            Gene_dict[gene_name].full_nt = \
            seq_dict[sp+"_"+chr_name][pos_start-1:pos_end].reverse_complement()

f = open(out_path+"YGOB_intron.fsa","w")
for gene_name in Gene_dict:
    sp = Gene_dict[gene_name].sp
    chr_id = Gene_dict[gene_name].chr_id
    starnd = Gene_dict[gene_name].strand
    pos = Gene_dict[gene_name].pos_str
    aa = str(Gene_dict[gene_name].aa)
    nt_seq = str(Gene_dict[gene_name].full_nt)
    f.write(">%s\t%s\t%s\t%s\t%s\t%s\n" % (gene_name,sp,chr_id,starnd,pos,aa))
    f.write("%s\n" % nt_seq)



