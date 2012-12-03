#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

#data_path = "/Users/kmeihua1017/Documents/G875data/data/"
#out_path = "/Users/kmeihua1017/Documents/G875data/upseq_95/"
data_path = "/Users/bingwang/zen/yeast_anno_pipe/YGOB_dataset/"
out_path = "/Users/bingwang/zen/get_up_out/"



class Gene():
    def __init__(self,gene_name):
        self.name = gene_name        

    def _get_seq(self,n,allow_overlap):
        try: 
            self.up_gene
        except:
            print "%s in sp %s has no up_gene?" % (self.name,self.sp)
        if not allow_overlap and self.up_gene != "":
            start = max(self.start-n,1,self.up_gene.end+1)
        else:
            start = max(self.start-n,1)
        if not allow_overlap and self.down_gene != "":
            end = min(self.down_gene.start-1,self.end+n,len(YGOB[self.sp][self.seqid]))
        else:
            end = min(self.end+n,len(YGOB[self.sp][self.seqid]))
        if self.strand == "+":
            self.up_seq = YGOB[self.sp][self.seqid][start-1:self.start-1]
            self.seq = YGOB[self.sp][self.seqid][self.start-1:self.end]
            self.down_seq = YGOB[self.sp][self.seqid][self.end:end]
        if self.strand == "-":
            self.up_seq = YGOB[self.sp][self.seqid][self.end:end].reverse_complement()
            self.seq = YGOB[self.sp][self.seqid][self.start-1:self.end].reverse_complement()
            self.down_seq = YGOB[self.sp][self.seqid][start-1:self.start-1].reverse_complement()

    def up_stream(self,n,allow_overlap=None):
        if allow_overlap == None:
            allow_overlap = False
        self._get_seq(n,allow_overlap)
        return self.up_seq

def gene_db_init():
    # input pillars.tab and all genome.tab file
    # output gene_db
    # for record in gene_db:
    #     record.name = gene_name
    #     record.homo_list = it's homolog genes
    #     record.sp = it's species
    #     record.cp = it's chromosome copies(valid for WGD species)
    #     record.seqid = chromosome_num
    #     record.strand = "+" or "-"
    #     record.start = ORF start
    #     record.end = ORF end
    #     record.up_gene = up stream gene name
    #     record.down_gene = down stream gene name
    #     record.descrip = description

    gene_db = {}
    pillar_file = data_path+"Pillars.tab"
    pillar_column =["Vpolyspora_1","Tphaffii_1","Tblattae_1","Ndairenensis_1",\
            "Ncastellii_1","Knaganishii_1","Kafricana_1","Cglabrata_1","Suvarum_1",\
            "Skudriavzevii_1","Smikatae_1","Scerevisiae_1","Ancestor","Zrouxii",\
            "Tdelbrueckii","Klactis","Egossypii","Ecymbalariae","Lkluyveri",\
            "Lthermotolerans","Lwaltii","Scerevisiae_2","Smikatae_2",\
            "Skudriavzevii_2","Suvarum_2","Cglabrata_2","Kafricana_2","Knaganishii_2",\
            "Ncastellii_2","Ndairenensis_2","Tblattae_2","Tphaffii_2","Vpolyspora_2"]


    f = open(pillar_file)
    for line in f:
        names = line.strip().split("\t")
        for i,name in enumerate(names):
            if name != "---":
                record = Gene(name)
                record.homo_list = [n for n in names if n != "---"]
                sp_info = pillar_column[i].split("_")
                if len(sp_info) == 1:
                    record.sp,record.cp = sp_info[0],1
                else:
                    record.sp,record.cp = sp_info[0],int(sp_info[1])
                gene_db[name] = record

    for genome in genome_list:
        genome_tab_file = data_path + genome + "_genome.tab"
        gene_db[genome] = []
        f = open(genome_tab_file)
        up = ""
        for line in f:
            name,strand,start,end,YGOB,seqid,short,coord,notes = line.split("\t")
            if name in gene_db and strand == "1":
                gene_db[genome].append(name)
                if up != "" and gene_db[up].seqid == int(seqid):
                    gene_db[up].down_gene = gene_db[name]
                    gene_db[name].up_gene = gene_db[up]
                else:
                    gene_db[name].up_gene = ""
                gene_db[name].strand = "+"
                gene_db[name].start = int(start)
                gene_db[name].end = int(end)
                gene_db[name].seqid = int(seqid)
                gene_db[name].descrip = notes.strip()
                gene_db[name].down_gene = ""
                up = name
        f.close()
        f = open(genome_tab_file)
        up = ""
        for line in f:
            name,strand,start,end,YGOB,seqid,short,coord,notes = line.split("\t")
            if name in gene_db and strand == "0":
                gene_db[genome].append(name)
                if up != "" and gene_db[up].seqid == int(seqid):
                    gene_db[up].down_gene = gene_db[name]
                    gene_db[name].up_gene = gene_db[up]
                else:
                    gene_db[name].up_gene = ""
                gene_db[name].strand = "-"
                gene_db[name].start = int(start)
                gene_db[name].end = int(end)
                gene_db[name].seqid = int(seqid)
                gene_db[name].descrip = notes.strip()
                gene_db[name].down_gene = ""
                up = name
        f.close()
    return gene_db

def load_seq():
    # input all sequence.fsa
    # return ygob dict
    # for sp in ygob:
    #     for seqid in ygob[sp]:
    #         ygob[sp][seqid] = seq
    ygob = {}
    for genome in genome_list:
        ygob[genome] = {}
        fsa_file = data_path + genome + "_sequence.fsa"
        for record in SeqIO.parse(fsa_file,"fasta"):
            seqid = int(record.id.split("_")[-1])
            ygob[genome][seqid] = record.seq
    return ygob

def return_record(name,n):
    if gene_db[name].sp in genome_list:
        id = name
        try: 
            sequence =gene_db[name].up_seq
        except:
            sequence = gene_db[name].up_stream(n)
        length = "up_%d_bps"%(len(sequence))
        chr_cp = "Copy:%d"%(gene_db[name].cp)
        gene_descrip = "Notes:%s"%(gene_db[name].descrip)
        # Define your own descrip here:
        descrip = "\t".join([length,chr_cp,gene_descrip])
    else:
        print name,
        print gene_db[name].sp,
        print " haven't initial"
    if len(sequence) == 0:
        print name,
        right = gene_db[name]
        up = gene_db[name].up_gene
        down = gene_db[name].down_gene
        strand = right.strand
        scaf_len = len(YGOB[gene_db[name].sp][gene_db[name].seqid])
        if strand == "+":
            if up != "" and up.end >= right.start:
                print "up overlap"
            elif up != "" and up.end + 1 == right.start:
                print "up continue"
            elif right.start == 1:
                print "start scaffold head"
            else:
                print name,gene_db[name].sp,gene_db[name].seqid,gene_db[name].strand,gene_db[name].start,gene_db[name].end
                try:
                    print "\t\t",up.name,up.seqid,up.strand,up.start,up.end
                except:
                    print "\t\t\t\"\""
                try:
                    print "\t\t",down.name,down.seqid,down.strand,down.start,down.end
                except:
                    print "\t\t\t\"\""
                print "\t\tscaffold length:",scaf_len
        else:
            if down != "" and down.start <= right.end:
                print "down overlap"
            elif down != "" and right.end + 1 == down.start:
                print "down continue"
            elif right.end == scaf_len:
                print "end at scaffold tail"
            else:
                print name,gene_db[name].sp,gene_db[name].seqid,gene_db[name].strand,gene_db[name].start,gene_db[name].end
                try:
                    print "\t\t",up.name,up.seqid,up.strand,up.start,up.end
                except:
                    print "\t\t\t\"\""
                try:
                    print "\t\t",down.name,down.seqid,down.strand,down.start,down.end
                except:
                    print "\t\t\t\"\""
                print "\t\tscaffold length:", scaf_len
        return ""
    return SeqRecord(sequence,id,name,descrip)

def write_homo_by_sp(name,n):
    # wirte all homo genes, file according to sp
    # write all genes for back ground freq calculation, file according to sp
    # write all genes frequency, file according to sp

    records_for_all = []
    records_for_list = []
    ####         frequency initial      ####
    freq = {}
    sum_chr_1 = 0
    sum_chr_2 = 0
    for chr_1 in "ATCG":
        freq[chr_1] = 0
        for chr_2 in "ATCG":
            freq[chr_1+chr_2] = 0
    ##########################
    for gene_name in gene_db[name]:
        record = return_record(gene_name,n)
        if record != "":
        ####        Calculating Frequency      ####
            for chr_1 in "ATCG":
                count_chr_1 = record.seq.count(chr_1)
                freq[chr_1] += count_chr_1
                sum_chr_1 += count_chr_1
                for chr_2 in "ATCG":
                    count_chr_2 = record.seq.count(chr_1+chr_2)
                    freq[chr_1+chr_2] += count_chr_2
                    sum_chr_2 += count_chr_2
        ################################
            records_for_all.append(record)
            if record.id in all_homo:
                records_for_list.append(record)
    ####       writting Frequency       ####
    freq_name = "%s%s_%d_genes_%d_up_freq.txt" % (out_path,name,len(records_for_all),n)
    f = open(freq_name,"w")
    for chr_1 in "ATCG":
        f.write("%s\t%.3f\n" % (chr_1, freq[chr_1] * 1. / sum_chr_1))
    for chr_1 in "ATCG":
        for chr_2 in "ATCG":
            f.write("%s\t%.3f\n" %(chr_1+chr_2, freq[chr_1+chr_2] * 1. / sum_chr_2))
    ########################################

    file_name = "%s%s_%d_genes_%d_up.fsa" % (out_path,name,len(records_for_all),n)
    SeqIO.write(records_for_all, file_name, "fasta")
    file_name_list = "%s%s_%d_homo_%d_up.fsa" % (out_path,name,len(records_for_list),n)
    SeqIO.write(records_for_list, file_name_list, "fasta")

def write_homo_by_genes(name,n):
    if name in gene_db: 
        gene_list = [homo for homo in gene_db[name].homo_list if gene_db[homo].sp in genome_list]
        write_homo_by_list(gene_list,name,n)
    else:
        print "%s are not in Scer pillar"%(name)

def write_homo_by_list(gene_list,list_name,n):
    records = []
    for gene in gene_list:
        record = return_record(gene,n)
        if record != "":
            records.append(record)
    file_name = "%s%s_%d_gene_%d_up.fsa"%(out_path,list_name,len(gene_list),n)
    SeqIO.write(records, file_name, "fasta")    

##############
#    main    #
##############

gene_list = ["YBR285W","YDL241W","YDR018C","YEL041W","YER121W","YFL052W","YFL054C","YFR017C","YGL117W",\
    "YGL146C","YGR043C","YGR205W","YHR029C","YHR033W","YHR210C","YIL024C","YIR016W","YJR115W",\
    "YLR312C","YML131W","YMR103C","YMR206W","YNL144C","YNR071C","YNR073C","YOR343C","YPL113C",\
    "YPL230W","YPR196W","YDR216W","YDR085C","YMR280C","YFR014C","YOR100C","YPR030W","YHR043C",\
    "YBL043W","YBR033W","YDR516C","YOR383C","YGR243W","YJL221C","YNR002C","YOR178C","YDR009W",\
    "YPL248C","YPR184W","YER054C","YEL011W","YKR058W","YFR015C","YHL032C","YIL155C","YPR005C",\
    "YER062C","YFR053C","YMR011W","YHR092C","YER065C","YPR006C","YMR081C","YKL217W","YCR091W",\
    "YPL054W","YGR244C","YGR289C","YBR298C","YBR299W","YBR297W","YKL093W","YFR030W","YLL061W",\
    "YDL079C","YDL079C","YDR277C","YOL104C","YBR001C","YPL134C","YPL171C","YDR406W","YIL107C",\
    "YAR071W","YBR093C","YBR296C","YLR327C","YBR050C","YDL194W","YER046W","YIL162W","SUC4",\
    "YGL096W","YKR098C","YER098W","YDR247W","YMR104C"]

genome_list = ["Ncastellii","Scerevisiae","Klactis","Smikatae","Skudriavzevii",\
            "Cglabrata","Knaganishii","Kafricana","Lkluyveri","Tdelbrueckii","Egossypii",\
                "Tphaffii","Zrouxii","Suvarum"]

gene_db = gene_db_init()
YGOB = load_seq()

all_homo = []
for gene in gene_list:
    if gene in gene_db:
        all_homo += [homo for homo in gene_db[gene].homo_list if gene_db[homo].sp in genome_list]
    else:
        print "%s not in pillar file"%(gene)
N = 1000

for name in gene_list:
    write_homo_by_genes(name,N)

for name in genome_list:
    write_homo_by_sp(name,N)

write_homo_by_list(all_homo, "all", N)










