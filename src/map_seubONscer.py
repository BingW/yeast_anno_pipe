#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"

from Bio import SeqIO
from bing import gff_parse
from collections import Counter

#Prepare: empty class Record()
class Record():
    pass

#Prepare: class Stack()
class Stack(list):
    def __init__(self,n):
        self.limit = n
    def push(self,value):
        if value == "" and len(self) == 0:
            pass
        else:
            self.append(value)
        if len(self) > self.limit:
            self.pop(0)
        if self.count("") >= 2:
            self.empty()
    def empty(self):
        del(self[:])

#Prepare: Input files
sp = "seub3000"
ref_file = "/Users/bingwang/zen/yeast_anno_pipe/SGD_data/scer.gff"
contig_file = "/Users/bingwang/zen/yeast_anno_pipe/%s.genome.txt.fasta"%(sp)
annotation_file = "/Users/bingwang/zen/yeast_anno_pipe/%s.annotation_final.txt"%(sp)
correspondances_file = "/Users/bingwang/zen/yeast_anno_pipe/%s.correspondances.txt"%(sp)

#Prepare: Output file
out_put_file = "/Users/bingwang/zen/yeast_anno_pipe/%sONscer.txt"%(sp)

#Prepare: numid2scaf dict
numid2scaf = {}
scaffold_list = []
f = open(correspondances_file)
for line in f:
    scaffold,numid,alpheid = line.strip().split("\t")
    numid2scaf[numid] = scaffold
    scaffold_list.append(scaffold)

#Prepare: scer_genes synteny list
#Prepare: scer2record dict
scer_genes = []
scer2record = {}
for record in gff_parse.gffIterator(open(ref_file)):
    if record.type == "gene" and record.attributes["orf_classification"][0] == "Dubious":
        #dubious.append(record.attributes["ID"][0])
        continue
    elif record.type == "gene" and not record.attributes["ID"][0].startswith("Q"):
        seqid = record.seqid
        start = int(record.start)
        end = int(record.end)
        record.id = record.attributes["ID"][0]
        scer_genes.append(record.id)
        scer2record[record.id] = record

#Prepare: scer2seub and seub2scer dicts
#Prepare: seub2record dict
scer2seub = {}
seub2scer = {}
seub2record = {}
seub_genes = []
seub_scaf2genes = {}
#seub_scaf2scer_chr = {}
for line in open(annotation_file):
    annos = line.rstrip().split("\t")
    seub,start,end,strand,scaf_num,scer = \
            annos[0],annos[2],annos[3],annos[1],annos[5],annos[8].split(" ")[1].split("|")
    record = Record()
    record.id = seub
    record.start = int(start)
    record.end = int(end)
    record.seqid = numid2scaf[scaf_num]
    record.scer = scer
    #if record.seqid not in seub_scaf2scer_chr:
    #    seub_scaf2scer_chr[record.seqid] = []
    record.strand = "+" if strand == "1" else "-"
    if record.seqid in seub_scaf2genes:
        seub_scaf2genes[record.seqid].append(record.id)
    else:
        seub_scaf2genes[record.seqid]=[record.id]
    seub_genes.append(record.id)
    seub2record[seub] = record
    seub2scer[record.id] = []
    for gene in scer:
        if gene in scer2record:
            seub2scer[seub].append(gene)
            #seub_scaf2scer_chr[record.seqid].append(scer2record[gene].seqid)
            if gene in scer2seub:
                scer2seub[gene].append(seub)
            else:
                scer2seub[gene] = [seub]

#for item in seub_scaf2scer_chr:
#    try:
#        seub_scaf2scer_chr[item] = Counter(seub_scaf2scer_chr[item]).most_common(1)[0][0]
#    except:
#        print item,Counter(seub_scaf2scer_chr[item]).most_common(1)

scer121seub = {}
scer121seub["00"] = "00"
seub121scer = {}
seub121scer["00"] = "00"
#Main step 1. find simplist 1 to 1 relation.
for scer in scer2seub:
    if len(scer2seub[scer]) == 1:
        seub = scer2seub[scer][0]
        if len(seub2scer[seub]) == 1 and scer == seub2scer[seub][0]:
            scer121seub[scer] = seub
            seub121scer[seub] = scer

for k in range(2):
    loop = 1
    #Main step 2. find scer 1: seub n relation.
    for scer in scer2seub:
        if len(scer2seub[scer]) > 0 and scer not in scer121seub:
            i = j = scer_genes.index(scer)
            up_scer = scer_genes[i-1]
            down_scer = scer_genes[j+1]
            while up_scer not in scer121seub:
                try:
                    i -= 1
                    up_scer = scer_genes[i-1]
                except:
                    up_scer = "00"
                    break
            while down_scer not in scer121seub:
                try:
                    j += 1
                    down_scer = scer_genes[j+1]
                except:
                    down_scer = "00"
                    break
            up_seub_seqid = scer121seub[up_scer].split("0")[1]
            down_seub_seqid = scer121seub[down_scer].split("0")[1]
            flag = 0
            for gene in scer2seub[scer]:
                if gene.split("0")[1] == up_seub_seqid or gene.split("0")[1] == down_seub_seqid:
                    seub = gene
                    flag += 1
            if flag == 1:
                scer121seub[scer] = seub
                seub121scer[seub] = scer
                if k == loop:
                    print "up:", up_scer, scer121seub[up_scer]
                    print "down:", down_scer, scer121seub[down_scer]
                    print "Set: ", scer, seub

    #Main step 3. find seub 1: scer n relation.
    for seub in seub2scer:
        if len(seub2scer[seub]) > 0 and seub not in seub121scer:
            i = j = seub_genes.index(seub)
            up_seub = seub_genes[i-1]
            down_seub = seub_genes[j+1]
            while up_seub not in seub121scer:
                try:
                    i -= 1
                    up_seub = seub_genes[i-1]
                except:
                    up_seub = "00"
                    break
            while down_seub not in seub121scer:
                try:
                    j += 1
                    down_seub = seub_genes[j+1]
                except:
                    down_seub = "00"
                    break
            up_scer_seqid = seub121scer[up_seub][1:3]
            down_scer_seqid = seub121scer[down_seub][1:3]
            flag = 0
            for gene in seub2scer[seub]:
                if gene[1:3] == up_scer_seqid or gene[1:3] == down_scer_seqid:
                    scer = gene
                    flag += 1
            if flag == 1:
                scer121seub[scer] = seub
                seub121scer[seub] = scer
                if k == loop:
                    print "up:", up_seub, seub121scer[up_seub]
                    print "down:", down_seub, seub121scer[down_seub]
                    print "Set: ", seub, scer 
            elif flag == 2 and seub2scer[seub][0][1:2] == seub2scer[seub][1][1:2] == up_scer_seqid == down_scer_seqid:
                print "Hello"
                seub_0_pos = int(seub2scer[seub][0][3:6])
                seub_1_pos = int(seub2scer[seub][1][3:6])
                up_pos = int(seub121scer[up_seub][3:6])
                down_pos = int(seub121scer[down_seub][3:6])
                if down_pos < seub_0_pos < up_pos or down_pos > seub_0_pos > up_pos:
                    scer = seub2scer[seub][0]
                    scer121seub[scer] = seub
                    seub121scer[seub] = scer
                else:
                    scer = seub2scer[seub][1]
                    scer121seub[scer] = seub
                    seub121scer[seub] = scer
                    if k == loop:
                        print "up:", up_seub, seub121scer[up_seub]
                        print "down:", down_seub, seub121scer[down_seub]
                        print "Set: ", seub, scer 

# make content use scer as reference
content = []
line_format = "%s\t"*9+"%s"
for scer in scer_genes:
    scer_r = scer2record[scer]
    seub = scer121seub[scer] if scer in scer121seub else ""
    if seub:
        seub_r = seub2record[seub]
        content.append(line_format%(scer_r.id,scer_r.seqid,scer_r.start,scer_r.end,scer_r.strand,\
                seub_r.id,seub_r.seqid,seub_r.start,seub_r.end,seub_r.strand))
    else:
        content.append(line_format%(scer_r.id,scer_r.seqid,scer_r.start,scer_r.end,scer_r.strand,\
                "","","","",""))

# calculate scaffold order,occurence and temp count
scaffold_order = [""]
scaffold_occurence = [] #nuber of genes in one scaf
occurence = 0
for line in content:
    now_scaffold = line.split("\t")[6]
    if now_scaffold:
        if now_scaffold != scaffold_order[-1]:
            scaffold_order.append(now_scaffold)
            scaffold_occurence.append(occurence)
            occurence = 1
        else:
            occurence += 1
scaffold_occurence.append(occurence)
scaffold_order.pop(0) #delete the first item which is ""
scaffold_occurence.pop(0) #delete the first item which is occurence of first ""
scaffold_count = Counter(scaffold_order) #one block count as 1

# This step aim at delete saffold interperation
scaffold_to_del = [scaf for i,scaf in enumerate(scaffold_order) if scaffold_occurence[i] == 1 and scaffold_order[i-1] == scaffold_order[i+1]]
print scaffold_to_del
for i,line in enumerate(content):
    if line.split("\t")[6] in scaffold_to_del:
        j = i-1
        while content[j].split("\t")[6] == "":
            j -= 1
        k = i+1
        while content[k].split("\t")[6] == "":
            k += 1
        if content[j].split("\t")[6] == content[k].split("\t")[6] != line.split("\t")[6]:
            del seub121scer[line.split("\t")[5]]
            content[i] = "\t".join(content[i].split("\t")[:5]) + "\t"*5
            print line,"->",content[i]

#print scaffold cover region > 2
print "These scaffold occured at 2 place:"
for i,scaf in enumerate(scaffold_order):
    if scaffold_order[i+1:].count(scaf) == 1:
        scaf_1_index = i
        scaf_2_index = i+1+scaffold_order[i+1:].index(scaf)
        if scaffold_occurence[scaf_1_index] > 1 and \
            scaffold_occurence[scaf_2_index] > 1 and \
            scaffold_occurence[scaf_1_index+1] > 1:
            print "\t",scaf

#print scaffold inerperation
print "There scaffold are abnormol interupted:"
for i,scaf in enumerate(scaffold_order):
    if scaffold_order[i+1:].count(scaf) == 1:
        scaf_1_index = i
        scaf_2_index = i+1+scaffold_order[i+1:].index(scaf)
        if scaffold_occurence[scaf_1_index] > 1 and \
            scaffold_occurence[scaf_2_index] > 1 and \
            scaffold_occurence[scaf_1_index+1] == 1:
            print "\t",scaf

#print single occurence scaffold
print "single occurence scaffold:"
for i,scaf in enumerate(scaffold_order):
    if scaffold_order.count(scaf) == 1 and scaffold_occurence[i] == 1:
        print "\t",scaf

#print no occurence scaffold
print "No occurence scaffold:"
for scaf in scaffold_list:
    if scaf not in scaffold_order:
        print "\t",scaf

single_occur_scaf = [scaf for scaf,occur in zip(scaffold_order,scaffold_occurence) if occur == 1 and scaffold_count[scaf] == 1]
no_occur_scaf = [scaf for scaf in scaffold_list if scaf not in scaffold_order]

#make content insert unmapped seub
#1. no homolog seub genes as a list
nohomo_seub = [seub for seub in seub2record.keys() if seub not in seub121scer]
print "total Seub sp genes:",len(nohomo_seub)

#make new content
new_content = []
for i,line in enumerate(content):
    new_content.append(line)
    #2. read content continuesly get this seub and get the nearlest two seubs if possible
    seub,scaf = line.split("\t")[5:7]
    if seub:
        prev_seub = ""
        j = i - 1
        while not prev_seub and j:
            prev_seub,prev_scaf = content[j].split("\t")[5:7]
            j -= 1
        next_seub = ""
        j = i + 1
        while not next_seub and j<len(content):
            next_seub,next_scaf = content[j].split("\t")[5:7]
            j += 1

        seub_pos = int(seub[seub[5:].index("0")+5:])
        if prev_scaf == scaf:
            prev_seub_pos = int(prev_seub[prev_seub[5:].index("0")+5:])
            if abs(prev_seub_pos - seub_pos) > 10:
            #3. if at same scaf but not continues, get direction, try find 'next seub' next to 'previous seub'
                suppose_pos = prev_seub_pos - 10
                while prev_seub_pos > suppose_pos > seub_pos:
                    range_scaf = prev_seub[prev_seub[5:].index("0")+5:]
                    suppose_seub = prev_seub.replace(str(range_scaf),"0"*(len(range_scaf)-len(str(suppose_pos)))+str(suppose_pos))
                    if suppose_seub in nohomo_seub:
                        new_seub = seub2record[suppose_seub]
                        new_line = line_format%("","","","","",new_seub.id,new_seub.seqid,new_seub.start,new_seub.end,new_seub.strand)
                        new_content.insert(-1,new_line)
                        nohomo_seub.remove(suppose_seub)
                        print "GapInsert: %s insert between %s and %s"%(suppose_seub,prev_seub,seub)
                    suppose_pos -= 10
                suppose_pos = prev_seub_pos + 10
                while prev_seub_pos < suppose_pos < seub_pos:
                    range_scaf = prev_seub[prev_seub[5:].index("0")+5:]
                    suppose_seub = prev_seub.replace(str(range_scaf),"0"*(len(range_scaf)-len(str(suppose_pos)))+str(suppose_pos))
                    if suppose_seub in nohomo_seub:
                        new_seub = seub2record[suppose_seub]
                        new_line = line_format%("","","","","",new_seub.id,new_seub.seqid,new_seub.start,new_seub.end,new_seub.strand)
                        new_content.insert(-1,new_line)
                        nohomo_seub.remove(suppose_seub)
                        print "GapInsert: %s insert between %s and %s"%(suppose_seub,prev_seub,seub)
                    suppose_pos += 10
            if next_scaf != scaf:
            #4. deal with scaf at the tail of one scaf
                suppose_pos = seub_pos + 10
                temp_pos = seub_pos
                while suppose_pos > seub_pos > prev_seub_pos and abs(suppose_pos - temp_pos) < 30:
                    range_scaf = seub[seub[5:].index("0")+5:]
                    suppose_seub = seub.replace(str(range_scaf),"0"*(len(range_scaf)-len(str(suppose_pos)))+str(suppose_pos))
                    if suppose_seub in nohomo_seub:
                        new_seub = seub2record[suppose_seub]
                        new_line = line_format%("","","","","",new_seub.id,new_seub.seqid,new_seub.start,new_seub.end,new_seub.strand)
                        new_content.append(new_line)
                        nohomo_seub.remove(suppose_seub)
                        temp_pos = suppose_pos
                        print "TailAppend: %s append after %s"%(suppose_seub,seub)
                    suppose_pos += 10
                suppose_pos = seub_pos - 10
                while suppose_pos < seub_pos < prev_seub_pos and abs(suppose_pos - temp_pos) < 30:
                    range_scaf = seub[seub[5:].index("0")+5:]
                    suppose_seub = seub.replace(str(range_scaf),"0"*(len(range_scaf)-len(str(suppose_pos)))+str(suppose_pos))
                    if suppose_seub in nohomo_seub:
                        new_seub = seub2record[suppose_seub]
                        new_line = line_format%("","","","","",new_seub.id,new_seub.seqid,new_seub.start,new_seub.end,new_seub.strand)
                        new_content.append(new_line)
                        nohomo_seub.remove(suppose_seub)
                        temp_pos = suppose_pos
                        print "TailAppend: %s append after %s"%(suppose_seub,seub)
                    suppose_pos -= 10
        elif next_scaf == scaf:
            #5. deal with scaf at the head of one scaf
            next_seub_pos = int(next_seub[next_seub[5:].index("0")+5:])
            suppose_pos = seub_pos + 10
            temp_pos = seub_pos
            insert_pos = -1
            while suppose_pos > seub_pos > next_seub_pos and abs(suppose_pos - temp_pos) < 100:
                range_scaf = seub[seub[5:].index("0")+5:]
                suppose_seub = seub.replace(str(range_scaf),"0"*(len(range_scaf)-len(str(suppose_pos)))+str(suppose_pos))
                if suppose_seub in nohomo_seub:
                    new_seub = seub2record[suppose_seub]
                    new_line = line_format%("","","","","",new_seub.id,new_seub.seqid,new_seub.start,new_seub.end,new_seub.strand)
                    new_content.insert(insert_pos,new_line)
                    nohomo_seub.remove(suppose_seub)
                    temp_pos = suppose_pos
                    insert_pos -= 1
                    print "HeadInsert: %s insert before %s"%(suppose_seub,seub)
                suppose_pos += 10
            suppose_pos = seub_pos - 10
            while suppose_pos < seub_pos < next_seub_pos and abs(suppose_pos - temp_pos) < 100:
                range_scaf = seub[seub[5:].index("0")+5:]
                suppose_seub = seub.replace(str(range_scaf),"0"*(len(range_scaf)-len(str(suppose_pos)))+str(suppose_pos))
                if suppose_seub in nohomo_seub:
                    new_seub = seub2record[suppose_seub]
                    new_line = line_format%("","","","","",new_seub.id,new_seub.seqid,new_seub.start,new_seub.end,new_seub.strand)
                    new_content.insert(insert_pos,new_line)
                    nohomo_seub.remove(suppose_seub)
                    temp_pos = suppose_pos
                    insert_pos -= 1
                    print "HeadInsert: %s insert before %s"%(suppose_seub,seub)
                suppose_pos -= 10

#exclude single_occur_scaf and no_occur_scaf
nohomo_seub = [seub for seub in nohomo_seub if seub2record[seub].seqid not in single_occur_scaf + no_occur_scaf]
print nohomo_seub
open(out_put_file,"w").write("\n".join(new_content))
        

'''
        if scaf == next_scaf and scaf in has_nohomo_scaf:
            if int(seub[seub[5:].index("0")+5:]) > int(next_seub[next_seub[5:].index("0")+5:]):
                range_small,range_big = int(next_seub[next_seub[5:].index("0")+5:]),int(seub[seub[5:].index("0")+5:])
                direction = 1 #seub > next_seub
            else:
                range_small,range_big = int(seub[seub[5:].index("0")+5:]),int(next_seub[next_seub[5:].index("0")+5:])
                direction = 0 #seub < next_seub
            for seq_num in range(range_small+10,range_big,10):
                range_scaf = seub[seub[5:].index("0")+5:]
                insert_seub = seub.replace(str(range_scaf),"0"*(len(range_scaf)-len(str(seq_num)))+str(seq_num))
                print range_small,range_big
                print insert_seub,
                if insert_seub in nohomo_seub:
                    print "^_^"
            #if - > 10:
            #    print range_small,range_big,range_big-range_small
            #    for insert_seub in nohomo_seub:
            #        if seub2record[insert_seub].seqid == scaf:
            #            print insert_seub,int(insert_seub[insert_seub[5:].index("0")+5:])
    #for seub in nohomo_seub:
    #    print seub,seub.split("0")[1],int(seub[seub[5:].index("0")+5:]),seub2record[seub].seqid,seub2record[seub].start,seub2record[seub].end,seub2record[seub].strand
'''




'''
dubious_homo_seub = [seub for seub in seub2record.keys() if seub2record[seub].scer != [""] and seub not in seub121scer]
all_seub = [seub for seub in seub2record.keys()]
print sorted(all_seub)
'''
'''
for seub in nohomo_seub:
    seub_pos = seub[seub[5:].index("0")+5:]
    small_pos = int(seub[seub[5:].index("0")+5:])-10
    small_nearby_seub = seub.replace(str(seub_pos),"0"*(len(seub_pos)-len(str(small_pos)))+str(small_pos))
    big_pos = int(seub[seub[5:].index("0")+5:])+10
    big_nearby_seub = seub.replace(str(seub_pos),"0"*(len(seub_pos)-len(str(big_pos)))+str(big_pos))
    small_nearby_seub,big_nearby_seub
    small_nearby_seub_index = [i for i,line in enumerate(content) if small_nearby_seub in line]
    if small_nearby_seub_index:
        small_nearby_seub_index = small_nearby_seub_index[0]
         content[small_nearby_seub_index-1].split("\t")[5]
        print content[small_nearby_seub_index].split("\t")[5]
        print content[small_nearby_seub_index+1].split("\t")[5]
    [j for j,line in enumerate(content) if big_nearby_seub in line]
    #break
'''
'''
has_nohomo_scaf = set([seub2record[seub].seqid for seub in nohomo_seub])
for i,line in enumerate(content[:-1]):
    seub,scaf = line.split("\t")[5:7]
    next_seub,next_scaf = content[i+1].split("\t")[5:7]
    if scaf == next_scaf and scaf in has_nohomo_scaf:
        if int(seub[seub[5:].index("0")+5:]) > int(next_seub[next_seub[5:].index("0")+5:]):
            range_small,range_big = int(next_seub[next_seub[5:].index("0")+5:]),int(seub[seub[5:].index("0")+5:])
            direction = 1 #seub > next_seub
        else:
            range_small,range_big = int(seub[seub[5:].index("0")+5:]),int(next_seub[next_seub[5:].index("0")+5:])
            direction = 0 #seub < next_seub
        for seq_num in range(range_small+10,range_big,10):
            range_scaf = seub[seub[5:].index("0")+5:]
            insert_seub = seub.replace(str(range_scaf),"0"*(len(range_scaf)-len(str(seq_num)))+str(seq_num))
            print range_small,range_big
            print insert_seub,
            if insert_seub in nohomo_seub:
                print "^_^"
        #if - > 10:
        #    print range_small,range_big,range_big-range_small
        #    for insert_seub in nohomo_seub:
        #        if seub2record[insert_seub].seqid == scaf:
        #            print insert_seub,int(insert_seub[insert_seub[5:].index("0")+5:])
#for seub in nohomo_seub:
#    print seub,seub.split("0")[1],int(seub[seub[5:].index("0")+5:]),seub2record[seub].seqid,seub2record[seub].start,seub2record[seub].end,seub2record[seub].strand
'''





'''
for scer in scer2seub:
    if len(scer2seub[scer]) > 0 and scer not in scer121seub:
        i = j = scer_genes.index(scer)
        up_scer = scer_genes[i-1]
        down_scer = scer_genes[j+1]
        while up_scer not in scer121seub:
            try:
                i -= 1
                up_scer = scer_genes[i-1]
            except:
                up_scer = "00"
                break
        while down_scer not in scer121seub:
            try:
                j += 1
                down_scer = scer_genes[j+1]
            except:
                down_scer = "00"
                break
        up_seub_seqid = scer121seub[up_scer].split("0")[1]
        down_seub_seqid = scer121seub[down_scer].split("0")[1]
        flag = 0
        for gene in scer2seub[scer]:
            if gene.split("0")[1] == up_seub_seqid or gene.split("0")[1] == down_seub_seqid:
                seub = gene
                flag += 1
        if flag == 1:
            scer121seub[scer] = seub
            seub121scer[seub] = scer
            print "up:", up_scer, scer121seub[up_scer]
            print "down:", down_scer, scer121seub[down_scer]
            print "Set: ", scer, seub

        #print "UP:",scer121seub[up_scer],scer121seub[up_scer].split("0")[1]
        #print "Down:",scer121seub[down_scer],scer121seub[down_scer].split("0")[1]
        #print "?:",scer2seub[scer]


for seub in seub2scer:
    if len(seub2scer[seub]) > 0 and seub not in seub121scer:
        i = j = seub_genes.index(seub)
        up_seub = seub_genes[i-1]
        down_seub = seub_genes[j+1]
        while up_seub not in seub121scer:
            try:
                i -= 1
                up_seub = seub_genes[i-1]
            except:
                up_seub = "00"
                break
        while down_seub not in seub121scer:
            try:
                j += 1
                down_seub = seub_genes[j+1]
            except:
                down_seub = "00"
                break
        up_scer_seqid = seub121scer[up_seub][1]
        down_scer_seqid = seub121scer[down_seub][1]
        flag = 0
        for gene in seub2scer[seub]:
            if gene[1] == up_scer_seqid or gene[1] == down_scer_seqid:
                scer = gene
                flag += 1
        if flag == 1:
            scer121seub[scer] = seub
            seub121scer[seub] = scer
            print "up:", up_seub, seub121scer[up_seub]
            print "down:", down_seub, seub121scer[down_seub]
            print "Set: ", seub, scer 

        #print "UP:",seub121scer[up_seub],seub121scer[up_seub][1]
        #print "Down:",seub121scer[down_seub],seub121scer[down_seub][1]
        #print "?:",seub2scer[seub]

'''

        


'''
content = []
#freq_stack = Stack(20)
for scer in scer_genes:
    seub_ids = scer2seub[scer.id] if scer.id in scer2seub else ""
    line = "%s\t"*9+"%s\n"
    if seub_ids:
        if len(seub_ids) == 1:
            seub = seub2record[seub_ids[0]]
            content.append(line%(scer.id,scer.seqid,scer.start,scer.end,scer.strand,\
                    seub.id,seub.seqid,seub.start,seub.end,seub.strand))
            #freq_stack.push(seub.seqid)
            flag = 'only'
        else:
            for seub_id in seub_ids:
                seub = seub2record[seub_id]
                #if seub_scaf2scer_chr[seub.seqid] == scer.seqid and 
                if seub.seqid == main_seqid:
                    if not (flag == 'no' or flag == 'only'):
                        conflicts.append(seub.seqid)
                    content.append(line%(scer.id,scer.seqid,scer.start,scer.end,scer.strand,\
                            seub.id,seub.seqid,seub.start,seub.end,seub.strand))
                    flag = 'yes'
                else:
                    content.append(("---"+line)%(scer.id,scer.seqid,scer.start,scer.end,scer.strand,\
                            seub.id,seub.seqid,seub.start,seub.end,seub.strand))
                    flag = 'no'
    else:
        #freq_stack.push("")
        content.append(line%(scer.id,scer.seqid,scer.start,scer.end,scer.strand,\
                "","","","",""))

print set(conflicts)
#for i,line in enumerate(content):
#    if not line.startswith("---"):
#        main_seqid = line.split("\t")[6]
#    else:
#        seub_seqid = line.split("\t")[6]

f = open(out_put_file,"w")
for line in content:
    f.write(line)
'''
