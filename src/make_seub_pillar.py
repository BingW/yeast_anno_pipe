#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"

#1. get scer 1:1 seub
#2. get scer 1:1 suva 1:1 anc
#3. map scer:suva:seub use 1:1:1
#4. insert suva specific genes, check seub gene use anc and blast
#5. insert seub specific genes
#6. output
seubONscer_file = "/Users/bingwang/zen/yeast_anno_pipe/seubONscer.txt"
suva_tab_file = "/Users/bingwang/zen/yeast_anno_pipe/YGOB_dataset/Suvarum_genome.tab"
pillar_file = "/Users/bingwang/zen/yeast_anno_pipe/YGOB_dataset/Pillars.tab"
write_file = "/Users/bingwang/zen/yeast_anno_pipe/Pillar_add_seub.txt"

content = []
scer2seub = {}
scer_seq = []
for line in open(seubONscer_file):
    record = line.split("\t")
    scer2seub[record[0]] = record[5]
    scer_seq.append(record[0])

'''
#2. get scer 1:1 suva
scer2suva = {}
suva2scer = {}
scer2suva[""] = ""
for i,line in enumerate(open(suva_tab_file)):
    content = line.split("\t")
    suva = content[0]
    if not "trna" in suva and not "CEN" in suva:
        scer = content[8][:content[8].find("(")].strip()
        homo_type = content[8][content[8].find("(")+1:content[8].find(")")]
        if len(scer) not in [0,7,9,12]:
            print len(scer),i,content[8]
        if homo_type not in ["REAL","HSP","PSEUDO","REPEAT","NOVEL","HYPO","INTER"]:
            print homo_type
        suva2scer[suva] = scer
        if scer != "":
            if scer in scer2suva:
                print "conflict:",scer2suva[scer],suva,"share same homo"
            else:
                scer2suva[scer] = suva
'''

content = []
for line in open(pillar_file):
    if 'Scer_' in line or "trna" in line:
        continue
    genes = line.split("\t")
    scer_1,scer_2,anc,suva_1,suva_2 = genes[11],genes[-12],genes[12],genes[8],genes[-9]
    dubious = ["YIL170W","YIR044C","YIL171W","YOL153C","YPL276W","YPL275W","YAR061W","YAR073W","YAR075W"]
    if scer_1 in dubious or scer_2 in dubious:
        continue
    seub_1 = "---" if scer_1 == "---" else scer2seub[scer_1]
    seub_2 = "---" if scer_2 == "---" else scer2seub[scer_2]
    if [scer_1,suva_1,seub_1,anc,scer_2,suva_2,seub_2].count('---') == 7:
        continue
    content.append("\t".join([scer_1,suva_1,seub_1,anc,scer_2,suva_2,seub_2]))

open(write_file,"w").write("\n".join(content))
