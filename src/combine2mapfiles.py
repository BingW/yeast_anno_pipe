# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"

suvaONscer = "/Users/bingwang/zen/yeast_anno_pipe/suvaONscer.txt"
seubONscer = "/Users/bingwang/zen/yeast_anno_pipe/seub3000ONscer.txt"
suva_seub_pair = "/Users/bingwang/zen/yeast_anno_pipe/suva_seub_pair.tab"

out_put_pillar = "/Users/bingwang/zen/yeast_anno_pipe/scer_suva_seub_pillar.tab"

scer_order = []
scer2block = {}
seub_order = []
seub2block = {}
scer121seub = {}
for line in open(seubONscer):
    line = line.replace("\n","").split("\t")
    scer,seub = line[0],line[5]
    if scer:
        scer_order.append(scer)
        scer2block[scer] = line[:5]
        scer121seub[scer] = seub
    if seub:
        seub_order.append(seub)
        seub2block[seub] = line[5:]

suva_order = []
scer121suva = {}
suva2block = {}
for line in open(suvaONscer):
    line = line.replace("\n","").split("\t")
    scer,suva = line[0],line[5]
    if scer:
        scer121suva[scer] = suva 
    if suva:
        suva_order.append(suva)
        suva2block[suva] = line[5:]

suva121seub = {}
for line in open(suva_seub_pair):
    line = line.replace("\n","").split("\t")
    suva,seub = line[0],line[1]
    suva121seub[suva] = seub

#write lines that:
#Scer Suva Seub
#  +   +    +
#  +   +    -  
#  +   -    +
#  +   -    -
#  -   +    -  Need to be checked whether - + - or - + +

content = []
seub_mapped = []
for line in open(suvaONscer):
    line = line.replace("\n","").split("\t")
    scer,suva = line[0],line[5]

    scer_block = scer2block[scer] if scer else [""]*5
    suva_block = suva2block[suva] if suva else [""]*5

    if scer and scer121seub[scer]:
        seub = scer121seub[scer]
        seub_mapped.append(seub)
        seub_block = seub2block[seub]
    else:
        seub_block = [""]*5
    content.append("\t".join(scer_block+suva_block+seub_block))

#check lines - + - that seub may has homolog with suva
for i,line in enumerate(content):
    line = line.split("\t")
    scer,suva,seub = line[0],line[5],line[10]
    if suva and not scer and suva in suva121seub:
        seub = suva121seub[suva]
        if seub in seub2block and seub not in seub_mapped:
            seub_block = seub2block[seub]
            seub_scaf = seub_block[1]
            j = i - 1
            while j>0 and not content[j].split("\t")[10]:
                j -= 1
            seub_up_scaf = content[j].split("\t")[11] if j > 0 else ""
            k = i + 1
            while k<len(content) and not content[k].split("\t")[10]:
                k += 1
            seub_down_scaf = content[k].split("\t")[11] if k < len(content) else ""
            if seub_scaf == seub_up_scaf or seub_scaf == seub_down_scaf:
                seub_mapped.append(seub)
                content[i] = "\t".join(line[:10] + seub_block)
            else:
                print suva,seub

#insert - - + based on seub's order
seub_order = iter(seub_order)
now_seub = seub_order.next()
new_content = []
for line in content:
    new_content.append(line)
    line = line.split("\t")
    scer,seub = line[0],line[10]
    if scer and seub:
        while seub != now_seub:
            seub_block = seub2block[now_seub]
            if now_seub not in seub_mapped:
                new_content.insert(-1,"\t".join([""]*10+seub_block))
            now_seub = seub_order.next()
        now_seub = seub_order.next()

open(out_put_pillar,"w").write("\n".join(new_content))


#TODO conflicts like continues - + +, + - + with same seub.
'''
#combine lines like:
#- + + and + - + where two seub genes are the same one
for i,line in enumerate(content):
    line = line.split("\t")
    scer,suva,seub = line[0],line[5],line[10]
    if seub and seub == content[i+1].split("\t")[10]:
        if scer and suva:
            print "some thing wrong @ line %d"%i
            print "\t".join(line)
            print content[i+1]
            continue
        elif scer:
            new_line = "\t".join(line[:5] + content[i+1].split("\t")[5:10] + line[10:])
        else:
            new_line = "\t".join(content[i+1].split("\t")[:5]+line[5:])
        print "\t".join(line)
        print content[i+1],"->"
        print new_line
'''

"""
content = []
seub_order = iter(seub_order)
now_seub = seub_order.next()
for line in open(suvaONscer):
    line = line.replace("\n","").split("\t")
    scer,suva = line[0],line[5]
    scer_block = scer2block[scer] if scer else [""]*5
    suva_block = suva2block[suva] if suva else [""]*5

    if scer and scer121seub[scer]:
        seub = scer121seub[scer]
        while seub != now_seub:
            seub_block = seub2block[now_seub]
            content.append("\t".join([""]*10+seub_block))
            now_seub = seub_order.next()
        seub_block = seub2block[seub]
        now_seub = seub_order.next()
    elif suva and suva in suva121seub:
        seub = suva121seub[suva]
        if seub in seub2block:
            seub_block = seub2block[seub]
        else:
            seub_block = [""]*5
    else:
        seub_block = [""]*5

    content.append("\t".join(scer_block+suva_block+seub_block))

open(out_put_pillar,"w").write("\n".join(content))
"""
