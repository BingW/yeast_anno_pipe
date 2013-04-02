#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"

map_file = "/Users/bingwang/zen/yeast_anno_pipe/seubONscer.txt"
content = []
for line in open(map_file):
    content.append(line.split("\t"))

headntail = []
scersp = []
gap = []
for i,line in enumerate(content):
    if line[6] == "":
        up_i = i - 1
        while up_i >= 0 and content[up_i][6] == "":
            up_i -= 1
        up_line = [] if up_i < 0 else content[up_i] 
        down_i = i + 1
        while down_i < len(content) and content[down_i][6] == "":
            down_i += 1
        down_line = [] if down_i == len(content) else content[down_i] 
        if up_line == [] or down_line == [] or not (line[1] == up_line[1] == down_line[1]):
            headntail.append(line)
        elif up_line[1] == down_line[1] == line[1]:
            if up_line[6] == down_line[6]:
                scersp.append(line)
            else:
                gap.append(line)

def print_region_list(region_list,region_type):
    count = 0
    for line in region_list:
        i = content.index(line)
        if content[i-1][5] != "":
            up_num = int(content[i-1][2])
            up_scaf = content[i-1][6]
            print "%s @ %s\t%s.."%(region_type,content[i-1][1],content[i-1][3]),
        if content[i+1][5] != "":
            down_num = int(content[i+1][3])
            down_scaf = content[i+1][6]
            print content[i+1][3],
            print "length:%d\t"%(down_num - up_num),
            print "up_scaf:%s\tdown_scaf:%s\t"%(up_scaf,down_scaf)
            count += 1
    print count

print_region_list(gap,"gap")
#print "coverage:", 1-((len(headntail)+len(gap)) * 1.0 / len(content))
#for gene in scersp:
#    print "\t".join(gene),
#print "total:",len(scersp)
