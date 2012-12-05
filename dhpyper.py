#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"

import os
import rpy2.robjects as robjects

work_path = "/Users/bingwang/zen/get_up_out/"
weight_ms_folder = work_path + "wg/"
sample_folder = work_path + "Sample/"
write_file = work_path + "hyper_test.txt"
meme_matrix_folder = work_path + "MEME_matrix/"

class Record():
    def __init__(self):
        pass

def mx_weight_parse(line):
    r = Record()
    r.name,r.strand,r.start,r.end,r.sequence,r.weight,r.pval \
            = line.strip().split("\t")
    r.weight = float(r.weight)
    r.pval = float(r.pval)
    return r

def r_dhyper(x,m,n,k):
    dhyper = robjects.r['dhyper']
    p = dhyper(x,m,n,k).r_repr()
    return p

sample_dict = {}
for sample in os.listdir(sample_folder):
    if sample.endswith("moreinfo.txt"):
        sp_name = sample[6:10]
        sample_dict[sp_name] = {}
        f = open(sample_folder+sample)
        f.next()
        for line in f:
            r = mx_weight_parse(line)
            sample_dict[sp_name][r.name] = r
        

bg_dict = {}
for bg in os.listdir(weight_ms_folder):
    if bg.endswith("_wg.txt"):
        sp_name = bg[14:18]
        bg_dict[sp_name] = {}
        f = open(weight_ms_folder+bg)
        f.next()
        for line in f:
            r = mx_weight_parse(line)
            bg_dict[sp_name][r.name] = r

f = open(write_file,"w")
f.write("%s\t%s\t%s\t%s\t%s\t%s\n"%\
        ("SP","P",">5.2_sample",">5.2_genome","total_sample","total_genome"))
for sp in sample_dict:
    total_genome = len(bg_dict[sp])
    total_sample = len(sample_dict[sp])
    sign_genome = len([r for r in bg_dict[sp] if bg_dict[sp][r].weight > 5.2])
    sign_sample = len([r for r in sample_dict[sp] if sample_dict[sp][r].weight > 5.2])
    n = total_genome - sign_genome
    p = r_dhyper(sign_sample,sign_genome,n,total_sample)
    f.write("%s\t%s\t%d\t%d\t%d\t%d\n"%\
            (sp,p,sign_sample,sign_genome,total_sample,total_genome))
