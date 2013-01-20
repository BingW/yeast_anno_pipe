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

def YGAPIterator(handle):
    #for record in YGAPIterator(handle):
    #    record.name = gene name
    #    record.strand = strand
    #    record.start = lowest coordinate
    #    record.end = highest coordinate
    #    record.ygob = YGOB browser status (all ON)
    #    record.seqid = chromosome/scaffold num
    #    record.shortname = gene short name
    #    record.coord = gene coordinates (including introns)
    #    record.orth = Ancestral Scer1|Scer2
    #    record.type = gene type(empty,PROTEIN,TKP/TY)
    #    record.pillar = YGOB pillar num
    #    record.tag = annotation tag(Cor/Fra/Man/NNN)
    #        Cor: successful frameshift correction
    #        Fra: should have frameshift but cannot correct
    #        Man: need to be manual attention(may untranslatable, overlap, strange length)
    #        NNN: N at begins or ends of the coding region
    #    record.anno = Annotation route(simple/multi/GETORF/DOGS)
    #        Simple: easily assigned to YGOB pillar
    #        Multi: gene has multigen family and its pillar may be unreliable
    #        GETORF: generate by GETORF
    #        DOGS: generate by DOGS
    class record():
        def __init__(self):
            pass
    while 1:
        line = handle.readline()
        if line.count("..") > 1:
            print "intron @:%s"%line
        if line.startswith("#"): continue
        if not line: return
        r = record()
        (r.name,r.strand,r.start,r.end,r.ygob,r.seqid,r.shortname,\
            r.coord,r.orth,r.type,r.pillar,r.tag,r.anno) = line.split("\t")
        r.anno = r.anno.strip()
        r.strand = "+" if r.strand == "1" else "-"
        #r.orth = [r.orth.split(" ")[0]] + r.orth.split(" ")[1].split("|")
        yield r

def run_blast():
    #first make 2 db
    #use blastn do the blast(!same strand)
    #needs 2 files (ygap_nc_fsa,devin_nc_fsa)
    #generate 4 files (ygap_nc_db,devin_nc_db,ygapVSdevin,devinVSygap)
    def make_db(fsa,db):
        os.system("makeblastdb -in %s -dbtype nucl -out %s"%(fsa,db))
    def run_blast(db,fsa,out_file):
        os.system("blastn -evalue 0.00001 -max_target_seqs 1 -strand plus "+\
            "-max_hsps_per_subject 1 " +\
            "-db %s -query %s -out %s "%(db,fsa,out_file) +\
            "-outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart " +\
            "qend sstart send \"")

    make_db(ygap_nc_fsa,ygap_nc_db)
    make_db(devin_nc_fsa,devin_nc_db)
    make_db(ygap_fsa,ygap_scaf_db)
    run_blast(devin_nc_db,ygap_nc_fsa,ygapVSdevin)
    run_blast(ygap_scaf_db,ygap_nc_fsa,ygapONscaf)
    run_blast(ygap_nc_db,devin_nc_fsa,devinVSygap)
    run_blast(ygap_scaf_db,devin_nc_fsa,devinONscaf)

    if REF:
        make_db(Ref_nc_fsa,Ref_nc_db)
        run_blast(Ref_nc_db,ygap_nc_fsa,ygapVSRef)  
        run_blast(ygap_nc_db,Ref_nc_fsa,RefVSygap)  
        run_blast(Ref_nc_db,devin_nc_fsa,devinVSRef)
        run_blast(devin_nc_db,Ref_nc_fsa,RefVSdevin)

def write_fsa(ygap,devin,refs = None):
    #write fsa file of their anontation 
    #input file is a list of records
    #generate 2 files:(ygap_nc_fsa, devin_nc_fsa)
    f = open(ygap_nc_fsa,"w")
    for r in ygap:
        line = ">%s__%s__%s__%s__%s__%s\n"%(r.name,r.seqid,r.coord,r.orth,r.tag,r.anno)
        f.write(line.replace(" ","_"))
        f.write("%s\n"%str(r.seq))

    f = open(devin_nc_fsa,"w")
    for r in devin:
        line = ">%s__%s__%s__%s\n"%("|".join(r.attributes["Gene"]),r.seqid,r.coord,r.orth)
        f.write(line.replace(" ","_"))
        f.write("%s\n"%str(r.seq))

    if REF:
        f = open(Ref_nc_fsa,"w")
        for r in refs:
            line = ">%s__%s__%s\n"%(r.id,r.seqid,r.coord)
            f.write(line)
            f.write("%s\n"%(str(r.seq)))

def get_seq(fsafile):
    #read chromosome sequsnece file 
    seq_dict = {}
    for record in SeqIO.parse(fsafile,"fasta"):
        try:
            record.id = int(record.id.split("_")[-1])
        except:
            pass
        seq_dict[record.id] = record.seq
    return seq_dict

def init():
    def get_scaf2seqid():
        #return secf2seqid dict
        #needs ygap_corr file
        scaf2seqid = {}
        f = open(ygap_corr)
        for line in f:
            scaf,seqid,chrname = line.split("\t")
            scaf = scaf.split("_")[-1]
            scaf2seqid[scaf] = seqid 
        return scaf2seqid
    #generate blast compare file 
    scaf2seqid = get_scaf2seqid()

    if REF:
        refs = []
        dubious = []
        for record in gff_parse.gffIterator(Ref_anno):
            if record.type == "gene" and record.attributes["orf_classification"][0] == "Dubious":
                dubious.append(record.attributes["ID"][0])
            elif record.type == "gene" and not record.attributes["ID"][0].startswith("Q"):
                seqid = record.seqid
                start = int(record.start)
                end = int(record.end)
                record.id = record.attributes["ID"][0]
                if record.strand == "+":
                    record.coord = "(%d..%d)"%(start,end)
                    record.seq = Ref_seq[seqid][start-1:end]
                elif record.strand == "-":
                    record.coord = "complement(%d..%d)"%(start,end)
                    record.seq = Ref_seq[seqid][start-1:end].reverse_complement()
                else:
                    print "Check strand:%s"%record.type
                refs.append(record)

    ygap = []
    for record in YGAPIterator(ygap_anno):
        if record.type == "PROTEIN" or record.type == "":
            seqid = int(record.seqid)
            start = int(record.start)
            end = int(record.end)
            if record.strand == "+":
                record.seq = ygap_seq[seqid][start-1:end]
            elif record.strand == "-":
                record.seq = ygap_seq[seqid][start-1:end].reverse_complement()
            else:
                print "Check strand:%s"%record.type
            ygap.append(record)

    devin = []
    for record in gff_parse.gffIterator(devin_gff):
        if record.type == "CDS" or record.type == "ORF":
            if REF:
                #exclude dubious genes
                flag = False
                for homo in record.attributes["SGD"]:
                    if homo in dubious:
                        flag = True
                if flag:
                    print "find %s is %s "%(record.attributes["Gene"],homo)
                    continue
            seqid = int(record.seqid)
            record.seqid = scaf2seqid[record.seqid]
            start = int(record.start)
            end = int(record.end)
            record.orth = "%s_%s"%(record.score,"|".join(record.attributes["SGD"]))
            if record.strand == "+":
                record.coord = "(%d..%d)"%(start,end)
                record.seq = devin_seq[seqid][start-1:end]
            elif record.strand == "-":
                record.coord = "complement(%d..%d)"%(start,end)
                record.seq = devin_seq[seqid][start-1:end].reverse_complement()
            else:
                print "Check strand:%s"%record.type
            devin.append(record)

    if REF:
        write_fsa(ygap,devin,refs)
    else:
        write_fsa(ygap,devin)

    run_blast()
    return 0

def blast_result_parse(queryVSdb):
    #q = query
    #s = seqdb
    matches = {}
    f = open(queryVSdb)
    for line in f:
        qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send = \
                line.strip().split("\t")
        pident = float(pident)
        if pident == 100:
            qgeneid = qseqid.split("__")[0]
            sgeneid = sseqid.split("__")[0]
            qcoord = qseqid[qseqid.find("(")+1:qseqid.find(")")]
            origin_qlength = int(qcoord.split("..")[1]) - int(qcoord.split("..")[0]) + 1
            scoord = sseqid[sseqid.find("(")+1:sseqid.find(")")]
            origin_slength = int(scoord.split("..")[1]) - int(scoord.split("..")[0]) + 1
            origin_length = origin_qlength if origin_qlength == origin_slength else 0
            matches[qgeneid] = {}
            if sgeneid not in matches:
                matches[sgeneid] = {}
            match_length = int(length)
            qstart = int(qstart)
            qend = int(qend)
            sstart = int(sstart)
            send = int(send)
            if match_length == origin_length:
                matches[qgeneid][sgeneid] = "exact match"
                matches[sgeneid][qgeneid] = "exact match"
            elif match_length == min(origin_qlength,origin_slength):
                if qstart == sstart == 1:
                    matches[qgeneid][sgeneid] = "disagree at stop"
                    matches[sgeneid][qgeneid] = "disagree at stop"
                elif qend == origin_qlength and send == origin_slength:
                    matches[qgeneid][sgeneid] = "disagree at start"
                    matches[sgeneid][qgeneid] = "disagree at start"
                elif min(qstart,sstart) == 1 and min(qend,send) == match_length:
                    matches[qgeneid][sgeneid] = "include"
                    matches[sgeneid][qgeneid] = "include"
                else:
                    print "Are there any other possiblity?Check %s and %s"%(qgeneid,sgeneid)
            elif match_length < min(origin_qlength,origin_slength):
                matches[qgeneid][sgeneid] = "overlap"
                matches[sgeneid][qgeneid] = "overlap"
            else:
                print "Are there any other possiblity?Check %s and %s"%(qgeneid,sgeneid)
    return matches

def combine_two_matches(matches_1,matches_2):
    #pairs[gene_1],[gene_2] = relation
    #    pairs are more restrict
    #    it checked competebality
    new_match = {}
    for name in matches_1:
        new_match[name] = {}
        for gene in matches_1[name]:
            new_match[name][gene] = matches_1[name][gene]
    for name in matches_2:
        if name in new_match:
            for gene in matches_2[name]:
                if gene not in new_match[name]:
                    new_match[name][gene] = matches_2[name][gene]
                elif new_match[name][gene] == matches_1[name][gene]:
                    pass
                else:
                    print "two matches have conflicts at %s to %s:"%(name,gene),
                    print "matches_1:%s\t\tmatches_2:%s"%(matches_1[name][gene],\
                            matches_2[name][gene])
        else:
            new_match[name] = {}
            for gene in matches_2[name]:
                new_match[name][gene] = matches_2[name][gene]
    pairs = {}
    for a_gene in new_match:
        pairs[a_gene] = {}
        if len(new_match[a_gene]) == 1:
            b_gene = new_match[a_gene].keys()[0]
            relation = new_match[a_gene][b_gene]
            if len(new_match[b_gene]) == 1:
                if a_gene not in new_match[b_gene]:
                    print "Is it possible? A 1:1 B but B 1:1 C"
                elif new_match[b_gene][a_gene] == relation: 
                    pairs[a_gene][b_gene] = relation
                else:
                    print "conflicts about relationship at %s to %s"%(a_gene,b_gene)
            else:
                if a_gene not in new_match[b_gene]:
                    print "Strange! A 1:1 B but B 1:n -A "
                elif new_match[b_gene][a_gene] != relation:
                    print "conflicts about relationship at %s to %s"%(a_gene,b_gene)
                else:
                    pairs[a_gene][b_gene] = relation
        else:
            for b_gene in new_match[a_gene]:
                relation = new_match[a_gene][b_gene]
                if len(new_match[b_gene]) == 1:
                    if a_gene not in new_match[b_gene]:
                        print "Is it possible? A n:1 B but B 1:1 C"
                    elif new_match[b_gene][a_gene] == relation:
                        pairs[a_gene][b_gene] = relation
                    else:
                        print "conflicts about relationship at %s to %s"%\
                                (a_gene,b_gene)
                else:
                    if a_gene not in new_match[b_gene]:
                        print "Strange! A n:1 B but B 1:n -A "
                    elif new_match[b_gene][a_gene] != relation:
                        print "conflicts about relationship at %s to %s"%\
                                (a_gene,b_gene)
                    else:
                        pairs[a_gene][b_gene] = relation
    return pairs

def write_summary(handle,pairs,a_dict,a_name,b_dict,b_name):

    print "%sVS%s:\t"%(a_name,b_name),
    #write_summary file and up 1000 stream sequences
    def write_up_down_fsa(handle, id_list, id_dict):
        for id in id_list:
            handle.write(">%s\tup:%d\tdown:%d\n"%(id,len(id_dict[id].up_1000),len(id_dict[id].down_1000)))
            handle.write("%s\n"%str(id_dict[id].up_1000))
            handle.write("%s\n"%str(id_dict[id].orf))
            handle.write("%s\n"%str(id_dict[id].down_1000))
    #def write_up_down_visual(handle, id_list, id_dict):
    def write_want_type(handle, want_type):
        lines = []
        a = []
        b = []
        record = []
        for gene_a in pairs:
            if gene_a in a_dict:
                for gene_b in pairs[gene_a]:
                    if gene_b in b_dict and pairs[gene_a][gene_b] == want_type:
                        a.append(gene_a)
                        b.append(gene_b)
                        record.append(gene_a)
                        record.append(gene_b)
                        if len(a_dict[gene_a].orf) == len(b_dict[gene_b].orf):
                            vslen = "="
                        elif len(a_dict[gene_a].orf) < len(b_dict[gene_b].orf):
                            vslen = "<"
                        else:
                            vslen = ">"
                        lines.append("%s\t%s\t%s\t%s\n"%(gene_a,vslen,gene_b,pairs[gene_a][gene_b]))
        handle.write(">>>%s: %d pairs(%s %d, %s %d)\n"%\
                (want_type,len(lines),a_name,len(set(a)),b_name,len(set(b))))
        print "%s:%d\t"%(want_type_short[want_type],len(set(b))),
        for line in lines:
            handle.write(line)
        return record
        #f = open("%s%s.%s.fsa"%(Out_Path,SP,want_type),"w")
        #write_up_down_fsa(f,set(a),a_dict)
        #write_up_down_fsa(f,set(b),b_dict)

    want_type_short = {"exact match":"E","disagree at stop":"De","disagree at start":"Ds",\
            "include":"I","overlap":"O"}
    record = []
    record += write_want_type(handle,"exact match")
    record += write_want_type(handle,"disagree at stop")
    record += write_want_type(handle,"disagree at start")
    record += write_want_type(handle,"include")
    record += write_want_type(handle,"overlap")
    a_none = [id for id in a_dict if id not in record]
    b_none = [id for id in b_dict if id not in record]
    handle.write(">>>none:\t(%s:%d)\t(%s:%d)\n"%(a_name,len(a_none),b_name,len(b_none)))
    print "F:%d\tM:%d\n"%(len(a_none),len(b_none))
    handle.write(">>>>%s:\n"%a_name)
    for line in a_none:
        handle.write("%s\n"%line)
    handle.write(">>>>%s:\n"%b_name)
    for line in b_none:
        handle.write("%s\n"%line)
    #f = open("%s%s.%s.fsa"%(Out_Path,SP,"Disagree"),"w")
    #write_up_down_fsa(f,set(a_none),a_dict)
    #write_up_down_fsa(f,set(b_none),b_dict)

def get_gene_dict(fsa_file,seq_dict):
    class gene():
        def __init__(self):
            pass
    gene_dict = {}
    id2descrip = {}
    f = open(fsa_file)
    for line in f:
        if line.startswith(">"):
            description = line[1:].strip()
            record = gene()
            record.id = line[1:].split("__")[0]
            seqid,coord = description.split("__")[1:3]
            try:
                seqid = int(seqid)
            except:
                pass
            orf_start,orf_end = coord[coord.find("(")+1:coord.find(")")].split("..")
            orf_start = int(orf_start)
            orf_end = int(orf_end)
            start = max(int(orf_start) - 1000,1)
            end = min(int(orf_end) + 1000,len(seq_dict[seqid]))
            record.start = orf_start
            record.end = orf_end
            record.seqid = seqid
            record.strand = "-" if "complement" in coord else "+"
            if "complement" in coord:
                record.up_1000 = seq_dict[seqid][orf_end:end].reverse_complement()
                record.orf = seq_dict[seqid][orf_start-1:orf_end].reverse_complement()
                record.down_1000 = seq_dict[seqid][start-1:orf_start-1].reverse_complement()
            else:
                record.up_1000 = seq_dict[seqid][start-1:orf_start-1]
                record.orf = seq_dict[seqid][orf_start-1:orf_end]
                record.down_1000 = seq_dict[seqid][orf_end:end]
            gene_dict[record.id] = record
    return gene_dict

def get_sub_dict(match_type,pairs,gene_dict):
    want_gene_dict = {}
    for gene_1 in pairs:
        if gene_1 in gene_dict:
            for gene_2 in pairs[gene_1]:
                if pairs[gene_1][gene_2] == match_type:
                    want_gene_dict[gene_1] = gene_dict[gene_1]
    return want_gene_dict


#============================================================================#
SP = "scer"
REF = True

#ref_in
Ref_anno = open("/Users/bingwang/zen/yeast_anno_pipe/SGD_data/scer.gff")
Ref_fsa = "/Users/bingwang/zen/yeast_anno_pipe/SGD_data/scer.fsa"
#ygap_in
ygap_anno = open("/Users/bingwang/zen/yeast_anno_pipe/YGAP_predict/%s.ygap.txt"%SP)
ygap_fsa = "/Users/bingwang/zen/yeast_anno_pipe/YGAP_predict/%s.genome.txt"%SP
ygap_corr = "/Users/bingwang/zen/yeast_anno_pipe/YGAP_predict/%s.ygap.corr"%SP
#devin_in
devin_gff = open("/Users/bingwang/zen/yeast_anno_pipe/SSS_dataset/%s.gff"%SP)
devin_fsa = "/Users/bingwang/zen/yeast_anno_pipe/SSS_dataset/%s.mergedscaf.fsa"%SP

#_ref_out
Ref_nc_fsa = "%s%s.ref.fsa"%(Out_Path,SP)
Ref_nc_db = "%s%s.ref.blastdb"%(Out_Path,SP)
RefVSygap = "%s%s.RefVSygap.xls"%(Out_Path,SP)
RefVSdevin = "%s%s.RefVSdevin.xls"%(Out_Path,SP)
dubious_file = "%s%s.dubious.tab"%(Out_Path,SP)
#ygap_out
ygap_nc_fsa = "%s%s.ygap.fsa"%(Out_Path,SP)
ygap_nc_db = "%s%s.ygap.blastdb"%(Out_Path,SP)
ygap_scaf_db = "%s%s.scaf.blastdb"%(Out_Path,SP)
ygapVSdevin = "%s%s.ygapVSdevin.xls"%(Out_Path,SP)
ygapVSRef = "%s%s.ygapVSRef.xls"%(Out_Path,SP)
ygapONscaf= "%s%s.ygapONscaf.xls"%(Out_Path,SP)
#devin_out
devin_nc_fsa = "%s%s.devin.fsa"%(Out_Path,SP)
devin_nc_db = "%s%s.devin.blastdb"%(Out_Path,SP)
devinVSygap = "%s%s.devinVSygap.xls"%(Out_Path,SP)
devinVSRef = "%s%s.devinVSRef.xls"%(Out_Path,SP)
devinONscaf= "%s%s.devinONscaf.xls"%(Out_Path,SP)
#summary_out
devin_VS_ygap_sum = "%s%s.devin_vs_ygap.summary"%(Out_Path,SP)
devin_VS_Ref_sum = "%s%s.devin_vs_ref.summary"%(Out_Path,SP)
ygap_VS_Ref_sum = "%s%s.ygap_vs_ref.summary"%(Out_Path,SP)
exact_VS_Ref_sum = "%s%s.exact_vs_ref.summary"%(Out_Path,SP)
dstart_VS_Ref_sum = "%s%s.dstart_vs_ref.summary"%(Out_Path,SP)
dstop_VS_Ref_sum = "%s%s.dstop_vs_ref.summary"%(Out_Path,SP)
include_VS_Ref_sum = "%s%s.include_vs_ref.summary"%(Out_Path,SP)
overlap_VS_Ref_sum = "%s%s.overlap_vs_ref.summary"%(Out_Path,SP)
ddevin_VS_Ref_sum = "%s%s.disdevin_vs_ref.summary"%(Out_Path,SP)
dygap_VS_Ref_sum = "%s%s.disygap_vs_ref.summary"%(Out_Path,SP)
all_VS_Ref_sum = "%s%s.united_vs_ref.summary"%(Out_Path,SP)
#sum_of_sum = "%s%s.sumofsum.summary"%(Out_Path,SP)

Ref_seq = get_seq(Ref_fsa)
ygap_seq = get_seq(ygap_fsa)
devin_seq = get_seq(devin_fsa)
##===============================================================================#
#init()
ygap_gene_dict = get_gene_dict(ygap_nc_fsa,ygap_seq)
devin_gene_dict = get_gene_dict(devin_nc_fsa,devin_seq)

devin_as_q = blast_result_parse(devinVSygap)
ygap_as_q = blast_result_parse(ygapVSdevin)
pairs_devin_ygap = combine_two_matches(devin_as_q, ygap_as_q)
write_summary(open(devin_VS_ygap_sum,"w"),pairs_devin_ygap,devin_gene_dict,"devin",ygap_gene_dict,"ygap")

if REF:
    Ref_gene_dict = get_gene_dict(Ref_nc_fsa,Ref_seq)

    devin_vs_ref = blast_result_parse(devinVSRef)
    ref_vs_devin = blast_result_parse(RefVSdevin)
    pairs_devin_ref = combine_two_matches(devin_vs_ref,ref_vs_devin)
    write_summary(open((devin_VS_Ref_sum),"w"),pairs_devin_ref,devin_gene_dict,"devin",Ref_gene_dict,"SGD")

    ygap_vs_ref = blast_result_parse(ygapVSRef)
    ref_vs_ygap = blast_result_parse(RefVSygap)
    pairs_ygap_ref = combine_two_matches(ygap_vs_ref,ref_vs_ygap)
    write_summary(open((ygap_VS_Ref_sum),"w"),pairs_ygap_ref,ygap_gene_dict,"ygap",Ref_gene_dict,"SGD")
    #write_want_type(handle,"exact match")
    #write_want_type(handle,"disagree at stop")
    #write_want_type(handle,"disagree at start")
    #write_want_type(handle,"include")
    #write_want_type(handle,"overlap")

    exact_match_dict = get_sub_dict("exact match",pairs_devin_ygap,ygap_gene_dict)
    write_summary(open((exact_VS_Ref_sum),"w"),pairs_ygap_ref,exact_match_dict,"exact match",Ref_gene_dict,"SGD")
    ddevin_dict = {gene: devin_gene_dict[gene] for gene in devin_gene_dict if gene not in pairs_devin_ygap}
    write_summary(open((ddevin_VS_Ref_sum),"w"),pairs_devin_ref,ddevin_dict,"ddevin",Ref_gene_dict,"SGD")
    dygap_dict = {gene: ygap_gene_dict[gene] for gene in ygap_gene_dict if gene not in pairs_devin_ygap}
    write_summary(open((dygap_VS_Ref_sum),"w"),pairs_ygap_ref,dygap_dict,"dygap",Ref_gene_dict,"SGD")
    all_devin_ygap_VS_Ref = {"exact match":[],"disagree at stop":[],"disagree at start":[],"include":[],"overlap":[]}
    set_1 = set([gene for gene in Ref_gene_dict if gene not in pairs_devin_ref])
    set_2 = set([gene for gene in Ref_gene_dict if gene not in pairs_ygap_ref])
    all_devin_ygap_VS_Ref["miss"] = set_1 & set_2
    for gene_1 in pairs_devin_ref:
        if gene_1 in Ref_gene_dict:
            for gene_2 in pairs_devin_ref[gene_1]:
                all_devin_ygap_VS_Ref[pairs_devin_ref[gene_1][gene_2]].append(gene_1)
    for gene_1 in pairs_ygap_ref:
        if gene_1 in Ref_gene_dict:
            for gene_2 in pairs_ygap_ref[gene_1]:
                all_devin_ygap_VS_Ref[pairs_ygap_ref[gene_1][gene_2]].append(gene_1)
    for item in all_devin_ygap_VS_Ref:
        print item,len(set(all_devin_ygap_VS_Ref[item]))
    print set([gene for gene in all_devin_ygap_VS_Ref["miss"] if not gene.startswith("Q")])


#=================================================================================#
#os.system("rm %s*.nhr"%Out_Path)
#os.system("rm %s*.nin"%Out_Path)
#os.system("rm %s*.nsq"%Out_Path)
