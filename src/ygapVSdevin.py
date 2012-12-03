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
    os.system("makeblastdb -in %s -dbtype nucl -out %s"%(ygap_nc_fsa,ygap_nc_db))
    os.system("makeblastdb -in %s -dbtype nucl -out %s"%(devin_nc_fsa,devin_nc_db))

    os.system("blastn -evalue 0.00001 -max_target_seqs 1 -strand plus "+\
        "-max_hsps_per_subject 1 " +\
        "-db %s -query %s -out %s "%(devin_nc_db,ygap_nc_fsa,ygapVSdevin) +\
        "-outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart " +\
        "qend sstart send \"")

    os.system("blastn -evalue 0.00001 -max_target_seqs 1 -strand plus "+\
        "-max_hsps_per_subject 1 " +\
        "-db %s -query %s -out %s "%(ygap_nc_db,devin_nc_fsa,devinVSygap) +\
        "-outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart " +\
        "qend sstart send \"")

def write_fsa(ygap,devin):
    #wirte fsa file of their anontation 
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

def get_seq(fsafile):
    #read chromosome sequsnece file 
    seq_dict = {}
    for record in SeqIO.parse(fsafile,"fasta"):
        record.id = int(record.id.split("_")[-1])
        seq_dict[record.id] = record.seq
    return seq_dict

def init():
    #generate blast compare file 
    scaf2seqid = get_scaf2seqid()

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

    write_fsa(ygap,devin)
    run_blast()
    return 

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
    for devin_gene in new_match:
        pairs[devin_gene] = {}
        if len(new_match[devin_gene]) == 1:
            ygap_gene = new_match[devin_gene].keys()[0]
            relation = new_match[devin_gene][ygap_gene]
            if len(new_match[ygap_gene]) == 1:
                if devin_gene not in new_match[ygap_gene]:
                    print "Is it possible? A 1:1 B but B 1:1 C"
                elif new_match[ygap_gene][devin_gene] == relation: 
                    pairs[devin_gene][ygap_gene] = relation
                else:
                    print "conflicts about relationship at %s to %s"%(devin_gene,ygap_gene)
            else:
                if devin_gene not in new_match[ygap_gene]:
                    print "Strange! A 1:1 B but B 1:n -A "
                elif new_match[ygap_gene][devin_gene] != relation:
                    print "conflicts about relationship at %s to %s"%(devin_gene,ygap_gene)
                else:
                    pairs[devin_gene][ygap_gene] = relation
        else:
            for ygap_gene in new_match[devin_gene]:
                relation = new_match[devin_gene][ygap_gene]
                if len(new_match[ygap_gene]) == 1:
                    if devin_gene not in new_match[ygap_gene]:
                        print "Is it possible? A n:1 B but B 1:1 C"
                    elif new_match[ygap_gene][devin_gene] == relation:
                        pairs[devin_gene][ygap_gene] = relation
                    else:
                        print "conflicts about relationship at %s to %s"%\
                                (devin_gene,ygap_gene)
                else:
                    if devin_gene not in new_match[ygap_gene]:
                        print "Strange! A n:1 B but B 1:n -A "
                    elif new_match[ygap_gene][devin_gene] != relation:
                        print "conflicts about relationship at %s to %s"%\
                                (devin_gene,ygap_gene)
                    else:
                        pairs[devin_gene][ygap_gene] = relation
    return pairs

def write_summary(handle,pairs,devin_none,ygap_none):
#write_summary file and up 1000 stream sequences
    def write_up_down_fsa(handle, id_list, id_dict):
        for id in id_list:
            handle.write(">%s\tup:%d\tdown:%d\n"%(id,len(id_dict[id].up_1000),len(id_dict[id].down_1000)))
            handle.write("%s\n"%str(id_dict[id].up_1000))
            handle.write("%s\n"%str(id_dict[id].orf))
            handle.write("%s\n"%str(id_dict[id].down_1000))
            
    def write_want_type(handle, want_type):
        lines = []
        devin = []
        ygap = []
        for gene_1 in pairs:
            if gene_1 in ygap_gene_dict:
                for gene_2 in pairs[gene_1]:
                    if pairs[gene_1][gene_2] == want_type:
                        ygap.append(gene_1)
                        devin.append(gene_2)
                        if len(ygap_gene_dict[gene_1].orf) == len(devin_gene_dict[gene_2].orf):
                            vslen = "="
                        elif len(ygap_gene_dict[gene_1].orf) < len(devin_gene_dict[gene_2].orf):
                            vslen = "<"
                        else:
                            vslen = ">"
                        lines.append("%s\t%s\t%s\t%s\n"%(gene_1,vslen,gene_2,pairs[gene_1][gene_2]))
        handle.write(">>>%s: %d pairs(devin %d, ygap %d)\n"%\
                (want_type,len(lines),len(set(devin)),len(set(ygap))))
        for line in lines:
            handle.write(line)
        f = open("%s%s.%s.fsa"%(Out_Path,SP,want_type),"w")
        write_up_down_fsa(f,set(devin),devin_gene_dict)
        write_up_down_fsa(f,set(ygap),ygap_gene_dict)

    write_want_type(handle,"exact match")
    write_want_type(handle,"disagree at stop")
    write_want_type(handle,"disagree at start")
    write_want_type(handle,"include")
    write_want_type(handle,"overlap")
    handle.write(">>>none:\t(devin:%d)\t(ygap:%d)\n"%(len(devin_none),len(ygap_none)))
    handle.write(">>>>ygap:\n")
    for line in ygap_none:
        handle.write("%s\n"%line)
    handle.write(">>>>devin:\n")
    for line in devin_none:
        handle.write("%s\n"%line)
    f = open("%s%s.%s.fsa"%(Out_Path,SP,"Disagree"),"w")
    write_up_down_fsa(f,set(devin_none),devin_gene_dict)
    write_up_down_fsa(f,set(ygap_none),ygap_gene_dict)

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
            seqid = int(seqid)
            orf_start,orf_end = coord[coord.find("(")+1:coord.find(")")].split("..")
            orf_start = int(orf_start)
            orf_end = int(orf_end)
            start = max(int(orf_start) - 1000,1)
            end = min(int(orf_end) + 1000,len(devin_seq[int(seqid)]))
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

def 
#============================================================================#
SP = "scer"
ygap_anno = open("/Users/bingwang/zen/yeast_anno_pipe/YGAP_predict/%s.ygap.txt"%SP)
ygap_fsa = "/Users/bingwang/zen/yeast_anno_pipe/YGAP_predict/%s.genome.txt"%SP
ygap_corr = "/Users/bingwang/zen/yeast_anno_pipe/YGAP_predict/%s.ygap.corr"%SP
ygap_nc_fsa = "%s%s.ygap.fsa"%(Out_Path,SP)
ygap_nc_db = "%s%s.ygap.blastdb"%(Out_Path,SP)
ygapVSdevin = "%s%s.ygapVSdevin.xls"%(Out_Path,SP)
devin_gff = open("/Users/bingwang/zen/yeast_anno_pipe/SSS_dataset/%s.gff"%SP)
devin_fsa = "/Users/bingwang/zen/yeast_anno_pipe/SSS_dataset/%s.mergedscaf.fsa"%SP
devin_nc_fsa = "%s%s.devin.fsa"%(Out_Path,SP)
devin_nc_db = "%s%s.devin.blastdb"%(Out_Path,SP)
devinVSygap = "%s%s.devinVSygap.xls"%(Out_Path,SP)
summary = "%s%s.summary"%(Out_Path,SP)
ygap_seq = get_seq(ygap_fsa)
devin_seq = get_seq(devin_fsa)
#===============================================================================#
#init()
ygap_gene_dict = get_gene_dict(ygap_nc_fsa,ygap_seq)
devin_gene_dict = get_gene_dict(devin_nc_fsa,devin_seq)

devin_as_q = blast_result_parse(devinVSygap)
ygap_as_q = blast_result_parse(ygapVSdevin)

pairs = combine_two_matches(devin_as_q, ygap_as_q)

devin_none = [id for id in devin_gene_dict if id not in pairs]
ygap_none = [id for id in ygap_gene_dict if id not in pairs]

write_summary(open(summary,"w"),pairs,devin_none,ygap_none)

os.system("rm %s*.nhr"%Out_Path)
os.system("rm %s*.nin"%Out_Path)
os.system("rm %s*.nsq"%Out_Path)
'''
for id in ygap_none:
    seqid,coord = id.split("__")[1:3]
    seqid = int(seqid)
    orf_start,orf_end = coord[coord.find("(")+1:coord.find(")")].split("..")
    orf_start = int(orf_start)
    orf_end = int(orf_end)
    start = max(orf_start - 1000,1)
    end = min(orf_end + 1000,len(devin_seq[seqid]))
    if "complement" in coord:
        f.write(">ygap:%s\tup:%d\tdown:%d\n"%(id,end-orf_end,orf_start-start))
        f.write(str(devin_seq[seqid][orf_end:end].reverse_complement())+"\n")
        f.write(str(devin_seq[seqid][orf_start-1:orf_end].reverse_complement())+"\n")
        f.write(str(devin_seq[seqid][start-1:orf_start-1].reverse_complement())+"\n")
    else:
        f.write(">ygap:%s\tup:%d\tdown:%d\n"%(id,orf_start-start,end-orf_end))
        f.write(str(devin_seq[seqid][start-1:orf_start-1])+"\n")
        f.write(str(devin_seq[seqid][orf_start-1:orf_end])+"\n")
        f.write(str(devin_seq[seqid][orf_end:end])+"\n")
    seq_dict = {}

'''
'''
pairs = {}
pair_type = {}
for name in matched_gene:
    pairs[name] = Pair()
    pair_type[name] = {}

for name in matched_gene:
    if name in ygapid_list:
        pairs[name].devin |= set([devin_as_q[name].list+ygap_as_q[name].list])
        for devin_name in pairs[name].devin:
            pair_type[name][devin_name] = #TODO relationship
            pairs[name].ygap |= set([devin_as_q[devin_name].list+ygap_as_q[devin_name].list]) 
    else:
        pairs[name].ygap |= set([devin_as_q[name].list+ygap_as_q[name].list])
        for ygap_name in pairs[name].ygap:
            pair_type[name][ygap_name] = #TODO relationship
            pairs[name].devin |= set([devin_as_q[ygap_name].list+ygap_as_q[ygap_name].list])

devin_match = []
devin_db_match = []
ygap_match = []
ygap_db_match = []
devin_overlap = []
ygap_overlap = []

f = open(devinVSygap)
for line in f:
    devinid,ygapid,pident,mismatch,gapopen,\
    devinstart,devinend,ygapstart,ygapend = line.strip().split("\t")
    pident = float(pident)
    devinlength = devinid[devinid.find("(")+1:devinid.find(")")]
    devinlength = int(devinlength.split("..")[1]) - int(devinlength.split("..")[0]) + 1
    ygaplength = ygapid[ygapid.find("(")+1:ygapid.find(")")]
    ygaplength = int(ygaplength.split("..")[1]) - int(ygaplength.split("..")[0]) + 1
    length = ygaplength if ygaplength == devinlength else 1
    if pident == 100 :
        if devinstart == ygapstart == "1" and devinend == ygapend == str(length):
            devin_match.append(devinid)
            ygap_db_match.append(ygapid)
        else:
            devin_overlap.append(devinid)

f = open(ygapVSdevin)
for line in f:
    ygapid,devinid,pident,mismatch,gapopen,\
    ygapstart,ygapend,devinstart,devinend = line.strip().split("\t")
    pident = float(pident)
    devinlength = devinid[devinid.find("(")+1:devinid.find(")")]
    devinlength = int(devinlength.split("..")[1]) - int(devinlength.split("..")[0]) + 1
    ygaplength = ygapid[ygapid.find("(")+1:ygapid.find(")")]
    ygaplength = int(ygaplength.split("..")[1]) - int(ygaplength.split("..")[0]) + 1
    length = ygaplength if ygaplength == devinlength else 1
    if pident == 100 :
        if devinstart == ygapstart == "1" and devinend == ygapend == str(length):
            ygap_match.append(ygapid)
            devin_db_match.append(devinid)
        else:
            ygap_overlap.append(ygapid)

'''
