#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copyleft Bing Wang
# LICENCES: GPL v3.0

#vesrion: 0.1.1
#modification:
#check two *.fsa files should appear in the command
#check if the out_put.fsa file exists, report error and exit program.


__docformat__ = "epytext en"

import os
import sys

def print_help():
    print "USAGE\n"
    print "bing_RBBH [-h/-help] [-e e_value] [-f] [-p] fasta_file_1 fasta_file_2 out_put_file\n"
    print "OPTIONAL ARGUMENTS\n"
    print "-h/help"
    print "\t print help info"
    print "-e e_value"
    print "\t set blast evalue cut off"
    print "-f "
    print "\t out put full information include blast result"
    print "-p "
    print "\t the two fasta file is protein sequence rather than nucle tide sequence\n"

def check_file():
    fsa_files = [f for f in sys.argv if (f.endswith(".fsa") or f.endswith(".fasta"))]
    if len(fsa_files) != 2:
        print "You must enter exactly 2 fasta files(must end with .fsa or .fasta)"
        sys.exit()

    elif (sys.argv[-1].endswith(".fsa") or sys.argv[-1].endswith(".fasta")):
        print "out_put_file not find"
        print_help()
        sys.exit()

    out_file_path = "/".join(sys.argv[-1].split("/")[:-1])
    out_file_path = os.getcwd() if out_file_path == "" else out_file_path
    out_put_file = sys.argv[-1].split("/")[-1]
    if not os.path.exists(out_file_path):
        os.system("mkdir out_file_path")
    if out_put_file in os.listdir(out_file_path):
        print "The out_put file exists, pleast delete it or rename it first"
        sys.exit()

    try:
        assert (sys.argv[-2].endswith(".fsa") or sys.argv[-2].endswith(".fasta"))
        assert (sys.argv[-3].endswith(".fsa") or sys.argv[-2].endswith(".fasta"))
        #add asserts here
    except:
        print_help()
        sys.exit()

    

def reciprocal_best_blastn(nt_file_a,nt_file_b,e_val):
    #first make 2 db
    #use blastn do the blast
    #generate 2 blast results
    def make_db(nt_file):
        os.system("makeblastdb -in %s -dbtype nucl -out %s"%(nt_file, nt_file+".db"))
    def run_blast(db,nt_file,out_file):
        os.system("blastn -evalue %f -max_target_seqs 1 "%(e_val)+\
            "-db %s -query %s -out %s "%(db,nt_file,out_file) +\
            "-outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart " +\
            "qend sstart send \"")

    make_db(nt_file_a)
    make_db(nt_file_b)
    run_blast(nt_file_a+".db",nt_file_b,nt_file_b+".out")
    run_blast(nt_file_b+".db",nt_file_a,nt_file_a+".out")
    os.system("rm %s"%(nt_file_a+".db*"))
    os.system("rm %s"%(nt_file_b+".db*"))

def reciprocal_best_blastp(p_file_a,p_file_b,e_val):
    #first make 2 db
    #use blastp do the blast
    #generate 2 blast results
    def make_db(p_file):
        os.system("makeblastdb -in %s -dbtype prot -out %s"%(p_file, p_file+".db"))
    def run_blast(db,p_file,out_file):
        os.system("blastp -evalue %f -max_target_seqs 1 "%(e_val)+\
            "-db %s -query %s -out %s "%(db,p_file,out_file) +\
            "-outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart " +\
            "qend sstart send \"")

    make_db(p_file_a)
    make_db(p_file_b)
    run_blast(p_file_a+".db",p_file_b,p_file_b+".out")
    run_blast(p_file_b+".db",p_file_a,p_file_a+".out")
    os.system("rm %s"%(p_file_a+".db*"))
    os.system("rm %s"%(p_file_b+".db*"))

def write_out_put(fsa_file_1,fsa_file_2,final_out):
    sp12sp2 = {}
    for line in open(fsa_file_1+".out"):
        sp1,sp2 = line.split("\t")[:2]
        sp12sp2[sp1] = sp2

    sp22sp1 = {}
    for line in open(fsa_file_2+".out"):
        sp2,sp1 = line.split("\t")[:2]
        sp22sp1[sp2] = sp1

    f = open(final_out,"w")
    for sp1 in sp12sp2.keys():
        sp2 = sp12sp2[sp1]
        if sp2 in sp22sp1 and sp22sp1[sp2] == sp1:
            f.write("%s\t%s\n"%(sp1,sp2))

def main():
    check_file()
    if "-h" in sys.argv or "-help" in sys.argv:
        print_help()
        sys.exit()

    fsa_file_1 = sys.argv[-3]
    fsa_file_2 = sys.argv[-2]
    final_out = sys.argv[-1]
    e_val = 0.00001

    if "-e" in sys.argv:
        e_val = float(sys.argv[sys.argv.index("-e") + 1])
    if "-p" in sys.argv:
        reciprocal_best_blastp(fsa_file_1,fsa_file_2,e_val)
    else:
        reciprocal_best_blastn(fsa_file_1,fsa_file_2,e_val)
    
    write_out_put(fsa_file_1,fsa_file_2,final_out)

    if "-f" in sys.argv:
        pass
    else:
        os.system("rm %s"%fsa_file_1+".out")
        os.system("rm %s"%fsa_file_2+".out")

main()
