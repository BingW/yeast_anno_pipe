#! /usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"


"""
This program will take two nt fasta file as input, 
use balst one to one pair.
>>> import compare2fsa
>>> compare2fsa.main(fsafile_1,fsa_file_2)
...
"""
from Bio import SeqIO
import os
os.chdir("/Users/bingwang/zen/yeast_anno_pipe/")
#TODO change dir need to be more elegant

def check_file(file_name):
    '''
    check if a file is fsa format
    >>> check_file(test_1.fsa)
    True
    >>> check_file(test_wrong.fsa)
    False
    '''
    pass

def db_construct(file_name):
    '''
    make a blast datebase
    >>> db_construct("/test/test_1.fsa")
    >>> os.path.isfile("/tests/test_1.fsa.db")
    True
    '''
    db_name = file_name + ".db"
    os.system("makeblastdb -in %s -dbtype nucl -out %s"%(fsa,db))

def run_blast(file_name,db_name):
    '''
    run blast program
    >>> run_blast("/test/test_1.fsa","/test/test_2.fsa.db")
    >>> os.path.isfile("/test/test_1.fsa.out")
    True
    '''
    blast_result_file = file_name + ".out"
    os.system("blastn -evalue 0.00001 -max_target_seqs 1 -strand plus "+\
        "-max_hsps_per_subject 1 " +\
        "-db %s -query %s -out %s "%(db_name,file_name,blast_result_file) +\
        "-outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart " +\
        "qend sstart send \"")

def main():
    check_file(file_1)
    db_construct(file_1)
    db_construct(file_2)
    blast(file_1,file_2)
    rm_db(file_1,file_2)
    blast_1 = read_blast(file_1)
    blast_2 = read_blast(file_2)
    pairs = pair(blast_1,blast_2)
    write_pairs(file_out)

if __name__ == "__main__":
    file_1 = open("tests/scer.devin.fsa")
    file_2 = open("tests/scer.ygap.fsa")
    file_out = open("test/scer/ygapVSdevin.tab")
    main(file_1,file_2,file_out)


