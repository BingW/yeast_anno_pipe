#!/usr/bin/env python
# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"

import os
import sys

#datapath = sys.argv[1]
datapath = "/Users/bingwang/zen/assemble_example/"

lib_dype_dict = {1:"pair end",2:"mate pair",3:"single"}
data_format_dict = {1:"two paired fastq files", 2:"two paired fasta files", 3:"single fastq file", \
        4:"single fasta file", 5:"paired fasta in single file", 6:"bam file"}

class Lib():
    def __init__(self,lib_type,lib_format):
        self.avg_ins=500 #average insert size
        self.asm_flags=3
        #1 only contig assembly
        #2 only scaffold assembly
        #3 both contig and scaffold assembly
        #4 only gap closure
        self.rd_len_cutoff=100 #use only first 100 bps of each read
        self.rank=1 #in which order the reads are used while scaffolding
        if lib_type != 2:
            self.reverse_seq= 0 #0 for pair end, 1 for mate pair
            self.pair_num_cutoff=3 # The minimum number for paired-end reads and mate-pair reads is 3 and 5 respectively.
            self.map_len=32
        elif lib_type == 2:
            self.reverse_seq= 0 #0 for pair end, 1 for mate pair
            self.pair_num_cutoff=3 # The minimum number for paired-end reads and mate-pair reads is 3 and 5 respectively.
            self.map_len=35
        if lib_format == 1:
            self.sequence = {"q1":"","q2":""}
        elif lib_format == 2:
            self.sequence = {"f1":"","f2":""}
        elif lib_format == 3:
            self.sequence = {"q":""}
        elif lib_format == 4:
            self.sequence = {"f":""}
        elif lib_format == 5:
            self.sequence = {"p":""}
        elif lib_format == 6:
            self.sequence = {"b":""}
        else:
            print "lib_format should be int druing [1-6]"

class SoapParameters():
    def __init__(self):
        self.max_rd_len = 100  #maximal read length
        self.s = datapath + "soap.config" #configFile: the config file of solexa reads
        self.o = "outputGraph" # outputGraph: prefix of output graph file name
        self.K = 23 # kmer(min 13, max 63/127): kmer size
        self.p = 8 # n_cpu: number of cpu for use
        self.a = 0 # initMemoryAssumption memory assumption initialized to avoid further reallocation, unit G
        self.d = 0 # KmerFreqCutoff: kmers with frequency no larger than KmerFreqCutoff will be deleted
        self.R = "OFF" #resolve repeats by reads
        self.D = 1 # EdgeCovCutoff: edges with coverage no larger than EdgeCovCutoff will be deleted
        self.M = 1 # mergeLevel(min 0, max 3): the strength of merging similar sequences during contiging
        self.m = 127 # max k when using multi kmer
        self.e = 0 # weight to filter arc when linearize two edges
        self.r = "OFF" # keep avalable reads (*.read)
        self.E = "OFF" # merge clean bubble before iterate
        self.f = "OFF" # output gap related reads in map step for using SRkgf to fill gap
        self.k = self.K # kmer_R2C(min 13, max 63): kmer size used for mapping read to contig
        self.F = "OFF" # fill gaps in scaffold
        self.u = "OFF" # un-mask contigs with high/low coverage before scaffolding
        self.w = "OFF" # keep contigs weakly connected to other contigs in scaffold
        self.G = 50 # gapLenDiff: allowed length difference between estimated and filled gap
        self.L = self.K+2 # minContigLen: shortest contig for scaffolding
        self.c = 0.1 # minContigCvg: minimum contig coverage (c*avgCvg), contigs shorter than 100bp with coverage smaller 
                     # than c*avgCvg will be masked before scaffolding unless -u is set
        self.C = 2 # maxContigCvg: maximum contig coverage (C*avgCvg), contigs with coverage larger than C*avgCvg or 
                   # contigs shorter than 100bp with coverage larger than 0.8*C*avgCvg will be masked before scaffolding
                   # unless -u is set
        self.b = 1.5 # insertSizeUpperBound: (b*avg_ins) will be used as upper bound of insert size for large insert size
                     # ( > 1000) when handling pair-end connections between contigs if b is set to larger than 1
        self.B = 0.6 # bubbleCoverage: remove contig with lower cvoerage in bubble structure if both contigs coverage are 
                     # smaller than bubbleCoverage*avgCvg
        self.N = 0 # genomeSize: genome size for statistics
        self.V = "OFF" # output visualization information of assembly
        self.lib = []

    def add_lib(self,library):
        self.lib.append(library)
    
    def write(self,handle):
        content = []
        content.append("#maximal read length")
        content.append("max_rd_len=%d"%self.max_rd_len)
        for lib_record in self.lib:
            content.append("[LIB]")
            content.append("#average insert size")
            content.append("avg_ins=%d"%lib_record.avg_ins)
            content.append("#assemble flags")
            content.append("#1 only contig assembly")
            content.append("#2 only scaffold assembly")
            content.append("#3 both contig and scaffold assembly")
            content.append("#4 only gap closure")
            content.append("asm_flags=%d"%lib_record.asm_flags)
            content.append("#use only first 100 bps of each read")
            content.append("rd_len_cutoff=%d"%lib_record.rd_len_cutoff)
            content.append("#in which order the reads are used while scaffolding")
            content.append("rank=%d"%lib_record.rank)
            content.append("#0 for pair end, 1 for mate pair")
            content.append("reverse_seq=%d"%lib_record.reverse_seq)
            content.append("#The minimum number for paired-end reads and mate-pair is 3 and 5 respectively")
            content.append("pair_num_cutoff=%d"%lib_record.pair_num_cutoff)
            content.append("#minimum aligned length to contigs for a reliable read location")
            content.append("map_len=%d"%lib_record.map_len)
            content.append("#Please modify these values to your file location")
            for seq_id in lib_record.sequence:
                content.append("%s=%s"%(seq_id,lib_record.sequence[seq_id]))
        content.append("# Below is for command line")
        content.append("# kmer(min 13, max 63/127): kmer size")
        content.append("-K=%d"%self.K)
        content.append("# n_cpu: number of cpu for use")
        content.append("-p=%d"%self.p)
        content.append("# initMemoryAssumption memory assumption initialized to avoid further reallocation, unit G")
        content.append("-a=%d"%self.a)
        content.append("# KmerFreqCutoff: kmers with frequency no larger than KmerFreqCutoff will be deleted")
        content.append("-d=%d"%self.d)
        content.append("# resolve repeats by reads")
        content.append("-R=%s"%self.R)
        content.append("# EdgeCovCutoff: edges with coverage no larger than EdgeCovCutoff will be deleted")
        content.append("-D=%d"%self.D)
        content.append("# mergeLevel(min 0, max 3): the strength of merging similar sequences during contiging")
        content.append("-M=%d"%self.M)
        content.append("# max k when using multi kmer")
        content.append("-m=%d"%self.m)
        content.append("# weight to filter arc when linearize two edges")
        content.append("-e=%d"%self.e)
        content.append("# keep avalable reads (*.read)")
        content.append("-r=%s"%self.r)
        content.append("# merge clean bubble before iterate")
        content.append("-E=%s"%self.E)
        content.append("# output gap related reads in map step for using SRkgf to fill gap")
        content.append("-f=%s"%self.f)
        content.append("# kmer_R2C(min 13, max 63): kmer size used for mapping read to contig")
        content.append("-k=%d"%self.k)
        content.append("# fill gaps in scaffold")
        content.append("-F=%s"%self.F)
        content.append("# un-mask contigs with high/low coverage before scaffolding")
        content.append("-u=%s"%self.u)
        content.append("# keep contigs weakly connected to other contigs in scaffold")
        content.append("-w=%s"%self.w)
        content.append("# gapLenDiff: allowed length difference between estimated and filled gap")
        content.append("-G=%d"%self.G)
        content.append("# minContigLen: shortest contig for scaffolding")
        content.append("-L=%d"%self.L)
        content.append("# minContigCvg: minimum contig coverage (c*avgCvg), contigs shorter than 100bp with coverage smaller")
        content.append("# than c*avgCvg will be masked before scaffolding unless -u is set")
        content.append("-c=%f"%self.c)
        content.append("# maxContigCvg: maximum contig coverage (C*avgCvg), contigs with coverage larger than C*avgCvg or")
        content.append("# contigs shorter than 100bp with coverage larger than 0.8*C*avgCvg will be masked before scaffolding")
        content.append("# unless -u is set")
        content.append("-C=%d"%self.C)
        content.append("# insertSizeUpperBound: (b*avg_ins) will be used as upper bound of insert size for large insert size")
        content.append("# ( > 1000) when handling pair-end connections between contigs if b is set to larger than 1")
        content.append("-b=%f"%self.b)
        content.append("# bubbleCoverage: remove contig with lower cvoerage in bubble structure if both contigs coverage are")
        content.append("# smaller than bubbleCoverage*avgCvg")
        content.append("-B=%f"%self.B)
        content.append("# genomeSize: genome size for statistics")
        content.append("-N=%d"%self.N)
        content.append("# output visualization information of assembly")
        content.append("-V=%s"%self.V)
        handle.write("\n".join(content))
        handle.close()
    
     def parse(self,handle):
         for line in handle:
             if not line.startswith("#") and "=" in line:
     
     def run_cmd(self):
         pass


def print_usage():
    s = """
    This program use SOAPdenove to assemble reads.
    
    1. To use this program you need Installed SOAPdenovo which includes SOAPdenovo127mer
    
    2. Simply type cmd in terminal

    $bing_assemble_soapdenovo.py <your work path>
    
    3. The program will generate a config file in your path, check and modify it. Then run the commendline again
    
    $bing_assemble_soapdenovo.py <your work path>
    """
    print s

def main():
    config = SoapParameters()

    if not os.path.exists(datapath):
        print_usage()
        sys.exit("\033[91m Error, path does not exists\033[0m")

    elif "project.config" not in os.listdir(datapath):
        total_lib_num = raw_input('Please tell me the how many libraries are there.[default:1]')
        try:
            total_lib_num = int(total_lib_num)
        except:
            total_lib_num = 1
        print "Library number set as %d\n"%(total_lib_num)
        
        for lib_id in range(total_lib_num):

            lib_type = raw_input("""
Please choose your lib type: 
[1]:pair end 
[2]:mate pair 
[3]:single.
[default:1]""")
            try:
                lib_type = int(lib_type)
            except:
                lib_type = 1
            print "The %d lib type set as %s\n"%(lib_id+1,lib_dype_dict[lib_type])

            data_format = raw_input("""
Please choose your lib data format: 
[1]:two paired fastq files
[2]:two paired fasta files
[3]:single fastq file 
[4]:single fasta file
[5]:paired fasta in single file
[6]:bam file
[default:1]""")
            try:
                data_format = int(data_format)
                assert data_format in range(1,7,1)
            except:
                data_format = 1
            print "The %d lib data type set as %s\n"%(lib_id+1,data_format_dict[data_format])

            lib_record = Lib(lib_type,data_format)
            config.add_lib(lib_record)
        config.write(open(datapath+"project.config","w"))
        print "project.config already generated"
        sys.exit("Please check config parameters in your path")
    else:
        total_lib_num = 3
        config.parse(open(datapath+"project.config"))
        run_cmd(config)

if __name__ == "__main__":
    main()
