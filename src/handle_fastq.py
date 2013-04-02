#coding: utf-8
import time
fastq_file = "/Users/bingwang/zen/yeast_anno_pipe/s_1_1_sequence.txt.FM1318.cut.fastq"
start = time.time()

f = open(fastq_file)
while 1:
    f.readline(100)
    pass

print time.time() - start

