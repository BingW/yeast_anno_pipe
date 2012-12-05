data = open("/Users/bingwang/zen/yeast_anno_pipe/SSS_dataset/Scer.gff").read().replace("\t",",")
g = open("/Users/bingwang/zen/yeast_anno_pipe/output/Scer_devin.csv","w").write(data)

data = open("/Users/bingwang/zen/yeast_anno_pipe/YGAP_predict/scer.YGAP.txt").read().replace("\t",",")
g = open("/Users/bingwang/zen/yeast_anno_pipe/output/Scer_ygap.csv","w").write(data)
