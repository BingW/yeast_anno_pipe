f = open("/Users/bingwang/zen/yeast_anno_pipe/YGOB_dataset/NT.fsa")
descrip = {}
for line in f:
    if line.startswith(">"):
        name = line.split(" ")[0][1:]
        description = line[line.find("}")+1:-1]
        descrip[name] = description
for record in descrip:
    print record,descrip[record]
