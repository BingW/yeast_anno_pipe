# coding:utf-8
# Author: bingwang
# Email: toaya.kase@gmail.com
# Copylight 2012-2012 Bing Wang
# LICENCES: GPL v3.0

__docformat__ = "epytext en"

suvaONscer = "/Users/bingwang/zen/yeast_anno_pipe/suvaONscer.txt"
seubONscer = "/Users/bingwang/zen/yeast_anno_pipe/seub3000ONscer.txt"
suva_seub_pair = "/Users/bingwang/zen/yeast_anno_pipe/suva_seub_pair.tab"

out_put_pillar = "/Users/bingwang/zen/yeast_anno_pipe/scer_suva_seub_pillar.tab"

class Gene():
    def __init__(self,block):
        self.name = block[0]
        self.seqid = block[1]
        self.start = int(block[2])
        self.end = int(block[3])
        self.strand = block[4]
        self.mapped = False
        self.conflict = ""
        self.block = "\t".join(block)

def remove_single_occurence_scaf(pillar):
    #remove seub and suva single occurence scaf
    #_pillar is hard copy of pillar
    _pillar = [[a,b,c] for [a,b,c] in pillar]
    suva_unmapped = []
    seub_unmapped = []
    for i,line in enumerate(_pillar):
        scer,suva,seub = line
        # remove suva's single occurence scaf
        if suva:
            suva_seqid = suva.split("0")[1]
            j = i-1
            up_suva = _pillar[j][1]
            while (not up_suva and j>0):
                up_suva = _pillar[j][1]
                j -= 1
            k = i +1
            down_suva = _pillar[k][1]
            while (not down_suva and k<len(_pillar)):
                down_suva = _pillar[k][1]
                k += 1
            if not up_suva:
                up_suva = "00"
            if not down_suva:
                down_suva = "00"
            if up_suva.split("0")[1] != suva.split("0")[1] != down_suva.split("0")[1]:
                _pillar[i][1] = ""
                suva_unmapped.append(suva)
                assert pillar[i] != _pillar[i]

        # remove seub's single occurence scaf
        if seub:
            seub_seqid = seub.split("0")[1]
            j = i-1
            up_seub = _pillar[j][2]
            while (not up_seub and j>0):
                up_seub = _pillar[j][2]
                j -= 1
            k = i +1
            down_seub = _pillar[k][2]
            while (not down_seub and k<len(_pillar)):
                down_seub = _pillar[k][2]
                k += 1
            if not up_seub:
                up_seub = "00"
            if not down_seub:
                down_seub = "00"
            if up_seub.split("0")[1] != seub.split("0")[1] != down_seub.split("0")[1]:
                _pillar[i][2] = ""
                seub_unmapped.append(seub)
                assert pillar[i] != _pillar[i]

    return _pillar,suva_unmapped,seub_unmapped

#-----get all dict initial------
scer_order = []
scer2record= {}
seub2record= {}
scer121seub = {}
seub121scer = {}
for line in open(seubONscer):
    line = line.replace("\n","").split("\t")
    #a block = [name,seqid,start,end,strand]
    scer,seub = line[0],line[5]
    if scer:
        scer_order.append(scer)
        scer2record[scer] = Gene(line[:5])
    if seub:
        seub2record[seub] = Gene(line[5:])
        if scer:
            scer121seub[scer] = seub
            seub121scer[seub] = scer

suva2record= {}
scer121suva = {}
suva121scer = {}
for line in open(suvaONscer):
    line = line.replace("\n","").split("\t")
    scer,suva = line[0],line[5]
    if suva:
        suva2record[suva] = Gene(line[5:])
        if scer:
            scer121suva[scer] = suva
            suva121scer[suva] = scer

suva121seub = {}
seub121suva = {}
for line in open(suva_seub_pair):
    line = line.replace("\n","").split("\t")
    suva,seub = line[0],line[1]
    suva121seub[suva] = seub
    seub121suva[seub] = suva

#---------Main start----------
pillar = []
#write 4 category's pillar:
#scer suva seub
#  +    +    +   possible conflict case "suva_seub_RBBH:suvaxx-seubxx VS YGAP anno:scerxx-suvaxx-seubxx" 
#  +    +    -   possible conflict case "suva_seub_RBBH:suvaxx-seubxx VS YGAP anno:scerxx-suvaxx-seub--"
#  +    -    +   possible conflict case "seub_suva_RBBH:seubxx-suvaxx VS YGAP anno:scerxx-suva---seubxx"
#  +    -    - 
for scer in scer_order:
    if scer in scer121suva:
        suva = scer121suva[scer]
        if scer in scer121seub:
            seub = scer121seub[scer]
            pillar.append([scer,suva,seub])
            seub2record[seub].mapped = True
            if suva in suva121seub and suva121seub[suva] != seub:
                conflict_info = "suva_seub_RBBH:%s-%s VS YGAP anno:%s-%s-%s"%(suva,suva121seub[suva],scer,suva,seub)
                print conflict_info
            if seub in seub121suva and seub121suva[seub] != suva:
                conflict_info = "suva_seub_RBBH:%s-%s VS YGAP anno:%s-%s-%s"%(seub121suva[seub],seub,scer,suva,seub)
                print conflict_info
        else:
            pillar.append([scer,scer121suva[scer],""])
            if suva in suva121seub:
                conflict_info = "suva_seub_RBBH:%s-%s VS YGAP anno:%s-%s-%s"%(suva,suva121seub[suva],scer,suva,"SEUB---")
                print conflict_info
        suva2record[scer121suva[scer]].mapped = True
    else:
        if seub in scer121seub:
            pillar.append([scer,"",scer121seub[scer]])
            if seub in seub121suva:
                conflict_info = "suva_seub_RBBH:%s-%s VS YGAP anno:%s-%s-%s"%(seub121suva[seub],seub,scer,"SUVA---",seub)
                print conflict_info
            seub2record[scer121seub[scer]].mapped = True
        else:
            pillar.append([scer,"",""])
    scer2record[scer].mapped = True

pillar,suva_unmap,seub_unmap = remove_single_occurence_scaf(pillar)
for suva in suva_unmap:
    suva2record[suva].mapped = False
for seub in seub_unmap:
    seub2record[seub].mapped = False

#insert suva specific (to scer) genes continues region:
#scer  suva  seub
# -     +     +
# -     +     -

new_pillar = []
for i,line in enumerate(pillar):
    #append a hard copy of line
    new_pillar.append([gene for gene in line])
    now_suva = line[1]
    if now_suva:
        j = i-1
        up_suva = pillar[j][1]
        while (not up_suva and j>0):
            up_suva = pillar[j][1]
            j -= 1
        k = i +1
        down_suva = pillar[k][1]
        while (not down_suva and k<len(pillar)):
            down_suva = pillar[k][1]
            k += 1
        if not up_suva:
            up_suva = "0000"
        if not down_suva:
            down_suva = "0000"
        up_scaf = up_suva.split("0")[1]
        up_pos = int(up_suva[-4:])
        now_scaf = now_suva.split("0")[1]
        now_pos = int(now_suva[-4:])
        down_scaf = down_suva.split("0")[1]
        down_pos = int(down_suva[-4:])
        assert now_suva == "0".join(["SUVA",now_scaf,"0"*(4-len(str(now_pos)))+str(now_pos)])
        if now_scaf == up_scaf == down_scaf:
            #in this case: 
            #A1  <-  upper
            #A3  <-  now
            #A..  <-  !fill these out
            #A5  <-  down
            if down_pos - now_pos > 10:
                while down_pos - now_pos > 10:
                    now_pos += 10
                    possible_gene = "0".join(["SUVA",now_scaf,"0"*(4-len(str(now_pos)))+str(now_pos)])
                    if possible_gene in suva2record:
                        if suva2record[possible_gene].mapped:
                            print "check @ %s, the possible up_stream gene %s is already mapped"%(down_suva,possible_gene)
                            break
                        else:
                            seub = "" if possible_gene not in suva121seub else suva121seub[possible_gene]
                            seub = seub if (seub in seub2record and not seub2record[seub].mapped) else ""
                            new_pillar.append(["",possible_gene,seub])
                            suva2record[possible_gene].mapped = True
                            if seub == "SEUB0AJ00200":
                                print seub2record[seub].mapped
                            if seub in seub2record:
                                seub2record[seub].mapped = True

            elif now_pos - down_pos > 10:
                while now_pos - down_pos > 10:
                    now_pos -= 10
                    possible_gene = "0".join(["SUVA",now_scaf,"0"*(4-len(str(now_pos)))+str(now_pos)])
                    if possible_gene in suva2record:
                        if suva2record[possible_gene].mapped:
                            print "check @ %s, the possible up_stream gene %s is already mapped"%(down_suva,possible_gene)
                            break
                        else:
                            seub = "" if possible_gene not in suva121seub else suva121seub[possible_gene]
                            seub = seub if (seub in seub2record and not seub2record[seub].mapped) else ""
                            new_pillar.append(["",possible_gene,seub])
                            suva2record[possible_gene].mapped = True
                            if seub in seub2record:
                                seub2record[seub].mapped = True
            
        elif (now_scaf == up_scaf and now_scaf != down_scaf):
            #in this case:
            #A3  <-  upper
            #A5  <-  now
            #A.. <-  ! fill these out
            #B.. <-  ! leave it for case 3
            #B2  <-  down
            if up_pos > now_pos:
                while now_pos > 100:
                    now_pos -= 10
                    possible_gene = "0".join(["SUVA",now_scaf,"0"*(4-len(str(now_pos)))+str(now_pos)])
                    if possible_gene in suva2record:
                        if suva2record[possible_gene].mapped:
                            break
                        else:
                            seub = "" if possible_gene not in suva121seub else suva121seub[possible_gene]
                            seub = seub if (seub in seub2record and not seub2record[seub].mapped) else ""
                            new_pillar.append(["",possible_gene,seub])
                            suva2record[possible_gene].mapped = True
                            if seub in seub2record:
                                seub2record[seub].mapped = True

            elif now_pos > up_pos:
                while now_pos < 9990:
                    now_pos += 10
                    possible_gene = "0".join(["SUVA",now_scaf,"0"*(4-len(str(now_pos)))+str(now_pos)])
                    if possible_gene in suva2record:
                        if suva2record[possible_gene].mapped:
                            break
                        else:
                            seub = "" if possible_gene not in suva121seub else suva121seub[possible_gene]
                            seub = seub if (seub in seub2record and not seub2record[seub].mapped)else ""
                            new_pillar.append(["",possible_gene,seub])
                            suva2record[possible_gene].mapped = True
                            if seub in seub2record:
                                seub2record[seub].mapped = True

        elif (now_scaf != up_scaf and now_scaf == down_scaf):
            #in this case:
            #A5  <-  upper
            #B.. <-  ! fill these out first
            #B2  <-  now
            #B.. <-  ! fill these out second
            #B4  <-  down
            if down_pos > now_pos:
                temp_pos = now_pos
                p = -1
                while temp_pos > 100:
                    temp_pos -= 10
                    possible_gene = "0".join(["SUVA",now_scaf,"0"*(4-len(str(temp_pos)))+str(temp_pos)])
                    if possible_gene in suva2record:
                        if suva2record[possible_gene].mapped:
                            break
                        else:
                            seub = "" if possible_gene not in suva121seub else suva121seub[possible_gene]
                            seub = seub if (seub in seub2record and not seub2record[seub].mapped) else ""
                            new_pillar.insert(p,["",possible_gene,seub])
                            suva2record[possible_gene].mapped = True
                            if seub in seub2record:
                                seub2record[seub].mapped = True
                            p -= 1

                while down_pos - now_pos > 10:
                    now_pos += 10
                    possible_gene = "0".join(["SUVA",now_scaf,"0"*(4-len(str(now_pos)))+str(now_pos)])
                    if possible_gene in suva2record:
                        if suva2record[possible_gene].mapped:
                            print "check @ %s, the possible up_stream gene %s is already mapped"%(down_suva,possible_gene)
                            break
                        else:
                            seub = "" if possible_gene not in suva121seub else suva121seub[possible_gene]
                            seub = seub if (seub in seub2record and not seub2record[seub].mapped) else ""
                            new_pillar.append(["",possible_gene,seub])
                            suva2record[possible_gene].mapped = True
                            if seub in seub2record:
                                seub2record[seub].mapped = True

            elif down_pos < now_pos:
                temp_pos = now_pos                
                p = -1
                while temp_pos < 9900:
                    temp_pos += 10
                    possible_gene = "0".join(["SUVA",now_scaf,"0"*(4-len(str(temp_pos)))+str(temp_pos)])
                    if possible_gene in suva2record:
                        if suva2record[possible_gene].mapped:
                            break
                        else:
                            seub = "" if possible_gene not in suva121seub else suva121seub[possible_gene]
                            seub = seub if (seub in seub2record and not seub2record[seub].mapped) else ""
                            new_pillar.insert(p,["",possible_gene,seub])
                            suva2record[possible_gene].mapped = True
                            if seub in seub2record:
                                seub2record[seub].mapped = True
                            p -= 1

                while now_pos - down_pos > 10:
                    now_pos -= 10
                    possible_gene = "0".join(["SUVA",now_scaf,"0"*(4-len(str(now_pos)))+str(now_pos)])
                    if possible_gene in suva2record:
                        if suva2record[possible_gene].mapped:
                            print "check @ %s, the possible up_stream gene %s is already mapped"%(down_suva,possible_gene)
                            break
                        else:
                            seub = "" if possible_gene not in suva121seub else suva121seub[possible_gene]
                            seub = seub if (seub in seub2record and not seub2record[seub].mapped) else ""
                            new_pillar.append(["",possible_gene,seub])
                            suva2record[possible_gene].mapped = True
                            if seub in seub2record:
                                seub2record[seub].mapped = True

        else:
            assert 1 == 0

pillar,suva_unmap,seub_unmap = remove_single_occurence_scaf(new_pillar)
for seub in seub_unmap:
    seub2record[seub].mapped = False

#insert seub specific (to scer) genes continues region:
#scer  suva  seub
# -     -     +

new_pillar = []
for i,line in enumerate(pillar):
    #append a hard copy of line
    new_pillar.append([gene for gene in line])
    now_seub = line[2]
    if now_seub:
        j = i-1
        up_seub = pillar[j][2]
        while (not up_seub and j>0):
            up_seub = pillar[j][2]
            j -= 1
        k = i +1
        down_seub = pillar[k][2]
        while (not down_seub and k<len(pillar)):
            down_seub = pillar[k][2]
            k += 1
        if not up_seub:
            up_seub = "0000"
        if not down_seub:
            down_seub = "0000"
        up_scaf = up_seub.split("0")[1]
        up_pos = int(up_seub[-4:])
        now_scaf = now_seub.split("0")[1]
        now_pos = int(now_seub[-4:])
        down_scaf = down_seub.split("0")[1]
        down_pos = int(down_seub[-4:])
        assert now_seub == "0".join(["SEUB",now_scaf,"0"*(4-len(str(now_pos)))+str(now_pos)])
        if now_scaf == up_scaf == down_scaf:
            #in this case: 
            #A1  <-  upper
            #A3  <-  now
            #A..  <-  !fill these out
            #A5  <-  down
            if down_pos - now_pos > 10:
                while down_pos - now_pos > 10:
                    now_pos += 10
                    possible_gene = "0".join(["SEUB",now_scaf,"0"*(4-len(str(now_pos)))+str(now_pos)])
                    if possible_gene in seub2record:
                        if seub2record[possible_gene].mapped:
                            print "check @ %s, the possible up_stream gene %s is already mapped"%(down_seub,possible_gene)
                            break
                        else:
                            new_pillar.append(["","",possible_gene])
                            seub2record[possible_gene].mapped = True

            elif now_pos - down_pos > 10:
                while now_pos - down_pos > 10:
                    now_pos -= 10
                    possible_gene = "0".join(["SEUB",now_scaf,"0"*(4-len(str(now_pos)))+str(now_pos)])
                    if possible_gene in seub2record:
                        if seub2record[possible_gene].mapped:
                            print "check @ %s, the possible up_stream gene %s is already mapped"%(down_seub,possible_gene)
                            break
                        else:
                            new_pillar.append(["","",possible_gene])
                            seub2record[possible_gene].mapped = True
            
        elif (now_scaf == up_scaf and now_scaf != down_scaf):
            #in this case:
            #A3  <-  upper
            #A5  <-  now
            #A.. <-  ! fill these out
            #B.. <-  ! leave it for case 3
            #B2  <-  down
            if up_pos > now_pos:
                while now_pos > 100:
                    now_pos -= 10
                    possible_gene = "0".join(["SEUB",now_scaf,"0"*(4-len(str(now_pos)))+str(now_pos)])
                    if possible_gene in seub2record:
                        if seub2record[possible_gene].mapped:
                            break
                        else:
                            new_pillar.append(["","",possible_gene])
                            seub2record[possible_gene].mapped = True

            elif now_pos > up_pos:
                while now_pos < 9990:
                    now_pos += 10
                    possible_gene = "0".join(["SUVA",now_scaf,"0"*(4-len(str(now_pos)))+str(now_pos)])
                    if possible_gene in seub2record:
                        if seub2record[possible_gene].mapped:
                            break
                        else:
                            new_pillar.append(["","",possible_gene])
                            seub2record[possible_gene].mapped = True

        elif (now_scaf != up_scaf and now_scaf == down_scaf):
            #in this case:
            #A5  <-  upper
            #B.. <-  ! fill these out first
            #B2  <-  now
            #B.. <-  ! fill these out second
            #B4  <-  down
            if down_pos > now_pos:
                temp_pos = now_pos
                p = -1
                while temp_pos > 100:
                    temp_pos -= 10
                    possible_gene = "0".join(["SEUB",now_scaf,"0"*(4-len(str(temp_pos)))+str(temp_pos)])
                    if possible_gene in seub2record:
                        if seub2record[possible_gene].mapped:
                            break
                        else:
                            new_pillar.append(["","",possible_gene])
                            seub2record[possible_gene].mapped = True
                            p -= 1

                while down_pos - now_pos > 10:
                    now_pos += 10
                    possible_gene = "0".join(["SEUB",now_scaf,"0"*(4-len(str(now_pos)))+str(now_pos)])
                    if possible_gene in seub2record:
                        if seub2record[possible_gene].mapped:
                            print "check @ %s, the possible up_stream gene %s is already mapped"%(down_seub,possible_gene)
                            break
                        else:
                            new_pillar.append(["","",possible_gene])
                            seub2record[possible_gene].mapped = True

            elif down_pos < now_pos:
                temp_pos = now_pos                
                p = -1
                while temp_pos < 9900:
                    temp_pos += 10
                    possible_gene = "0".join(["SEUB",now_scaf,"0"*(4-len(str(temp_pos)))+str(temp_pos)])
                    if possible_gene in seub2record:
                        if seub2record[possible_gene].mapped:
                            break
                        else:
                            new_pillar.insert(p,["","",possible_gene])
                            seub2record[possible_gene].mapped = True
                            p -= 1

                while now_pos - down_pos > 10:
                    now_pos -= 10
                    possible_gene = "0".join(["SEUB",now_scaf,"0"*(4-len(str(now_pos)))+str(now_pos)])
                    if possible_gene in seub2record:
                        if seub2record[possible_gene].mapped:
                            print "check @ %s, the possible up_stream gene %s is already mapped"%(down_seub,possible_gene)
                            break
                        else:
                            new_pillar.append(["","",possible_gene])
                            seub2record[possible_gene].mapped = True
        else:
            assert 1 == 0

content = []
for line in new_pillar:
    scer_block = "\t"*4 if line[0] not in scer2record else scer2record[line[0]].block
    suva_block = "\t"*4 if line[1] not in suva2record else suva2record[line[1]].block
    seub_block = "\t"*4 if line[2] not in seub2record else seub2record[line[2]].block
    content.append("\t".join([scer_block,suva_block,seub_block]))
open(out_put_pillar,"w").write("\n".join(content))
