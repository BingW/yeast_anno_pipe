#class I
#scer_chrI and suva_chr_1
seub_order =  [139, 91, 95,217,176,172]
seub_strand = ["-","-","+","-","-","-"]
#scer_chrIII and suva_chr_3 [1190-1560,100-1110],[+,+]
#scaf_29 is continues when suva_chr_3 reverse, but only one gene's position
seub_order =  [170, 96,132, 29,147, 53,236,103]
seub_strand = ["-","-","+","-","-","-","-","+"]
#scer_chrV and suva_chr_5 [120-500,640-520,660-2670],[+,-,+]
#scaf_15 observed the same during suva's 120-720's region
seub_order =  [ 15, 97, 25,193,180,181, 60,213]
seub_strand = ["+","-","-","-","-","-","+","-"]
#scer_chrIX and suva_chr_9
seub_order =  [244, 74,191,117,116,111]
seub_strand = ["+","+","+","+","+","+"]
#scer_chrXI and suva_chr_11
seub_order =  [146,109,110, 16, 41,232,206, 5]
seub_strand = ["-","-","+","+","+","-","+","+"]
#scer_chrXII and suva_chr_10 #strange chr no.
seub_order =  [134, 31, 59,179, 24, 42, 17,157, 65,190,149, 89,200,125, 46, 84,220]
seub_strand = ["-","+","-","-","-","+","-","-","-","-","-","+","-","-","-","+","+"]
#scer_chrXIII and suva_chr_13
seub_order =  [152, 32, 64,216, 1 ,107,143,173, 21, 14, 57,162, 23,158]
seub_strand = ["-","+","-","+","-","-","+","+","-","+","+","+","-","+"]
#scer_chrXIV and suva_chr_14
#seub_scaf_104 and suva_chr_14(2830 - 2790) small reverse
seub_order =  [ 78,121, 45, 79, 83,104,154,227,119]
seub_strand = ["-","+","-","+","+","-","-","-","+"]

#class II
#scer_chrVII and suva_chr_7
#seub_scaf_215 and seub_chr_199 appear > 1, no more evidence for exchange between chromosomes
#exclude scaf_215, keep scaf_199
seub_order =  [234,140,186, 50,228, 9 , 51,197, 68,241,199,148,243,239, 72, 85,198,224,150]
seub_strand = ["+","+","+","+","-","+","-","+","+","-","+","-","+","+","+","+","-","+","+"]

#scer_chrXVI and suva_chr_16
#seub Scaffold_199 appears twice, no more evidence for exchange between chromosome
#since this one is shor, delete this one
seub_order  = [184,163, 6 , 33, 69, 43,214, 52, 81, 67, 94,128, 28,106,233, 63,192, 56,155,230]
seub_strand = ["-","+","+","+","-","-","-","+","-","-","-","-","+","-","-","+","+","-","+","+"]

#class III

#scer_chrII & scer_chrIV VS suva_scaf_2 & suva_scaf_4
#both of the watershed also appear in seub, so trust suva's seperation
#no scaffolds appears in other place
seub_chr2   = [ 34,118,183, 11,212, 49,204, 61, 2 ,159,171, 38,177,237,196,156, 87,144, 27, 70,137,102,122,115,120,114, 37]
seub_strand = ["-","-","+","-","-","+","-","+","+","-","-","-","-","+","+","+","-","-","+","-","-","+","+","-","-","+","+"] 
seub_chr4   = [127,135, 86,165,166, 8 , 47,226,138, 35, 12, 88,218,168, 82, 30,169, 22, 92, 10, 39, 58, 73]
seub_strand = ["-","-","-","+","+","+","-","-","+","+","-","-","-","-","-","+","+","-","+","+","-","-","-"]

#scer_chrX & scer_chrVI VS suva_scaf_12 & suva_scaf_6
#scaf_126 & scaf_40 inter exchange
#scaf_108 and scaf_215 continues where suva parted trust scer's genome
seub_chr6   = [ 20,242,209,126, 40,112,108]
seub_strand = ["+","-","+","-","-","-","-"]
seub_chr12  = [ 90,219,240,215,160, 54, 77,238, 13, 99]
seub_strand = ["+","-","+","+","+","+","+","+","+","+"]

#scer_chrXV & scer_chrVIII VS suva_scaf_15 & suva_scaf_8
#scaf_208 & scaf_188 recipical part, trust suva's partition
#scaf_203 I just keep chr_15's tail's part
seub_chr8  =  [ 3 , 48,188,225,164, 76,205,185,189,207,195]
seub_strand = ["-","-","-","-","-","-","+","-","-","+","+"]

seub_chr15 =  [100, 4 ,223,124,208, 44,123,131,101,211, 55,194,203]
seub_strand = ["+","-","+","-","+","+","-","-","+","-","-","-","+"]

