ori = read.table("/your/path/original-pvalue-GPWAS-output.txt",sep=' ',head=F) # P-value generated from original phenotype data in GPWAS
perm = read.table("/your/path/merged-pvalue-GPWAS-outout.txt",sep=' ',head=F) # P-value generated from permuted phenotype data in GPWAS
ori_p = sort(ori$V3)
perm_p = sort(perm$V3)
start = 0.1 # start of p-value threshold
loop_number = 30 # picking a loop number for iterating p-value threshold
step = 2 # picking any integer number, the smaller value the narrower range of p-value threshold you can choose
print ("P-value threshold; FDR; Filtered Gene Number")
for (i in c(1:30)){ 
  thres = step**(-1*i)
  ori_count = 0
  perm_count = 0
  for (o in ori_p){
    if (o < thres){
      ori_count = ori_count + 1
    }
  }
  for (p in perm_p){
    if (p < thres){
      perm_count = perm_count + 1
    }
  }
  print (paste(thres,fdr,ori_count,sep=';'))
}
