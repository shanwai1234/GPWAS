args = commandArgs(trailingOnly=TRUE)
ori = read.table(args[1],sep=' ',head=F) # P-value generated from original phenotype data in GPWAS
perm = read.table(args[2],sep=' ',head=F) # P-value generated from permuted phenotype data in GPWAS
ori_p = sort(ori$V3)
perm_p = sort(perm$V3)
loop_number = 30 # picking a loop number for iterating p-value threshold
step = 2 # picking any integer number, the smaller value the narrower range of p-value threshold you can choose
print ("P-value threshold; FDR; Filtered Gene Number")
for (i in c(1:loop_number)){
  thres = step**(-1*i) # iterating different p-value threshold 
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
  fdr = (perm_count/as.integer(args[3]))/ori_count
  print (paste(thres,fdr,ori_count,sep=';'))
}
