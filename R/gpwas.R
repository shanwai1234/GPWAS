#' Genome-Phenome Wide Association Study
#' @param ingeno Input genotype file name/directory.
#' @param inpheno Input phenotype file name/directory.
#' @param inpc Input folder with PCA parsed population structure covariance. If n number of chromosomes, n number of separate files should be included, as SNPs on each chromosome is excluded for performing PCA once.
#' @param g A list of specific gene that needs to analysis. By default the model will run for all of genes detected in the input genotype file.
#' @param gp Output file name/directory for selected phenotypes with every gene as well as p value of each selected phenotype.
#' @param gv Output file name/directory of terminated p value for each gene.
#' @param R Number of iteration for scanning all of input phenotypes with one specific gene.
#' @param pc Number of principle components that need to be included as covariances. Default is 3.
#' @param selectIn P value threshold for considering a phenotype is significant to be included into the model. Default is 0.01.
#' @param selectOut P value threshold for considering a phenotype is significant to be removed from the model. Default is 0.01.
#'
#' @import MASS leaps
#' @return Two files with significance level per gene and selected phenotypes per gene
#' @export
#'
#' @examples gpwas(ingeno,inpheno,inpc,g,gp,gv,R=35)
gpwas = function(ingeno,inpheno,inpc,g,gp,gv,R=num,pc=3,selectIn=0.01,selectOut=0.01){
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)
  mygeno = read.table(ingeno,sep='\t',head=T, check.names = FALSE)
  if (missing(g)){
    mygeno = mygeno
  } else {
    gfile = read.table(g,header=T)
    gset = as.character(gfile$Candidate)
    mygeno = mygeno[mygeno$Gene %in% gset, ]
  }
  mypheno = read.table(inpheno, head = TRUE)
  gene = unique(mygeno$Gene)
  R = R
  glist = c()
  gvalue = c()

  plist = c()
  for (i in c(1:pc)){
    plist = c(plist,paste('PC[, ',toString(i),']',sep=''))
  }
  init_pc = paste(plist,collapse= " + ")
  med_pc = paste(init_pc,' + ',sep='')

  for (g in gene){
    print (g)
    glist = c(glist,g)
    tp = mygeno[mygeno$Gene==g,]
    th = tp$SNP[1]
    chr = unlist(strsplit(as.character(th),'_'))[1]
    chrom = gsub('S','',chr)
    PC = read.table(paste(inpc,'/Population_structure_exclude_chrom_',chrom,'.txt',sep=''),sep=' ',head=T)
    PC = PC[,c(1:pc)]
    #--------------- Data processing ---------------
    A = tp[,-1]
    B1 = mypheno[, -1]
    M = dim(B1)
    C = names(A)
    C1 = intersect(C, as.character(mypheno[, 1]))
    N = length(C1)
    DataGene = matrix(0, N, dim(A)[1])
    DataPheno = matrix(0, N, M[2])
    #-------------- Missing data cleaning -------------
    for (i in 1 : N){
      index1 = which(C == C1[i])
      DataGene[i, ] = A[, index1]
    }

    for (j in 1 : N){
      index2 = which(as.character(mypheno[, 1]) == C1[j])
      DataPheno[j, ] = as.numeric(B1[index2, ])
    }

    DataGene1 = data.frame(DataGene, row.names = C1)
    names(DataGene1) = as.character(A[, 1])
    DataPheno1 = data.frame(DataPheno, row.names = C1)
    names(DataPheno1) = names(B1)
    snp_number = length(unique(A$SNP))
    SNPname = unique(A$SNP)

    #--------------- Model Selection Running ---------------
    #----------Processing genes with multiple SNPs----------------
    if (snp_number > 1){
      #---------Internal collenarity control per gene------------
      Q1 = cov(DataGene)
      E1 = eigen(Q1)
      E1value = E1$values
      E1score = E1$vectors
      E2 = which(E1value < 10^(-5))
      if (length(E2) == 0) SNPIndex = c(1 : length(SNPname))
      if (length(E2) > 0){
        keep = c()
        tempBefore = c()
        for (i in length(E2) : 1){
          temp = E1score[, E2[i]]
          tempNow = setdiff(which(abs(temp) > 10^(-5)), tempBefore)
          if (length(tempNow) > 0) keep = c(keep, tempNow[1])
          tempBefore = union(which(abs(temp) > 10^(-5)), tempBefore)
        }
        SNPIndex = c(setdiff(c(1 : length(E1value)), tempBefore), keep)
      }
      M1 = dim(DataGene)
      indexNA = which(is.na(colSums(DataPheno)) == TRUE)
      if (length(indexNA)==0){
        DataPheno2 = DataPheno
        DataPheno3 = DataPheno1
      }
      else{
        DataPheno2 = DataPheno[, -indexNA]
        DataPheno3 = DataPheno1[, -indexNA]
      }
      names3 = names(DataPheno3)
      M2 = dim(DataPheno2)

      SNPname = paste0("SNP", 1 : M1[2])
      DataGene2 = data.frame(DataGene[, SNPIndex])
      names(DataGene2) = SNPname[SNPIndex]
      Covariate = c()
      fmla0 = as.formula(paste("cbind(", paste(SNPname[SNPIndex], collapse= ","), ") ~", init_pc))
      fit = lm(formula = fmla0, data = DataGene2)
      thresholdIn = selectIn
      thresholdOut = selectOut

      for (rep in 1 : R){
        candidate = setdiff(c(1 : M2[2]), Covariate)
        pv = c()

        #--- add in ---
        for (i in 1 : length(candidate)){
          tempIndex = c(Covariate, candidate[i])
          tempData = DataPheno2[, tempIndex]
          tempNames = names3[tempIndex]
          data1 = data.frame(tempData)
          names(data1) = tempNames
          data2 = cbind(DataGene2, data1)
          fmla1 = as.formula(paste("cbind(", paste(SNPname[SNPIndex], collapse= ","), ") ~ ", med_pc, paste(tempNames, collapse= "+")))
          fit1 = lm(formula = fmla1, data = data2)
          # jointly using multi-SNP determing each phenotype p value using F test
          out1 = anova(fit1)$"Pr(>F)"
          out2 = out1[-c(2,3,4,length(out1))]
          temp1 = length(out2)
          if (temp1 <= length(tempIndex)){ pv[i] = 1}
          if (temp1 > length(tempIndex)){ pv[i] = out2[temp1]}
        }
        # determining if the most significant phenotype pass the threshold
        index1 = order(pv)[1]
        if (pv[index1] < thresholdIn){
          tempIndex0 = c(Covariate, candidate[index1])
          tempData0 = DataPheno2[, tempIndex0]
          tempNames0 = names3[tempIndex0]
          data10 = data.frame(tempData0)
          names(data10) = tempNames0
          data20 = cbind(DataGene2, data10)
          fmla = as.formula(paste("cbind(", paste(SNPname[SNPIndex], collapse= ","), ") ~ ", med_pc, paste(tempNames0, collapse= "+")))
          fit = lm(formula = fmla, data = data20)
          Covariate = c(Covariate, candidate[index1])
          cat("iteration =", rep, ", Cov =", tempNames0, "\n")
        }
        if (pv[index1] >= thresholdIn){ cat("iteration =", rep, "No new add-in", "\n") }

        #--- leave out ---
        out01 = anova(fit)$"Pr(>F)"
        out02 = out01[-c(2,3,4,length(out01))]
        if (length(out02) == 1) pv0 = 0
        if (length(out02) >= 2) pv0 = out02[-1]
        index2 = order(pv0)[length(pv0)]
        pvmax = pv0[index2]
        if (pvmax > thresholdOut){
          Covariate = Covariate[-index2]
          tempData01 = DataPheno2[, Covariate]
          tempNames01 = names3[Covariate]
          data11 = data.frame(tempData01)
          names(data11) = tempNames01
          data21 = cbind(DataGene2, data11)
          fmla = as.formula(paste("cbind(", paste(SNPname[SNPIndex], collapse= ","), ") ~ ", med_pc, paste(tempNames01, collapse= "+")))
          fit = lm(formula = fmla, data = data21)
          cat("iteration =", rep, ", Cov =", tempNames0, "\n")
        }
        if (pvmax <= thresholdOut){ cat("iteration =", rep, "No leave-out", "\n")}
      }
      data0 = cbind(DataGene2, PC)
      fmla0 = as.formula(paste("cbind(", paste(SNPname[SNPIndex], collapse= ","), ") ~ ", init_pc))
      fit0 = lm(formula = fmla0, data = data0)
      anova.compare = anova(fit, fit0)
      pp = anova.compare$`Pr(>F)`[2]
      gvalue = c(gvalue,pp)
    }

    #----------Processing genes with single SNP----------------
    else{
      M1 = dim(DataGene)
      #--------Missing data cleaning --------------------------
      indexNA = which(is.na(colSums(DataPheno)) == TRUE)
      if (length(indexNA)==0){
        DataPheno2 = DataPheno
        DataPheno3 = DataPheno1
      }
      else{
        DataPheno2 = DataPheno[, -indexNA]
        DataPheno3 = DataPheno1[, -indexNA]
      }
      names3 = names(DataPheno3)
      M2 = dim(DataPheno2)

      SNP = DataGene[, 1]
      Covariate = c()
      thresholdIn = selectIn
      thresholdOut = selectOut

      for (rep in 1 : R){
        candidate = setdiff(c(1 : M2[2]), Covariate)
        pv = c()

        #--- add in ---
        for (i in 1 : length(candidate)){
          tempIndex = c(Covariate, candidate[i])
          tempData = DataPheno2[, tempIndex]
          tempNames = names3[tempIndex]
          data1 = data.frame(tempData)
          names(data1) = tempNames
          data2 = cbind(SNP, data1)
          fmla1 = as.formula(paste("SNP ~ ", med_pc, paste(tempNames, collapse= "+")))
          fit1 = lm(formula = fmla1, data = data2)
          out1 = summary(fit1)$coefficients
          temp1 = dim(out1)[1]
          if (temp1 <= length(tempIndex + 3)){ pv[i] = 1}
          if (temp1 > length(tempIndex + 3)){ pv[i] = out1[temp1, 4]}
        }
        index1 = order(pv)[1]
        if (pv[index1] < thresholdIn){
          tempIndex0 = c(Covariate, candidate[index1])
          tempData0 = DataPheno2[, tempIndex0]
          tempNames0 = names3[tempIndex0]
          data10 = data.frame(tempData0)
          names(data10) = tempNames0
          data20 = cbind(SNP, data10)
          fmla = as.formula(paste("SNP ~ ", med_pc, paste(tempNames0, collapse= "+")))
          fit = lm(formula = fmla, data = data20)
          Covariate = c(Covariate, candidate[index1])
          cat("iteration =", rep, ", Cov =", tempNames0, "\n")
        }
        if (pv[index1] >= thresholdIn){ cat("iteration =", rep, "No new add-in", "\n") }

        #--- leave out ---
        out0 = summary(fit)$coefficients
        out01 = out0[-c(1 : 4), ]
        if (dim(out0)[1] == 4) pv0 = 0
        if (dim(out0)[1] == 5) pv0 = out01[4]
        if (dim(out0)[1] > 5) pv0 = out01[, 4]
        index2 = order(pv0)[length(pv0)]
        pvmax = pv0[index2]
        if (pvmax > thresholdOut){
          Covariate = Covariate[-index2]
          tempData01 = DataPheno2[, Covariate]
          tempNames01 = names3[Covariate]
          data11 = data.frame(tempData01)
          names(data11) = tempNames01
          data21 = cbind(SNP, data11)
          fmla = as.formula(paste("SNP ~ ", med_pc, paste(tempNames01, collapse= "+")))
          fit = lm(formula = fmla, data = data21)
          cat("iteration =", rep, ", Cov =", tempNames01, "\n")
        }
        if (pvmax <= thresholdOut){ cat("iteration =", rep, "No leave-out", "\n") }
      }
      data0 = cbind(DataGene, PC)
      colnames(data0)[1] = "SNP"
      fmla0 = as.formula(paste("SNP~ ", init_pc))
      fit0 = lm(formula = fmla0, data = data0)
      anova.compare = anova(fit, fit0)
      pp = anova.compare$`Pr(>F)`[2]
      gvalue = c(gvalue,pp)
    }
    kk = anova(fit)
    tk = cbind(rownames(kk),kk$`Pr(>F)`)
    myl = length(rownames(kk))
    left = rep(g,myl)
    d = cbind(left,tk)
    write.table(d,file=gp,append=T,col.names = FALSE)
    gd = cbind(glist,gvalue)
  }
  write.table(gd,file=gv,append=T,col.names = FALSE)
}
