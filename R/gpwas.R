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
gpwas = function(ingeno,inpheno,inpc,gp,gv,R=num,pc=3,g=NULL,selectIn=0.01,selectOut=0.01){
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)
  mygeno = read.table(ingeno,sep='\t',head=T, check.names = FALSE) # 'mygeno' is the SNP data for all the genes
  if (missing(g)){
    mygeno = mygeno
  } else {
    gfile = read.table(g,header=T)
    gset = as.character(gfile$Candidate)
    mygeno = mygeno[mygeno$Gene %in% gset, ]
  }
  mypheno = read.table(inpheno, head = TRUE) # 'mypheno' is the phenotype data
  gene = unique(mygeno$Gene) # 'gene' is the names for all the genes
  glist = c()
  gvalue = c()

  plist = c()
  for (i in c(1:pc)){
    plist = c(plist,paste('PC[, ',toString(i),']',sep=''))
  }
  init_pc = paste(plist,collapse= " + ") # equation expression for PC scores in the model
  med_pc = paste(init_pc,' + ',sep='')

  for (g in gene){ # iterating all candidate genes
    #--------------- Data pre-processing for the gth gene ---------------
    print (g)
    glist = c(glist,g) # the vector of the gene names
    tp = mygeno[mygeno$Gene==g,] # SNPs for the gene in the gth iteration
    th = tp$SNP[1]
    chr = unlist(strsplit(as.character(th),'_'))[1]
    chrom = gsub('S','',chr)
    PC = read.table(paste(inpc,'/exclude-chr',chrom,'.txt',sep=''),sep=' ',head=T)
    PC = PC[,c(1:pc)]

    A = tp[,-1] # delete the first column in the dataset 'tp'
    B1 = mypheno[, -1] # delete the first column in the dataset 'mypheno'
    M = dim(B1)
    C = names(A)
    C1 = intersect(C, as.character(mypheno[, 1])) # obtain the common genotypes shared by genotype data 'tp' and the phenotype data 'mypheno'
    N = length(C1)
    DataGene = matrix(0, N, dim(A)[1])
    DataPheno = matrix(0, N, M[2])

    #-------------- Prepare the genotype data 'DataGene' and the phenotype data 'DataPheno' with the common genotypes in 'C1' for the gth gene -------------
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

    #--------------- Model Selection and Inference ---------------
    #---------- Genes with multiple SNPs----------------
    if (snp_number > 1){
      #---------Internal collenarity control for the SNPs of each gene------------
      Q1 = cov(DataGene)
      E1 = eigen(Q1) # eigen decomposition for the covariance of the SNPs data
      E1value = E1$values
      E1score = E1$vectors
      E2 = which(E1value < 10^(-5)) # find the very small eigenvalues which indicate collenarity
      if (length(E2) == 0) SNPIndex = c(1 : length(SNPname))
      if (length(E2) > 0){
        keep = c()
        tempBefore = c()
        for (i in length(E2) : 1){ # start with the smallest eigenvalue smaller than 10^(-5)
          temp = E1score[, E2[i]] # get the corresponding eigenvector
          tempNow = setdiff(which(abs(temp) > 10^(-5)), tempBefore)
          if (length(tempNow) > 0) keep = c(keep, tempNow[1]) # only keep the first component of the eigenvectors (corresponding to those very small eigenvalues)
          # with absolute value larger than 10^(-5)
          tempBefore = union(which(abs(temp) > 10^(-5)), tempBefore)
        }
        SNPIndex = c(setdiff(c(1 : length(E1value)), tempBefore), keep) # 'SNPIndex' contains the final SNPs going into the step-wise selection model
      }
      #--------- delete the missing values ---------
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
      #--------- Initial model with SNPs as responses and PC scores as covariates ---------
      SNPname = paste0("SNP", 1 : M1[2])
      DataGene2 = data.frame(DataGene[, SNPIndex]) # SNP data with the SNPs in the selected 'SNPIndex' set
      names(DataGene2) = SNPname[SNPIndex]
      fmla0 = as.formula(paste("cbind(", paste(SNPname[SNPIndex], collapse= ","), ") ~", init_pc))
      fit = lm(formula = fmla0, data = DataGene2)

      #--------- Step-wise model selection ---------
      Covariate = c()
      thresholdIn = selectIn # threshold for selection of phenotypes in stepwise selection procedure
      thresholdOut = selectOut # threshold for drop-out of phenotypes in stepwise selection procedure

      for (rep in 1 : R){ # iterating the set number of iterations R
        candidate = setdiff(c(1 : M2[2]), Covariate) # candidate phenotype covariates in the (rep)th repetition
        pv = c()
        #--- add in ---
        for (i in 1 : length(candidate)){# Add each one of the candidate covariates into the current model; calculate its pvalue
          tempIndex = c(Covariate, candidate[i])
          tempData = DataPheno2[, tempIndex]
          tempNames = names3[tempIndex]
          data1 = data.frame(tempData)
          names(data1) = tempNames
          data2 = cbind(DataGene2, data1)
          fmla1 = as.formula(paste("cbind(", paste(SNPname[SNPIndex], collapse= ","), ") ~ ", med_pc, paste(tempNames, collapse= "+"))) # add one covariate 'candidate[i]'
          fit1 = lm(formula = fmla1, data = data2) # fit the model with 'candidate[i]'
          out1 = anova(fit1)$"Pr(>F)"
          out2 = out1[-c(2:(pc+1),length(out1))]
          temp1 = length(out2)
          if (temp1 <= length(tempIndex)){ pv[i] = 1}
          if (temp1 > length(tempIndex)){ pv[i] = out2[temp1]} # 'pv[i]' is the pvalue associated with 'candidate[i]'
        }
        index1 = order(pv)[1] # order the pvalues for the candidate covariates
        if (pv[index1] < thresholdIn){ # if the smallest candidate pvalue is less than 'thresholdIn', add it into the current model
          tempIndex0 = c(Covariate, candidate[index1])
          tempData0 = DataPheno2[, tempIndex0]
          tempNames0 = names3[tempIndex0]
          data10 = data.frame(tempData0)
          names(data10) = tempNames0
          data20 = cbind(DataGene2, data10)
          fmla = as.formula(paste("cbind(", paste(SNPname[SNPIndex], collapse= ","), ") ~ ", med_pc, paste(tempNames0, collapse= "+")))
          fit = lm(formula = fmla, data = data20) # update the model
          Covariate = c(Covariate, candidate[index1]) # update the covariates in the model
          cat("iteration =", rep, ", Cov =", tempNames0, "\n")
        }
        if (pv[index1] >= thresholdIn){ cat("iteration =", rep, "No new add-in", "\n") } # if the smallest candidate pvalue is larger than 'thresholdIn',
        # stop adding new varaible into the model

        #--- leave out ---
        out01 = anova(fit)$"Pr(>F)" # pvalues for all the covariates currently in the model
        out02 = out01[-c(2:(pc+1),length(out01))] # delete the pvalues corresponding to the PC scores
        if (length(out02) == 1) pv0 = 0
        if (length(out02) >= 2) pv0 = out02[-1]
        index2 = order(pv0)[length(pv0)]
        pvmax = pv0[index2] # get the maximum pvalue
        if (pvmax > thresholdOut){ # if the maximum pvalue is larger than 'thresholdOut', delete the least significant variable in the model
          Covariate = Covariate[-index2] # delete the least significant variable in the model
          tempData01 = DataPheno2[, Covariate]
          tempNames01 = names3[Covariate]
          data11 = data.frame(tempData01)
          names(data11) = tempNames01
          data21 = cbind(DataGene2, data11)
          fmla = as.formula(paste("cbind(", paste(SNPname[SNPIndex], collapse= ","), ") ~ ", med_pc, paste(tempNames01, collapse= "+")))
          fit = lm(formula = fmla, data = data21) # re-fit the model
          cat("iteration =", rep, ", Cov =", tempNames0, "\n")
        }
        if (pvmax <= thresholdOut){ cat("iteration =", rep, "No leave-out", "\n")} # if the maximum pvalue is less than 'thresholdOut',
        # keep all the current variables
      }
      data0 = cbind(DataGene2, PC)
      fmla0 = as.formula(paste("cbind(", paste(SNPname[SNPIndex], collapse= ","), ") ~ ", init_pc))
      fit0 = lm(formula = fmla0, data = data0) # compare the final model from step-wise selection and the initial model with only the PC scores
      anova.compare = anova(fit, fit0) # Anova analysis for this comparison
      pp = anova.compare$`Pr(>F)`[2] # pvalue for testing whether some of the phenotype variables in the final model are nonzero
      gvalue = c(gvalue,pp) # gather such pvalues for each gene
    }

    #---------- Genes with single SNP (the codes are very similar to the multiple SNPs case, we ignore some of the comments on the codes) ----------
    else{
      M1 = dim(DataGene)
      #---------- delete the missing values ----------
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
      #--------- Initial model with the SNP as response and PC scores as covariates ---------
      SNP = DataGene[, 1]
      single_SNP = as.formula(paste('SNP', paste(1, init_pc, sep=' + '), sep=' ~ '))
      fit = lm(single_SNP)

      #--------- Step-wise model selection ---------
      Covariate = c()
      thresholdIn = selectIn # threshold for selection of phenotypes in stepwise selection procedure
      thresholdOut = selectOut # threshold for drop-out of phenotypes in stepwise selection procedure

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
          if (temp1 <= length(tempIndex + pc)){ pv[i] = 1}
          if (temp1 > length(tempIndex + pc)){ pv[i] = out1[temp1, 4]}
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
        out01 = out0[-c(1 : (pc+1)), ]
        if (dim(out0)[1] == (pc+1)) pv0 = 0
        if (dim(out0)[1] == (pc+2)) pv0 = out01[4]
        if (dim(out0)[1] > (pc+2)) pv0 = out01[, 4]
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
    write.table(d,file=gp,append=T,col.names = FALSE) # generating the file about each inspected gene with the selected phenotypes
    gd = cbind(glist,gvalue) # the data frame with the names of genes and their associated pvalues 'pp'
  }
  write.table(gd,file=gv,append=T,col.names = FALSE) # generating the file containing the p-values for each gene by GPWAS
}
