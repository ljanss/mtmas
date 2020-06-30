# Multi-Trait Marker Selection (MTMAS).
# - LDinfo: matrix with LD information (variance-covariances) between the marker genotypes, with
#   snpnames on rows
# - Ntrait: number of traits, should match the number of gwasfiles (one per trait)
# - gwasfiles: vector with names of gwas output files, the minimum information in the gwas summary files should
#   be columns with the matching names "SNP" (marker names matching the names in the LDinfo file),
#   "P.value" and "effect" (the 'beta' or marker effect). The effects should be standardised to the scale of the
#   analysed response variable by dividing by the standard deviation of the responses (y, trait values) used in the GWAS.
# - PvalCutoff: cutoff for taking markers from the GWAS outputs
# - markerPenalty: vector with range of markerPenalties to apply, it depends on number of traits and weights
#   what level of penatly is useful, this can require some testing.
# - popSize, maxiter: parameters for the genetic algorithm: population size (number of solutions) to use per
#   iteration, and the maximum number of iterations to run.

mtmas <- function(LDinfo, Ntrait, gwasfiles, traitWeights,
                   PvalCutoff, markerPenalty, popSize, maxiter) {

   # Read files with GWAS outputs and build a vector to indicate which genotypes to
   # keep from the geno input. Collect betas and p-values for all traits.
   # The GWAS files have SNPs in different order, all GWAS results are reordered
   # to the order in the LDinfo file by matching to 'snpnames'.
   snpnames = rownames(LDinfo)
   geno_keep = rep(FALSE,ncol(geno))
   betas = matrix(0,ncol(geno),Ntrait)
   pvalues = matrix(1,ncol(geno),Ntrait)
   for(trait in 1:Ntrait) {
      gwas_result = read.table(gwasfiles[trait],header=TRUE,sep=",")
      snp_row_select = gwas_result[,"P.value"] < PvalCutoff
      geno_keep[match(gwas_result[snp_row_select,"SNP"],snpnames)] = TRUE
      betas[match(gwas_result$SNP,snpnames),trait] = gwas_result$effect
      pvalues[match(gwas_result$SNP,snpnames),trait] = gwas_result$P.value
   }
   geno = geno[,geno_keep]
   betas = betas[geno_keep,]
   colnames(betas) = paste("beta",seq(1:Ntrait),sep="")
   pvalues = pvalues[geno_keep,]
   colnames(pvalues) = paste("pvalue",seq(1:Ntrait),sep="")
   snpnames = snpnames[geno_keep]
   Nsnp = ncol(geno)
   Nid = nrow(geno)
   cat(paste("There are",Nsnp,"markers selected for the Multi-Trait Marker Selection\n"))

   gtVar = var(geno,na.rm=TRUE,use="pair")
   if (sum(is.na(gtVar))>0) {   # this can happen if there are many NA in the markers,
      gtVar[is.na(gtVar)] = 0   # just by chance, some covariances cannot be computed
   }
   smallest_eval = min(eigen(gtVar)$values)
   if(smallest_eval < 0) {
      gtVar = gtVar + diag(-smallest_eval,nrow(gtVar),ncol(gtVar))
      cat(paste(-smallest_eval,"was added to the LD matrix to make it semi-positive definite\n"))
   }

   betas_squared = betas^2

   # Function to compute explained variance from selected markers per trait.
   # The 'x' coming in the function is a vector with 0/1 to indicate selected markers.
   markerExplVar <- function(x) {
      select = (x==1)
      if (sum(select)==0) return(0)
      if (sum(select)==1) return(gtVar[select,select]*betas_squared[select,])
      Vcond = diag(chol(gtVar[select,select],pivot=FALSE, tol=0.001))^2
      return(t(betas_squared[select,])%*%Vcond)
   }

   library(GA)

   garun = list()
   solutions = matrix(0,Nsnp,length(markerPenalty))
   
   for(penalty in 1:length(markerPenalty)) {
      gares = ga(type="binary",fitness=function(x){sum(traitWeights*markerExplVar(x))-markerPenalty[penalty]*sum(x)},
               nBits=Nsnp,popSize=popSize, maxiter=maxiter, monitor=FALSE, keepBest=TRUE)
      if (penalty==1)
      garun = c(garun,gares)
      solutions[,penalty] = gares@solution
      cat(paste("Optimization for marker penalty",markerPenalty[penalty],"selected",sum(gares@solution),"markers\n"))
   }
   
   colnames(solutions) = paste("solution",seq(1:length(markerPenalty)),sep="")
   NmarkerSelect = apply(solutions,2,sum)
   explainedVar = apply(solutions,2,markerExplVar)
   markerInfo = data.frame(snpnames, pvalues, betas, solutions)
   out = list(MarkerInfo=markerInfo, MarkerPenalty=markerPenalty, NmarkerSelect=NmarkerSelect, ExplainedVar=explainedVar,
           gtVar=gtVar, GAruns=garun)
   # the ga outputs from every run are now not included in the return output

   return(out)
   
}
