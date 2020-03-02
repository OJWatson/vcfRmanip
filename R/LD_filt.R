
##--------------------------------------------------------------------
## LINKAGE DISEQUILIBRIUM
##--------------------------------------------------------------------

#' @title cormat
#' @description Calculates genetic autocorrelation among variants.
#' @details Expects a matrix of allele frequencies with rows as loci and columns as samples
#'
#' @author Bob Verity & Nick Brazeau
#'
#' @export
#'

# Calculate pairwise correlation matrix between observations. Input is a matrix with n rows corresponding to n multivariate measurements, output is a n-by-n correlation matrix. NA values are imputed as the mean.

corMat <- function(m) {
  tol <- 1e-9 # tolerance for denominator if AF become the same
  n <- nrow(m) # number of loci
  c <- matrix(NA,n,n)
  diag(c) <- 1

  # precompute these
  xi_min_mui <- vector("list", n)
  sum_squared_xi_min_mui <- vector("list", n)
  for (i in 1:n) {

    x1 <- unlist(m[i,])
    mu1 <- mean(x1,na.rm=TRUE)
    xi_min_mui[[i]] <- x1-mu1
    sum_squared_xi_min_mui[[i]] <- sum((x1-mu1)^2,na.rm=TRUE)

  }

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {

      c[i,j] <- c[j,i] <- sum(xi_min_mui[[i]]*xi_min_mui[[j]],na.rm=TRUE)/sqrt( sum_squared_xi_min_mui[[i]]*sum_squared_xi_min_mui[[j]] + tol)

    }
  }
  return(c)
}



#' @title genautocorr setup
#' @description Calculates genetic autocorrelation for later linkage disequilibrium filtering.
#' @details From an object of class \code{vcfR}, calculate the genetic autocorrelation as the estimate of linkage disequilibrium by genetic distance.
#' @param vcffile A variant call file (vcf) path. This VCF will be converted to an object of class \code{vcfR}.
#'
#' @author Bob Verity & Nick Brazeau
#'
#' @export
#'

genautocorr <- function(vcffile = NULL, vcfR = NULL, biallelicsnps=TRUE){


  # -----------------------------------------------------
  # Read and check input
  #------------------------------------------------------
  if(is.null(vcffile)){
    if(class(vcfR) != "vcfR"){
      stop("vcfR object must be of class vcfR")
    }
    vcf <- vcfR
  } else{
    vcf <- vcfR::read.vcfR(file=vcffile, verbose=T) # read vcf
  }
  if(biallelicsnps == FALSE){
    stop("Must take in biallelic SNPs")
  }
  vcf <-vcfR::extract.indels(vcf, return.indels = F) # subset to SNPs
  vcf <- vcf[vcfR::is.biallelic(vcf)] # subset to biallelic

  # extract the within sample allele frequencies
  ad <- vcfR::extract.gt(vcf, element = "AD")
  refad <- masplit(ad, record=1, sort=0, decreasing = 0)
  altad <- masplit(ad, record=2, sort=0, decreasing = 0)

  NRAF <- altad/(altad + refad)
  NRAF[is.nan(NRAF)] <- NA # the 0,0 are returning Nans


  # extract distances via positions from vcfR object
  CHROM <- vcfR::getCHROM(vcf)
  POS <- vcfR::getPOS(vcf)
  vcfpos <- cbind.data.frame(CHROM, POS)
  vcfdf <- cbind.data.frame(vcfpos, NRAF)

  vcflist <- split(vcfdf, vcfdf$CHROM)

  # -----------------------------------------------------
  # calculations
  #------------------------------------------------------

  cormatgendistwrapper <- function(vcfdf_fromlist){

    # get correlation matrix. NA values are imputed as the mean
    df1 <- vcfdf_fromlist[, !colnames(vcfdf_fromlist) %in% c("CHROM", "POS")]
    c <- corMat(df1)
    # get distance between SNPs. This can be extracted from the row names of the vcf
    df2 <- vcfdf_fromlist[, colnames(vcfdf_fromlist) %in% c("POS")]
    gendist <- as.matrix(dist(df2))

    ret <- list(vcfAF = vcfdf_fromlist, corMat=c, gendist=gendist)
    return(ret)
  }


  retlist <- parallel::mclapply(vcflist, cormatgendistwrapper)

  return(retlist)

}



#' @title vcfR2LDfiltered
#' @description Filtering an object of class \code{vcfR} for linkage
#'   disequilibrium via genetic autocorrelation.
#' @param vcffile A variant call file (vcf) path. This VCF will be converted
#'   to an object of class \code{vcfR}.
#' @param vcfR A variant call file objet (vcf) path.
#' @param genautocorrresult Genetic autocorrelation object produced
#'   by \code{genautocorr}
#' @param threshDist Genetic distance to filter at.
#' @param threshR2 Genetic autocorrelation to filter at
#' @param random Are sites being filtered chosen at random from sites that need
#'   to be filtered. Default = FALSE
#' @param bial
#'
#'
#' @export
#'

vcfR2LDfiltered <- function(vcffile = NULL, vcfR = NULL,
                            genautocorrresult=NULL,
                            threshDist=NULL,
                            threshR2 = NULL,
                            random = FALSE){


  # -----------------------------------------------------
  # Read and check input
  #------------------------------------------------------
  if(is.null(threshDist) && is.null(threshR2)){
    stop("Must specify one of linkage disequilibrium threshold (see help and tutorial).")
  }
  if(!is.null(threshDist) && !is.null(threshR2)){
    message("Both threshDist and threshR2 are specified. threshDist will be used")
    threshR2 <- NULL
  }
  if(is.null(genautocorrresult)){
    stop("Must specify a genetic autocorrelation results object using the genautocorr function")
  }

  if(is.null(vcffile)){
    if(class(vcfR) != "vcfR"){
      stop("vcfR object must be of class vcfR")
    }
    vcf <- vcfR
  } else{
    vcf <- vcfR::read.vcfR(file=vcffile, verbose=T) # read vcf
  }

  vcf <-vcfR::extract.indels(vcf, return.indels = F) # subset to SNPs
  vcf <- vcf[vcfR::is.biallelic(vcf)] # subset to biallelic

  vcfdf <- cbind.data.frame(vcf@fix, vcf@gt)
  vcflist <- split(vcfdf, f=factor(vcfdf$CHROM))
  if(length(vcflist) != length(genautocorrresult)){
    stop("The number of chromosomes in the vcfR objec and the results from the genetic autocorrelation analysis differ.")
  }
  # -----------------------------------------------------
  # Filter based on distance
  #------------------------------------------------------

  filter_distance <- function(vcflist, genautocorrresult){
    gendist <- as.matrix(genautocorrresult$gendist)
    diag(gendist) <- Inf	# block self-comparison
    while (any(gendist<threshDist)) {
      w <- which(gendist<threshDist, arr.ind=TRUE)
      if(random) {
        pos <- sample(length(w),1,FALSE)
      } else {
        pos <- 1
      }
      vcflist <- vcflist[-w[pos, pos],]
      gendist <- gendist[-w[pos, pos],-w[pos, pos]]

    }
    return(vcflist)
  }


  # -----------------------------------------------------
  # Filter based on autocorr
  #------------------------------------------------------

  filter_autocorr <- function(vcflist, genautocorrresult){
    corMat <- as.matrix(genautocorrresult$corMat)
    diag(corMat) <- Inf	# block self-comparison
    while (any(corMat>threshR2)) {
      w <- which(corMat>threshR2, arr.ind=TRUE)
      if(random) {
        pos <- sample(length(w),1,FALSE)
      } else {
        pos <- 1
      }
      vcflist <- vcflist[-w[pos, pos],]
      corMat <- corMat[-w[pos, pos],-w[pos, pos]]

    }
    return(vcflist)
  }

  if(is.null(threshR2)) {
    updatedvcflist <- parallel::mcmapply(filter_distance, vcflist, genautocorrresult, SIMPLIFY = F)
    fm <- paste("##Filtered for genetic autocorrelation with a threshold distance of", threshDist)
  } else {
    updatedvcflist <- parallel::mcmapply(filter_autocorr, vcflist, genautocorrresult, SIMPLIFY = F)
    fm <- paste("##Filtered for genetic autocorrelation with a maximum r2 of", threshR2)
  }
  updatedvcfdf <- do.call(rbind, updatedvcflist)

  # -----------------------------------------------------
  # Return to vcfR object
  #------------------------------------------------------
  fix <- as.matrix(updatedvcfdf[,1:8])
  gt <- as.matrix(updatedvcfdf[,9:ncol(updatedvcfdf)])
  meta <- append(vcf@meta, fm)

  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)

  return(newvcfR)

}



