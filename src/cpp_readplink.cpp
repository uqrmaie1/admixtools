/*
 *
 * This is based on Gad Abraham's "plink2R" package:
 * https://github.com/gabraham/plink2R/blob/master/plink2R/
 *
 */

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <sstream>
#include <string>

using namespace Rcpp;


#define PACK_DENSITY 4
#define PLINK_NA 3

/* 3 is 11 in binary, we need a 2 bit mask for each of the 4 positions */
#define MASK0 3	  /* 3 << 2 * 0 */
#define MASK1 12  /* 3 << 2 * 1 */
#define MASK2 48  /* 3 << 2 * 2 */
#define MASK3 192 /* 3 << 2 * 3 */


/*
 *                   plink BED           sparsnp
 * minor homozyous:  00 => numeric 0     10 => numeric 2
 * heterozygous:     10 => numeric 2     01 => numeric 1
 * major homozygous: 11 => numeric 3     00 => numeric 0
 * missing:          01 => numeric 1     11 => numeric 3
 *
 *
 * http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml says,
 * The bytes in plink are read backwards HGFEDCBA, not GHEFCDAB, but we read
 * them forwards as a character (a proper byte)
 *
 * By default, plink usage dosage of the *major* allele, since allele A1 is
 * usually the minor allele and the code "1" refers to the second allele A2,
 * so that "11" is A2/A2 or major/major.
 *
 * We always use minor allele dosage, to be consistent with the output from
 * plink --recodeA which used minor allele dosage by default.
 *
 * out: array of genotypes
 * in: array of packed genotypes (bytes)
 * n: number of bytes in input
 *
 */
void decode_plink(unsigned char *out,
                  const unsigned char *in, const unsigned int n) {
  unsigned int i, k;
  unsigned char tmp, geno;
  unsigned int a1, a2;

  for(i = 0 ; i < n ; ++i) {
    tmp = in[i];
    k = PACK_DENSITY * i;

    /* geno is interpreted as a char, however a1 and a2 are bits for allele 1 and
     * allele 2. The final genotype is the sum of the alleles, except for 01
     * which denotes missing.
     */
    geno = (tmp & MASK0);
    a1 = !(geno & 1);
    a2 = !(geno >> 1);
    out[k] = (geno == 1) ? 3 : a1 + a2;
    k++;

    geno = (tmp & MASK1) >> 2;
    a1 = !(geno & 1);
    a2 = !(geno >> 1);
    out[k] = (geno == 1) ? 3 : a1 + a2;
    k++;

    geno = (tmp & MASK2) >> 4;
    a1 = !(geno & 1);
    a2 = !(geno >> 1);
    out[k] = (geno == 1) ? 3 : a1 + a2;
    k++;

    geno = (tmp & MASK3) >> 6;
    a1 = !(geno & 1);
    a2 = !(geno >> 1);
    out[k] = (geno == 1) ? 3 : a1 + a2;
  }
}


// [[Rcpp::export]]
NumericMatrix cpp_read_plink(String bedfile, int nsnp, int nind, IntegerVector indvec,
                             int first, int last, bool transpose = false, bool verbose = true) {

  int val;
  long long len, bytespersnp;
  int readsnps = last - first;
  int headersize = 3;

  std::ifstream in(bedfile.get_cstring(), std::ios::in | std::ios::binary);

  if(!in) {
    stop("Genotype file not found!");
  }
  in.seekg(0, std::ifstream::end);
  // file size in bytes
  len = (long long)in.tellg();
  bytespersnp = (long long)ceil((double)nind / PACK_DENSITY);
  if((len-headersize) != nsnp*bytespersnp) {
    stop("Unexpected bed file size!");
  }

  int nindused = 0;
  int* blockused = new int[bytespersnp];
  for(int i = 0 ; i < bytespersnp; i++) {
    blockused[i] = 0;
  }
  for(int i = 0 ; i < nind; i++) {
    if(indvec[i] == 1) {
      nindused++;
      blockused[i/PACK_DENSITY] = 1;
    }
  }

  NumericMatrix geno(transpose?nindused:readsnps, transpose?readsnps:nindused);
  std::fill(geno.begin(), geno.end(), NA_REAL);

  in.seekg(3+first*bytespersnp, std::ifstream::beg);
  unsigned char* tmp = new unsigned char[bytespersnp + 1];
  tmp[bytespersnp] = '\0';
  unsigned char tmpi;

  // Allocate more than the sample size since data must take up whole bytes
  unsigned char* tmp2 = new unsigned char[bytespersnp * PACK_DENSITY + 1];
  tmp2[bytespersnp * PACK_DENSITY] = '\0';

  int k, tmp3;
  for(int j = 0 ; j < readsnps; j++) {
    //for(unsigned int j = 0 ; j < 3; j++) {
    if(verbose && j % 1000 == 0) Rcout << "\r" << j/1000 << "k SNPs read...";

    // read raw genotypes
    in.read((char*)tmp, sizeof(char) * bytespersnp);

    for(int l = 0; l < bytespersnp; l++) {
      if(!blockused[l]) continue;

      tmpi = tmp[l];
      k = PACK_DENSITY * l;

      /* geno is interpreted as a char, however a1 and a2 are bits for allele 1 and
       * allele 2. The final genotype is the sum of the alleles, except for 11
       * which denotes missing.
       */
      tmp3 = (tmpi & MASK0);      tmp2[k] = tmp3 < 2 ? tmp3 + 2 : 3 - tmp3;
      tmp3 = (tmpi & MASK1) >> 2; tmp2[k+1] = tmp3 < 2 ? tmp3 + 2 : 3 - tmp3;
      tmp3 = (tmpi & MASK2) >> 4; tmp2[k+2] = tmp3 < 2 ? tmp3 + 2 : 3 - tmp3;
      tmp3 = (tmpi & MASK3) >> 6; tmp2[k+3] = tmp3 < 2 ? tmp3 + 2 : 3 - tmp3;
    }

    int c = 0;
    if(!transpose) {
      for(int i = 0; i < nind; i++) {
        if(!indvec[i]) continue;
        val = (double)tmp2[i];
        if(val != 3) geno(j, c) = val;
        c++;
      }
    } else {
      for(int i = 0; i < nind; i++) {
        if(!indvec[i]) continue;
        val = (double)tmp2[i];
        if(val != 3) geno(c, j) = val;
        c++;
      }
    }
  }
  if(verbose) Rcout << std::endl;

  delete[] tmp;
  delete[] tmp2;
  delete[] blockused;
  in.close();

  return geno;
}


// [[Rcpp::export]]
NumericVector cpp_plink_ploidy(String genofile, int nsnp, int nind, IntegerVector indvec, int ntest = 1000) {

  int val, k, tmp3;
  long long bytespersnp, len;
  int headersize = 3;
  ntest = std::min(ntest, nsnp);
  std::ifstream in(genofile.get_cstring(), std::ios::in | std::ios::binary);
  if(!in) {
    stop("Genotype file not found!");
  }
  in.seekg(0, std::ifstream::end);
  len = (long long)in.tellg();
  bytespersnp = (long long)ceil((double)nind / PACK_DENSITY);
  if((len-headersize) != nsnp*bytespersnp) {
    stop("Unexpected bed file size!");
  }

  NumericVector ploidy(nind);
  ploidy.fill(1.0);
  int* blockused = new int[bytespersnp];
  for(int i = 0 ; i < bytespersnp; i++) {
    blockused[i] = 0;
  }
  for(int i = 0 ; i < nind; i++) {
    if(indvec[i] != -1) {
      blockused[i/PACK_DENSITY] = 1;
    }
  }
  unsigned char* tmp = new unsigned char[bytespersnp];
  unsigned char tmpi;
  unsigned char* tmp2 = new unsigned char[bytespersnp * PACK_DENSITY];

  in.seekg(headersize, std::ifstream::beg);
  for(int j = 0 ; j < ntest; j++) {
    in.read((char*)tmp, sizeof(char) * bytespersnp);

    for(int l = 0 ; l < bytespersnp; ++l) {
      if(!blockused[l]) continue;
      tmpi = tmp[l];
      k = PACK_DENSITY * l;
      tmp3 = (tmpi & MASK0);      tmp2[k] = tmp3 < 2 ? tmp3 + 2 : 3 - tmp3;
      tmp3 = (tmpi & MASK1) >> 2; tmp2[k+1] = tmp3 < 2 ? tmp3 + 2 : 3 - tmp3;
      tmp3 = (tmpi & MASK2) >> 4; tmp2[k+2] = tmp3 < 2 ? tmp3 + 2 : 3 - tmp3;
      tmp3 = (tmpi & MASK3) >> 6; tmp2[k+3] = tmp3 < 2 ? tmp3 + 2 : 3 - tmp3;
    }
    for(int i = 0; i < nind; i++) {
      if(indvec[i] == -1) continue;
      val = (double)tmp2[i];
      if(val == 1) ploidy(i) = 2.0;
    }
  }

  delete[] tmp;
  delete[] tmp2;
  delete[] blockused;
  in.close();

  return ploidy;
}

// [[Rcpp::export]]
List cpp_plink_to_afs(String genofile, int nsnp, int nind, IntegerVector indvec,
                          int first, int last, IntegerVector ploidy,
                          bool transpose, bool verbose) {

  int val, k, tmp3, pop;

  long long len, bytespersnp;
  int headersize = 3;
  int readsnps = std::min(nsnp, last?last:nsnp) - first;

  std::ifstream in(genofile.get_cstring(), std::ios::in | std::ios::binary);

  if(!in) {
    Rcout << "Error reading file " << genofile.get_cstring() << std::endl;
    throw std::runtime_error("io error");
  }
  in.seekg(0, std::ifstream::end);
  // file size in bytes
  len = (long long)in.tellg();
  bytespersnp = (long long)ceil((double)nind / PACK_DENSITY);
  if((len-headersize) != nsnp*bytespersnp) {
    stop("Unexpected bed file size!");
  }

  int numpop = max(indvec)+1;

  int* blockused = new int[bytespersnp];
  for(int i = 0 ; i < bytespersnp; i++) {
    blockused[i] = 0;
  }
  for(int i = 0 ; i < nind; i++) {
    if(indvec[i] != -1) {
      blockused[i/PACK_DENSITY] = 1;
    }
  }

  //set ploidy
  unsigned char* tmp = new unsigned char[bytespersnp];
  unsigned char tmpi;
  // Allocate more than the sample size since data must take up whole bytes
  unsigned char* tmp2 = new unsigned char[bytespersnp * PACK_DENSITY];

  NumericMatrix afs(transpose?numpop:readsnps, transpose?readsnps:numpop);
  std::fill(afs.begin(), afs.end(), NA_REAL);
  NumericMatrix counts(transpose?numpop:readsnps, transpose?readsnps:numpop);
  std::fill(counts.begin(), counts.end(), 0.0);

  in.seekg(headersize+first*bytespersnp, std::ifstream::beg);

  NumericVector altalleles(numpop);
  NumericVector observedalleles(numpop);

  for(int j = 0 ; j < readsnps; j++) {
    if(verbose && j % 1000 == 0) Rcout << "\r" << j/1000 << "k SNPs read...";

    // read raw genotypes
    in.read((char*)tmp, sizeof(char) * bytespersnp);

    for(int l = 0 ; l < bytespersnp; l++) {
      if(!blockused[l]) continue;

      tmpi = tmp[l];
      k = PACK_DENSITY * l;

      /* geno is interpreted as a char, however a1 and a2 are bits for allele 1 and
       * allele 2. The final genotype is the sum of the alleles, except for 11
       * which denotes missing.
       */
      tmp3 = (tmpi & MASK0);      tmp2[k] = tmp3 < 2 ? tmp3 + 2 : 3 - tmp3;
      tmp3 = (tmpi & MASK1) >> 2; tmp2[k+1] = tmp3 < 2 ? tmp3 + 2 : 3 - tmp3;
      tmp3 = (tmpi & MASK2) >> 4; tmp2[k+2] = tmp3 < 2 ? tmp3 + 2 : 3 - tmp3;
      tmp3 = (tmpi & MASK3) >> 6; tmp2[k+3] = tmp3 < 2 ? tmp3 + 2 : 3 - tmp3;
    }

    altalleles.fill(0);
    observedalleles.fill(0);
    for(int i = 0; i < nind; i++) {
      if(indvec[i] == -1) continue;
      val = (double)tmp2[i];
      if(val != 3) {
        pop = indvec(i);
        altalleles(pop) += val / (3.0-ploidy(i));
        observedalleles(pop) += ploidy(i);
      }
    }

    if(!transpose) {
      for(int i = 0; i < numpop; i++) {
        if(observedalleles(i) > 0) {
          afs(j, i) = altalleles(i)/observedalleles(i);
          counts(j, i) = observedalleles(i);
        }
      }
    } else {
      for(int i = 0; i < numpop; i++) {
        if(observedalleles(i) > 0) {
          afs(i, j) = altalleles(i)/observedalleles(i);
          counts(i, j) = observedalleles(i);
        }
      }
    }
  }
  if(verbose) Rcout << std::endl;

  delete[] tmp;
  delete[] tmp2;
  delete[] blockused;
  in.close();

  return Rcpp::List::create(_["afs"] = afs, _["counts"] = counts);
}

