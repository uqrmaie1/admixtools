// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <sstream>
#include <string>

using namespace Rcpp;

#define PACK_DENSITY 4

/* 3 is 11 in binary, we need a 2 bit mask for each of the 4 positions */
#define MASK0 3	  /* 3 << 2 * 0 */
#define MASK1 12  /* 3 << 2 * 1 */
#define MASK2 48  /* 3 << 2 * 2 */
#define MASK3 192 /* 3 << 2 * 3 */


// [[Rcpp::export]]
NumericMatrix cpp_read_packedancestrymap(String genofile, int nsnp, int nind, IntegerVector indvec,
                                         int first, int last, bool transpose = false, bool verbose = true) {

  int val;
  long long len, bytespersnp;
  int readsnps = last - first;

  std::ifstream in(genofile.get_cstring(), std::ios::in | std::ios::binary);

  if(!in) {
    Rcout << "Error reading file " << genofile.get_cstring() << std::endl;
    throw std::runtime_error("io error");
  }
  in.seekg(0, std::ifstream::end);
  // file size in bytes
  len = (long long)in.tellg();
  bytespersnp = len/(nsnp+1);

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

  // char* header = new char[bytespersnp];
  // in.seekg(0, std::ifstream::beg);
  // in.read((char*)header, bytespersnp);
  // Rcout << "header " << header << std::endl;

  in.seekg((first+1)*bytespersnp, std::ifstream::beg);
  char* tmp = new char[bytespersnp + 1];
  tmp[bytespersnp] = '\0';
  char tmpi;

  // Allocate more than the sample size since data must take up whole bytes
  char* tmp2 = new char[bytespersnp * PACK_DENSITY + 1];
  tmp2[bytespersnp * PACK_DENSITY] = '\0';

  int k;
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
      tmp2[k] = (tmpi & MASK3) >> 6;
      tmp2[k+1] = (tmpi & MASK2) >> 4;
      tmp2[k+2] = (tmpi & MASK1) >> 2;
      tmp2[k+3] = (tmpi & MASK0);
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
NumericVector cpp_packedancestrymap_ploidy(String genofile, int nsnp, int nind, IntegerVector indvec, int ntest = 1000) {

  int val, k;
  long long bytespersnp, len;

  ntest = std::min(ntest, nsnp);
  std::ifstream in(genofile.get_cstring(), std::ios::in | std::ios::binary);
  if(!in) {
    Rcout << "Error reading file " << genofile.get_cstring() << std::endl;
    throw std::runtime_error("io error");
  }
  in.seekg(0, std::ifstream::end);
  len = (long long)in.tellg();
  bytespersnp = len/(nsnp+1);

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

  in.seekg(bytespersnp, std::ifstream::beg);
  for(int j = 0 ; j < ntest; j++) {
    in.read((char*)tmp, sizeof(char) * bytespersnp);

    for(int l = 0 ; l < bytespersnp; ++l) {
      if(!blockused[l]) continue;
      tmpi = tmp[l];
      k = PACK_DENSITY * l;
      tmp2[k] = (tmpi & MASK3) >> 6;
      tmp2[k+1] = (tmpi & MASK2) >> 4;
      tmp2[k+2] = (tmpi & MASK1) >> 2;
      tmp2[k+3] = (tmpi & MASK0);
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
List cpp_packedancestrymap_to_afs(String genofile, int nsnp, int nind, IntegerVector indvec,
                                      int first, int last, IntegerVector ploidy,
                                      bool transpose, bool verbose) {
  //same arguments as cpp_read_packedancestrymap, except indvec assigns populations. indiv not used: -1

  int val, k, pop;

  long long len, bytespersnp;
  int readsnps = std::min(nsnp, last?last:nsnp) - first;

  std::ifstream in(genofile.get_cstring(), std::ios::in | std::ios::binary);

  if(!in) {
    Rcout << "Error reading file " << genofile.get_cstring() << std::endl;
    throw std::runtime_error("io error");
  }
  in.seekg(0, std::ifstream::end);
  // file size in bytes
  len = (long long)in.tellg();
  bytespersnp = len/(nsnp+1);

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

  in.seekg((first+1)*bytespersnp, std::ifstream::beg);

  NumericVector altalleles(numpop);
  NumericVector observedalleles(numpop);

  for(int j = 0 ; j < readsnps; j++) {
    if(verbose && j % 1000 == 0) Rcout << "\r" << j/1000 << "k SNPs read...";

    // read raw genotypes
    in.read((char*)tmp, sizeof(char) * bytespersnp);

    for(int l = 0 ; l < bytespersnp; ++l) {
      if(!blockused[l]) continue;

      tmpi = tmp[l];
      k = PACK_DENSITY * l;

      /* geno is interpreted as a char, however a1 and a2 are bits for allele 1 and
       * allele 2. The final genotype is the sum of the alleles, except for 11
       * which denotes missing.
       */
      tmp2[k] = (tmpi & MASK3) >> 6;
      tmp2[k+1] = (tmpi & MASK2) >> 4;
      tmp2[k+2] = (tmpi & MASK1) >> 2;
      tmp2[k+3] = (tmpi & MASK0);
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


// [[Rcpp::export]]
NumericMatrix cpp_read_eigenstrat(String genofile, int nsnp, int nind, IntegerVector indvec,
                                  int first, int last, bool transpose = false, bool verbose = true) {

  int val;
  long long len, bytespersnp;
  int readsnps = last - first;

  std::ifstream in(genofile.get_cstring(), std::ios::in | std::ios::binary);

  if(!in) {
    Rcout << "Error reading file " << genofile.get_cstring() << std::endl;
    throw std::runtime_error("io error");
  }
  in.seekg(0, std::ifstream::end);
  // file size in bytes
  len = (long long)in.tellg();
  bytespersnp = len/nsnp;

  int nindused = 0;
  for(int i = 0 ; i < nind; i++) {
    if(indvec[i] == 1) {
      nindused++;
    }
  }

  NumericMatrix geno(transpose?nindused:readsnps, transpose?readsnps:nindused);
  std::fill(geno.begin(), geno.end(), NA_REAL);
  in.seekg(first*bytespersnp, std::ifstream::beg);
  std::string tmp;
  for(int j = 0 ; j < readsnps; j++) {
    int c = 0;
    if(verbose && j % 1000 == 0) Rcout << "\r" << j/1000 << "k SNPs read...";
    std::getline(in, tmp);
    if(!transpose) {
      for(int i = 0; i < nind; i++) {
        if(!indvec[i]) continue;
        val = tmp[i]-'0';
        if(val != 9) geno(j, c) = val;
        c++;
      }
    } else {
      for(int i = 0; i < nind; i++) {
        if(!indvec[i]) continue;
        val = tmp[i]-'0';
        if(val != 9) geno(c, j) = val;
        c++;
      }
    }
  }
  if(verbose) Rcout << std::endl;
  in.close();
  return geno;
}

