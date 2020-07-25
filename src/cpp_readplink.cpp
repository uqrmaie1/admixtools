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
List cpp_read_plink_afs(String bedfile, const NumericVector indvec, const NumericVector indvec2,
                        bool adjust_pseudohaploid, bool verbose) {
  // indvec: assignes each individual to population; indvec2: which individuals to keep

  NumericMatrix X, afmat, countmat, scanmat;
  unsigned int N, p;
  unsigned long long len, filesize;
  unsigned int np, nsnps;

  N = indvec.length();
  p = 0;
  nsnps = 0;

  int val, k, pop, npop, nind;
  npop = max(indvec);
  nind = indvec2.length();
  std::ifstream in(bedfile.get_cstring(), std::ios::in | std::ios::binary);

  if(!in) {
    std::cerr << "[read_afs] Error reading file " << bedfile.get_cstring() << std::endl;
    throw std::runtime_error("io error");
  }
  in.seekg(0, std::ifstream::end);
  // file size in bytes, ignoring first 3 bytes (2byte magic number + 1byte mode)
  len = (long)in.tellg() - 3;
  // size of packed data, in bytes, per SNP
  np = (long)ceil((double)N / PACK_DENSITY);
  nsnps = len / np;
  in.seekg(3, std::ifstream::beg);
  unsigned char* tmp = new unsigned char[np];

  // Allocate more than the sample size since data must take up whole bytes
  unsigned char* tmp2 = new unsigned char[np * PACK_DENSITY];
  afmat = NumericMatrix(nsnps, npop);
  countmat = NumericMatrix(nsnps, npop);

  // if(verbose) std::cout << ">>> Detected BED file: " << this->bedfile <<
  //     " with " << len << " bytes, " << N << " samples, " << nsnps
  //              << " SNPs." << std::endl;
  NumericVector popafs(npop), popsum(npop), poptot(npop);

  // determine ploidy
  NumericVector ploidy(nind);
  ploidy.fill(1);
  int nscan = 1000;
  for(unsigned int j = 0; j < nscan; j++) {
    if(j == nsnps) break;
    in.read((char*)tmp, sizeof(char) * np);
    decode_plink(tmp2, tmp, np);
    for(unsigned int i = 0; i < nind; i++) {
      k = indvec2(i)-1;
      val = (double)tmp2[k];
      if(val == 1 || !adjust_pseudohaploid) ploidy(i) = 2;
      //ploidy(i) = 2;
    }
  }

  in.seekg(3, std::ifstream::beg);
  for(unsigned int j = 0 ; j < nsnps ; j++) {
    if(verbose && j % 1000 == 0) Rcout << "\r" << j/1000 << "k SNPs read...";
    popsum = NumericVector(npop);
    poptot = NumericVector(npop);

    // read raw genotypes
    in.read((char*)tmp, sizeof(char) * np);

    // decode the genotypes
    decode_plink(tmp2, tmp, np);
    for(unsigned int i = 0; i < nind; i++) {
      k = indvec2(i)-1;
      pop = indvec(k)-1;
      val = (double)tmp2[k];
      if(val == PLINK_NA) val = 0;
      else poptot(pop) += ploidy(i);
      popsum(pop) += val/(3-ploidy(i));
    }

    for(unsigned int i = 0; i < npop; i++) {
      if(poptot(i) > 0) popafs(i) = popsum(i)/poptot(i);
      else popafs(i) = NA_REAL;
    }
    afmat.row(j) = popafs;
    countmat.row(j) = poptot;
  }
  if(verbose) Rcout << std::endl;

  delete[] tmp;
  delete[] tmp2;

  in.close();

  return Rcpp::List::create(_["afmat"] = afmat, _["countmat"] = countmat);
}


// [[Rcpp::export]]
NumericMatrix cpp_read_plink(String bedfile, int nsnp, int nind, IntegerVector indvec,
                             int first, int last, bool transpose = false, bool verbose = true) {

  int val;
  long len, bytespersnp;
  int readsnps = last - first;

  std::ifstream in(bedfile.get_cstring(), std::ios::in | std::ios::binary);

  if(!in) {
    Rcout << "Error reading file " << bedfile.get_cstring() << std::endl;
    throw std::runtime_error("io error");
  }
  in.seekg(0, std::ifstream::end);
  // file size in bytes
  len = (long)in.tellg();
  bytespersnp = (long)ceil((double)nind / PACK_DENSITY);

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
  char* tmp = new char[bytespersnp + 1];
  tmp[bytespersnp] = '\0';
  char tmpi;

  // Allocate more than the sample size since data must take up whole bytes
  char* tmp2 = new char[bytespersnp * PACK_DENSITY + 1];
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



