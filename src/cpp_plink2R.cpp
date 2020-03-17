/*
 *
 * This is based on Gad Abraham's "plink2R" package:
 * https://github.com/gabraham/plink2R/blob/master/plink2R/
 * Modified to return per group allele frequencies rather than raw genotypes.
 *
 */

#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen]]

using namespace Rcpp;
using namespace Eigen;

//#include <fcntl.h>
//#include <stdexcept>
//#pragma once

//#include <iostream>
//#include <fstream>
//#include <vector>
//#include <algorithm>

//#include <Eigen/Core>
//#include <Eigen/Dense>
//#include <Eigen/Eigen>
//#include <Eigen/SVD>
//#include <Eigen/Sparse>

#define PACK_DENSITY 4
#define PLINK_NA 3

/* 3 is 11 in binary, we need a 2 bit mask for each of the 4 positions */
#define MASK0 3	  /* 3 << 2 * 0 */
#define MASK1 12  /* 3 << 2 * 1 */
#define MASK2 48  /* 3 << 2 * 2 */
#define MASK3 192 /* 3 << 2 * 3 */

#define BUFSIZE 100

class Data {
public:

  MatrixXd X, afmat, countmat, scanmat;
  unsigned int N, p;
  unsigned long long len, filesize;
  unsigned int np, nsnps;
  unsigned char* data;
  const char *bedfile, *famfile, *bimfile;
  bool verbose;
  NumericVector indvec, indvec2;

  Data(const char* bedfile, const char* famfile, bool verbose);
  Data(const char* bedfile, const NumericVector indvec, const NumericVector indvec2, bool verbose);
  ~Data();
  void read_bed();
  void read_afs();
};



Data::Data(const char* bedfile, const NumericVector indvec, const NumericVector indvec2, bool verbose)
{
  srand48(time(NULL));
  N = indvec.length();
  p = 0;
  nsnps = 0;
  this->bedfile = bedfile;
  this->indvec = indvec;
  this->indvec2 = indvec2;
  this->verbose = verbose;
}

Data::~Data() {}

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
                  const unsigned char *in, const unsigned int n)
{
  unsigned int i, k;
  unsigned char tmp, geno;
  unsigned int a1, a2;

  for(i = 0 ; i < n ; ++i)
  {
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


// Expects PLINK BED in SNP-major format
void Data::read_afs()
{
  int val, k, pop, npop, nind;
  npop = max(this->indvec);
  //nindtot = this->indvec.length();
  nind = this->indvec2.length();
  std::ifstream in(this->bedfile, std::ios::in | std::ios::binary);

  if(!in)
  {
    std::cerr << "[read_afs] Error reading file "
              << this->bedfile << std::endl;
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
  afmat = MatrixXd(nsnps, npop);
  countmat = MatrixXd(nsnps, npop);

  // if(verbose) std::cout << ">>> Detected BED file: " << this->bedfile <<
  //     " with " << len << " bytes, " << N << " samples, " << nsnps
  //              << " SNPs." << std::endl;
  VectorXd popafs(npop), popsum(npop), poptot(npop);

  // determine ploidy
  VectorXd ploidy = VectorXd::Ones(nind);
  int nscan = 1000;
  for(unsigned int j = 0; j < nscan; j++) {
    if(j == nsnps) break;
    in.read((char*)tmp, sizeof(char) * np);
    decode_plink(tmp2, tmp, np);
    for(unsigned int i = 0; i < nind; i++) {
      k = this->indvec2(i)-1;
      val = (double)tmp2[k];
      if(val == 1) ploidy(i) = 2;
    }
  }

  in.seekg(3, std::ifstream::beg);
  for(unsigned int j = 0 ; j < nsnps ; j++)
  {
    if(verbose && j % 1000 == 0) Rcout << "\r" << j/1000 << "k SNPs read...";
    popsum = VectorXd::Zero(npop);
    poptot = VectorXd::Zero(npop);

    // read raw genotypes
    in.read((char*)tmp, sizeof(char) * np);

    // decode the genotypes
    decode_plink(tmp2, tmp, np);
    for(unsigned int i = 0; i < nind; i++) {
      k = this->indvec2(i)-1;
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
}

// [[Rcpp::export]]
List read_plink_afs_cpp(String bedfile, const NumericVector indvec, const NumericVector indvec2, bool verbose)
{
  // indvec: assignes each individual to population; indvec2: which individuals to keep
  Data data(bedfile.get_cstring(), indvec, indvec2, verbose);
  data.read_afs();
  return List::create(data.afmat, data.countmat);
}


// Data::Data(const char* bedfile, const char* famfile, bool verbose)
// {
//   srand48(time(NULL));
//   N = 0;
//   p = 0;
//   nsnps = 0;
//   this->bedfile = bedfile;
//   this->famfile = famfile;
//   this->verbose = verbose;
//   if(verbose)
//     std::cout << "bedfile: " << bedfile
//               << " famfile: " << famfile << std::endl;
//     this->mask = VectorXd::Ones(N).array() == 1.0;
// }

// Expects PLINK BED in SNP-major format
// void Data::read_bed()
// {
//   if(verbose)
//     std::cout << ">>> Reading BED file '" << this->bedfile << "'" << std::endl;
//   std::ifstream in(this->bedfile, std::ios::in | std::ios::binary);
//
//   if(!in)
//   {
//     std::cerr << "[read_bed] Error reading file "
//               << this->bedfile << std::endl;
//     throw std::runtime_error("io error");
//   }
//
//   in.seekg(0, std::ifstream::end);
//
//   // file size in bytes, ignoring first 3 bytes (2byte magic number + 1byte mode)
//   len = (unsigned int)in.tellg() - 3;
//
//   // size of packed data, in bytes, per SNP
//   np = (unsigned int)ceil((double)N / PACK_DENSITY);
//   nsnps = len / np;
//   in.seekg(3, std::ifstream::beg);
//
//   unsigned char* tmp = new unsigned char[np];
//
//   // Allocate more than the sample size since data must take up whole bytes
//   unsigned char* tmp2 = new unsigned char[np * PACK_DENSITY];
//   X = MatrixXd(N, nsnps);
//
//   if(verbose)
//     std::cout << ">>> Detected BED file: " << this->bedfile <<
//       " with " << len << " bytes, " << N << " samples, " << nsnps
//                << " SNPs." << std::endl;
//       VectorXd tmp3(N);
//
//       // The Rcpp code in read_plink will take care of converting PLINK_NA to
//       // NA_REAL
//       for(unsigned int j = 0 ; j < nsnps ; j++)
//       {
//         // read raw genotypes
//         in.read((char*)tmp, sizeof(char) * np);
//
//         // decode the genotypes
//         decode_plink(tmp2, tmp, np);
//
//         for(unsigned int i = 0 ; i < N ; i++)
//           tmp3(i) = (double)tmp2[i];
//         X.col(j) = tmp3;
//       }
//
//       p = X.cols();
//
//       delete[] tmp;
//       delete[] tmp2;
//
//       in.close();
// }

// NumericMatrix read_plink_cpp(String bedfile, String famfile, bool verbose)
// {
//   if(verbose) Rcout << "[read_plink] bedfile: " << bedfile.get_cstring() <<
//       " famfile:" << famfile.get_cstring() << std::endl;
//
//   Data data(bedfile.get_cstring(), famfile.get_cstring(), verbose);
//   data.read_bed();
//
//   NumericMatrix X(wrap(data.X));
//   NumericVector X2 = wrap(ifelse(X == PLINK_NA, NA_REAL, X));
//   NumericMatrix X3(X.nrow(), X.ncol(), X2.begin());
//   return X3;
// }

