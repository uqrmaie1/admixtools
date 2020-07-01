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

  int val, bytespersnp;
  long len;
  int readsnps = std::min(nsnp, last?last:nsnp) - first;

  std::ifstream in(genofile.get_cstring(), std::ios::in | std::ios::binary);

  if(!in) {
    Rcout << "Error reading file " << genofile.get_cstring() << std::endl;
    throw std::runtime_error("io error");
  }
  in.seekg(0, std::ifstream::end);
  // file size in bytes
  len = (long)in.tellg();
  bytespersnp = len/(nsnp+1);

  // size of packed data, in bytes, per SNP
  //np = (long)ceil((double)nind / PACK_DENSITY);

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
  unsigned char* tmp = new unsigned char[bytespersnp];
  unsigned char tmpi;

  // Allocate more than the sample size since data must take up whole bytes
  unsigned char* tmp2 = new unsigned char[bytespersnp * PACK_DENSITY];

  int k;
  for(int j = 0 ; j < readsnps; j++) {
    //for(unsigned int j = 0 ; j < 3; j++) {
    if(verbose && j % 1000 == 0) Rcout << "\r" << j/1000 << "k SNPs read...";

    // read raw genotypes
    in.read((char*)tmp, bytespersnp);

    for(int l = 0 ; l < bytespersnp ; ++l) {
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





// // [[Rcpp::export]]
// List cpp_packedancestrymap_to_aftable(String genofile, int nind, int nsnp,
//                                       IntegerVector popind, bool verbose) {
//
//   int val, bytespersnp, np;
//   long len;
//
//   int npop = max(popind);
//   NumericMatrix afs(nsnp, npop);
//   NumericMatrix counts(nsnp, npop);
//   std::fill(afs.begin(), afs.end(), 0.0);
//   std::fill(counts.begin(), counts.end(), 0);
//
//   std::ifstream in(genofile.get_cstring(), std::ios::in | std::ios::binary);
//
//   if(!in) {
//     Rcout << "Error reading file " << genofile.get_cstring() << std::endl;
//     throw std::runtime_error("io error");
//   }
//   in.seekg(0, std::ifstream::end);
//   // file size in bytes
//   len = (long)in.tellg();
//   bytespersnp = len/(nsnp+1);
//   // size of packed data, in bytes, per SNP
//   np = (long)ceil((double)nind / PACK_DENSITY);
//
//   in.seekg(bytespersnp, std::ifstream::beg);
//   unsigned char* tmp = new unsigned char[np];
//   unsigned char tmpi;
//   int* seekind = new int[np];
//   double indperpop[npop];
//   for(int i = 0; i < npop; i++) {
//     indperpop[i] = 0.0;
//   }
//   for(int i = 0; i < nind; i++) {
//     indperpop[popind(i)-1]++;
//     //Rcout << "indperpop " << indperpop[popind(i)-1] << " popind(i) " <<  popind(i) << std::endl;
//   }
//
//   int seeknum = 0;
//   int induse = 0;
//   for(int i = 0; i < nind; i++) {
//     if(popind(i) != 0) {
//       induse++;
//       if(seeknum == 0 || i/4 != seekind[seeknum-1]) {
//         seekind[seeknum] = i/4;
//         seeknum++;
//       }
//     }
//   }
//   int* popinduse = new int[induse];
//   int c = 0;
//   for(int i = 0; i < nind; i++) {
//     if(popind(i) != 0) {
//       popinduse[c] = popind(i);
//       c++;
//     }
//   }
//   for(int i = 0; i < seeknum; i++) {
//     Rcout << seekind[i] << " ";
//   }
//   Rcout << std::endl;
//
//   unsigned char* tmp2 = new unsigned char[seeknum * PACK_DENSITY];
//
//   //in.seekg(bytespersnp, std::ifstream::beg);
//   for(int j = 0; j < nsnp; j++) {
//     if(verbose && j % 1000 == 0) Rcout << "\r" << j/1000 << "k SNPs read...";
//     int offset = (j+1) * bytespersnp;
//
//     int k, p, pop, nonmiss;
//     for(int l = 0 ; l < seeknum; l++) {
//       p = seekind[l];
//       in.seekg(offset + p, std::ifstream::beg);
//       in.read((char*)tmp, 1);
//       tmpi = tmp[0];
//       k = PACK_DENSITY * l;
//
//       tmp2[k] = (tmpi & MASK3) >> 6;
//       tmp2[k+1] = (tmpi & MASK2) >> 4;
//       tmp2[k+2] = (tmpi & MASK1) >> 2;
//       tmp2[k+3] = (tmpi & MASK0);
//     }
//     for(int i = 0 ; i < induse; i++) {
//       pop = popinduse[i];
//       //Rcout << "k+q " << k+q << std::endl;
//       //Rcout << "pop " << pop << std::endl;
//       if(pop != 0) {
//         val = (double)tmp2[i];
//         nonmiss = val != 3;
//         //if(j == 0) Rcout << "indperpop[pop-1] " << indperpop[pop-1] << " pop " << pop << std::endl;
//         counts(j, pop-1) += nonmiss;
//         if(nonmiss) afs(j, pop-1) += val/2.0/indperpop[pop-1];
//       }
//     }
//
//   }
//   if(verbose) Rcout << std::endl;
//
//   delete[] tmp;
//   delete[] tmp2;
//
//   in.close();
//
//   return Rcpp::List::create(_["afs"] = afs, _["counts"] = counts);
// }


