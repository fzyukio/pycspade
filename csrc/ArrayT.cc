#include <cstdio>
#include <cerrno>
#include "ArrayT.h"

ArrayT::ArrayT (int sz, int npart){
   totSize = sz;
   theSize = 0;
   lastPos = 0;
   theFlg = 0;
   //theIncr = incr;
   theArrayT = NULL;
   offset = new long[npart];
   for (int i=0; i < npart; i++) offset[i]=0;
   if (sz > 0){
      theArrayT =  (int *) malloc (totSize*sizeof(int));
      //theArrayT = new int [totSize];
      if (theArrayT == NULL){
         throw std::runtime_error("memory:: ArrayT");
      }
   }
}

ArrayT::~ArrayT(){
   if (theArrayT) {
      free(theArrayT);
      //delete [] theArrayT;
   }
   delete [] offset;
   theArrayT = NULL;
}






