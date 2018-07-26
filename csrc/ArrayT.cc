#include <cstdio>
#include <cerrno>
#include "ArrayT.h"

ArrayT::ArrayT (int sz, int npart){
   totSize = sz;
   theSize = 0;
   lastPos = 0;
   theFlg = 0;
   theArray = nullptr;
   offset = new long[npart];
   for (int i=0; i < npart; i++) offset[i]=0;
   if (sz > 0){
      theArray =  (int *) malloc (totSize*sizeof(int));
      if (theArray == nullptr){
         throw runtime_error("memory:: ArrayT");
      }
   }
}

ArrayT::~ArrayT(){
   if (theArray) {
      free(theArray);
   }
   delete [] offset;
   theArray = nullptr;
}






