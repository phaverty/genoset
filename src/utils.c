#include <R.h>
#include <Rinternals.h>
#include "genoset.h"

// Take an SEXP and a pointer to a char array, fill out the char array 
//   with 0 (not NA) or 1 (NA) by looking at the values in vec
// Then we can compute on it without branching or polymorphism to account 
//   for type, Julia DataArrays style.
// This burns a bit of memory, but hopefully avoiding branching code for
//   NA checking will save a lot of time.
// Use the resulting NAs in a mathematical context rather than branching:
//   y += na[i] * x;
//   y += (! na[i]) * x;
// The compiler will pipeline the math for you, branches would prevent that.
// Argh, just noticed that Hadley thought of doing an isna C function, but his 
//   is inside out and returns INTSXP.
//  This could alternatively be implemented as vec<bool>. Probably very similar speed, but smaller.
inline void isNA(SEXP vec, char* na) {
  if (TYPEOF(vec) == REALSXP) {
    double* vec_p = REAL(vec);    
    for (int i=0; i < LENGTH(vec); i++) {
      na[i] = isnan(vec_p[i]);
    }
  } else if (TYPEOF(vec) == INTSXP) {
    int* vec_p = INTEGER(vec);
    for (int i=0; i < LENGTH(vec); i++) {
      na[i] = vec_p[i] == NA_INTEGER;
    }
  } else if (TYPEOF(vec) == LGLSXP) {
    int* vec_p = INTEGER(vec);
    for (int i=0; i < LENGTH(vec); i++) {
      na[i] = vec_p[i] == NA_LOGICAL;
    }
  } else if (TYPEOF(vec) == STRSXP) {
    for (int i=0; i < LENGTH(vec); i++) {
      na[i] = STRING_ELT(vec, i) == NA_STRING; 
    }
  } else {
    error("vec must contain either 'integer', 'logical' or 'character' or 'numeric' values");
  }
}

// Like the above isNA, but returns the number of NAs.
//   This may be used to decide proceed with a simpler algorithm 
//   if == 0 or the sum itself may be used to simplify subsequent
//   code, like a mean.
inline int numNA(SEXP vec, char* na) {
  int num_na = 0;
  if (TYPEOF(vec) == REALSXP) {
    double* vec_p = REAL(vec);    
    for (int i=0; i < LENGTH(vec); i++) {
      na[i] = ISNA(vec_p[i]);
      num_na += na[i];
    }
  } else if (TYPEOF(vec) == INTSXP) {
    int* vec_p = INTEGER(vec);
    for (int i=0; i < LENGTH(vec); i++) {
      na[i] = vec_p[i] == NA_INTEGER;
      num_na += na[i];
    }
  } else if (TYPEOF(vec) == LGLSXP) {
    int* vec_p = INTEGER(vec);
    for (int i=0; i < LENGTH(vec); i++) {
      na[i] = vec_p[i] == NA_LOGICAL;
      num_na += na[i];
    }
  } else if (TYPEOF(vec) == STRSXP) {
    for (int i=0; i < LENGTH(vec); i++) {
      na[i] = STRING_ELT(vec, i) == NA_STRING; 
      num_na += na[i];
    }
  } else {
    error("vec must contain either 'integer', 'logical' or 'character' or 'numeric' values");
  }
  return(num_na);
}

// Take widths, like from an Rle, compute start and end index for each run (1-based).
inline void widthToStartEnd(int* width, size_t* start, size_t* end, size_t n) {
  start[0] = 1;
  end[0] = width[0];
  for (int i=1; i < n; i++) {
    start[i] = end[i-1] + 1;
    end[i] = end[i-1] + width[i];
  }
}

// Take widths, like from an Rle, compute start index for each run (1-based).
inline void widthToStart(int* width, size_t* start, size_t n) {
  start[0] = 1;
  for (int i=1; i < n; i++) {
    start[i] = start[i-1] + width[i-1];
  }
}

/***********************
Searching
***********************/
// Fun fact: if you do --values before you pass it in, you get 1-based indices back, a la findInterval.
// To take advantage of a searching for a previous value in a sorted array, start 'low' at that location

inline size_t leftBound(size_t* values, size_t low, size_t high, size_t query) {
  // Right bound likely close to previous (low)
  size_t probe = low + 1;
  size_t jump = 2;
  //printf("cast pre low: %u, probe%u, high:%u, values[probe]: %u, query: %u\n", low, probe, high, values[probe], query);
  while (probe <= high && values[probe] <= query) {
    low = probe;
    probe += jump;
    jump = jump << 1;
  }
  high = probe > high ? high + 1 : probe;
  probe = ((high-low) >> 1) + low; // Hack to avoid overflow 
  while (low < probe) {
    if (values[probe] > query) {
      high = probe;
    } else {
      low = probe;
    }
    probe = ((high-low) >> 1) + low; // Hack to avoid overflow 
  }
  return(low);
}

// Like leftBound above, but works on just pointers. Not any faster. Convenient if you don't need the indices.
inline int* leftBoundPointer(int* low, int* high, int query) {
  // Right bound likely close to previous
  int* probe = low + 1;
  int jump = 2;
  while (probe <= high && *probe <= query) {
    low = probe;
    probe += jump;
    jump = jump << 1;
  }
  high = probe > high ? high + 1 : probe;
  // Now binary search for closest left bound
  probe = ((high-low) >> 1) + low; // Hack to avoid overflow 
  while (low < probe) {
    if (*probe > query) {
      high = probe;
    } else {
      low = probe;
    }
    probe = ((high-low) >> 1) + low; // Hack to avoid overflow 
  }
  return(low);
}
