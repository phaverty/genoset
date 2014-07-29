#include <R.h>
#include <Rinternals.h>

/* TODO: break out redundant code in binary_bound and binary_bound_by_chr for ease of reading/maintenance */

SEXP binary_bound(SEXP starts, SEXP stops, SEXP positions) {
  int query_index, probe, left, right, low, high, jump;

  // Get what we need from input
  int* starts_p = INTEGER(starts);
  int* stops_p = INTEGER(stops);
  int* positions_p = INTEGER(positions);
  int num_queries = LENGTH(starts);
  int num_positions = LENGTH(positions);
  --positions_p;  // Hack for 1-based indices

  // Make results matrix
  int num_protected = 0;
  SEXP bounds, dimnames, colnames;
  PROTECT(bounds = allocMatrix(INTSXP, num_queries, 2)); num_protected++;
  PROTECT(dimnames = allocVector(VECSXP, 2)); num_protected++;
  PROTECT(colnames = allocVector(STRSXP, 2)); num_protected++;
  SET_VECTOR_ELT(dimnames, 0, R_NilValue);
  SET_STRING_ELT(colnames, 0, mkChar("left"));
  SET_STRING_ELT(colnames, 1, mkChar("right"));
  SET_VECTOR_ELT(dimnames, 1, colnames);
  setAttrib(bounds, R_DimNamesSymbol, dimnames);
  int *bounds_p = INTEGER(bounds);

  // Initialize
  low = 1;
  high = num_positions; /* Set high to off right end */
  
  for (query_index=0; query_index < num_queries; query_index++) {
   
    /* Left bound */
    left = starts_p[query_index];

    /* If data unsorted, current target may be anywhere left of low for previous target, just start at left */
    if (low > 0 && left < positions_p[low]) {
      low = 1;
    }
    
    /* Right bound likely close to high from previous gene */
    for (jump=1; ; jump*=2) {
      high += jump;
      if (high >= num_positions) {
	high = num_positions;
	break;
      }
      if (left < positions_p[high]) {  /* Note difference to similar code resetting high after first binary search */
	break;
      }
      low = high;
    }

    /* Now binary search for right bound */
    while (high - low > 1) {
      probe = (high + low) >> 1;
      if (positions_p[probe] > left) {
        high = probe;
      } else {
	low = probe;
      }
    }
    bounds_p[query_index] = low;
    
    /* Right bound */
    right = stops_p[query_index];

    /* Right bound likely close to left bound, relative to length of positions, so expand exponentially, a la findInterval */
    for (jump=1; ; jump*=2) {
      high += jump;
      if (high >= num_positions) {
	high = num_positions;
	break;
      }
      if (right <= positions_p[high]) {
	break;
      }
      low = high;
    }

    /* Now binary search for right bound */
    while (high - low > 1) {
      probe = (high + low) >> 1;
      if (positions_p[probe] < right) {
        low = probe;
      } else {
        high = probe;
      }
    }
    bounds_p[query_index + num_queries] = high; /* Right bound goes in second column of bounds */
  
    low = bounds_p[query_index]; /* Reset low to left end of this query to start next query */
    
  } /* End foreach query */
  UNPROTECT(num_protected);
  return(bounds);
}

SEXP binary_bound_by_chr(SEXP nquery, SEXP query_chr_indices, SEXP query_starts, SEXP query_ends, SEXP query_names, SEXP subject_chr_indices, SEXP subject_starts, SEXP subject_ends) {

  /* Pointers into data in arg objects */
  double *p_query_chr_indices = REAL(query_chr_indices);
  double *p_subject_chr_indices = REAL(subject_chr_indices);
  int *p_query_starts = INTEGER(query_starts);
  int *p_query_ends = INTEGER(query_ends);
  int *p_subject_starts = INTEGER(subject_starts);
  int *p_subject_ends = INTEGER(subject_ends);
  p_query_starts--; p_query_ends--;  // -- gives 1-based indexing
  p_subject_starts--; p_subject_ends--;  // -- gives 1-based indexing

  /* Make results matrix, "bounds" */
  int num_protected = 0;
  SEXP bounds, dimnames, rownames, colnames;
  PROTECT(bounds = allocMatrix(INTSXP, INTEGER(nquery)[0], 2)); num_protected++;
  PROTECT(dimnames = allocVector(VECSXP, 2)); num_protected++;
  PROTECT(rownames = allocVector(STRSXP, INTEGER(nquery)[0])); num_protected++;
  PROTECT(colnames = allocVector(STRSXP, 2)); num_protected++;
  SET_VECTOR_ELT(dimnames, 0, rownames);
  SET_STRING_ELT(colnames, 0, mkChar("left"));
  SET_STRING_ELT(colnames, 1, mkChar("right"));
  SET_VECTOR_ELT(dimnames, 1, colnames);
  setAttrib(bounds, R_DimNamesSymbol, dimnames);
  int *p_bounds = INTEGER(bounds);

  int low, high, probe, jump;  /* markers for binary search */
  int init_low, init_high;  /* index bounds for search in chr */
  int left, right;  /* values to search for */

  /* markers for row in two cols of result matrix, bounds */
  int left_bound_index = 0;
  int right_bound_index = INTEGER(nquery)[0];

  /* foreach chr */
  for (int chr_start_index=0, chr_end_index=nrows(query_chr_indices); chr_start_index < nrows(query_chr_indices); chr_start_index++, chr_end_index++) {
    /* binary_bound all queries in chr*/
    /* How best to handle being off each end of chr?  Test upfront?  Test at end? */

    init_low = low = p_subject_chr_indices[chr_start_index];
    init_high = high = p_subject_chr_indices[chr_end_index];

    /* This loop should work on query_index in range specified by query_chr_indices for this chr, left and right bound indices must increment every loop but be set at the top */
    for (int query_index=p_query_chr_indices[chr_start_index]; query_index <= p_query_chr_indices[chr_end_index]; query_index++, left_bound_index++, right_bound_index++) {
      SET_STRING_ELT(rownames,left_bound_index, STRING_ELT(query_names,query_index-1));
      left = p_query_starts[query_index];
      right = p_query_ends[query_index];

      // Check for being completely off either end
      if (right <= p_subject_starts[init_low]) {
	p_bounds[right_bound_index] = init_low;
	p_bounds[left_bound_index] = init_low;
	continue;
      }
      if (left >= p_subject_ends[init_high]) {
	p_bounds[right_bound_index] = init_high;
	p_bounds[left_bound_index] = init_high;
	continue;
      }

      /* Left bound */
      /* If data unsorted, current target may be anywhere left of low for previous target, just start at left */
      if (left < p_subject_starts[low]) {
	low = init_low;
      }

      /* Right bound likely close to high from previous gene */
      for (jump=1; ; jump*=2) {
	high += jump;
	if (high >= init_high) {
	  high = init_high;
	  break;
	}
	if (left < p_subject_starts[high]) {  /* Note difference to similar code resetting high after first binary search */
	  break;
	}
	low = high;
      }

      /* Now binary search for left bound */
      while (high - low > 1) {
	probe = (high + low) >> 1;

	if (p_subject_starts[probe] > left) {
	  high = probe;
	} else {
	  low = probe;
	}
      }
      p_bounds[left_bound_index] = low;

      if (right == p_subject_ends[low]) { p_bounds[right_bound_index] = low;  continue; } /* Hack for query range length 1 */

      /* Right bound */
      /* Right bound likely close to left bound, relative to length of positions, so expand exponentially, a la findInterval */
      for (jump=1; ; jump*=2) {
	high += jump;
	if (high >= init_high) {
	  high = init_high;
	  break;
	}
	if (right <= p_subject_ends[high]) {
	  break;
	}
	low = high;
      }

      /* Now binary search for right bound */
      while (high - low > 1) {
	probe = (high + low) >> 1;
	if (p_subject_ends[probe] < right) {
	  low = probe;
	} else {
	  high = probe;
	}
      }
      p_bounds[right_bound_index] = high; /* Right bound goes in second column of bounds */
  
      low = p_bounds[left_bound_index]; /* Reset low to left end of this query to start next query */
    
    } /* End foreach query in chr*/
  } /* End foreach chr */
  UNPROTECT(num_protected);
  return(bounds);
}
