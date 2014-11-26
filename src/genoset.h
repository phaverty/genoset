#ifndef GENOSET_H
#define GENOSET_H

// utils.c
void isNA(SEXP vec, char* na);
int numNA(SEXP vec, char* na); // Would rather overload (C++ only) with int isNA(SEXP vec, char* na);
void widthToStartEnd(int* width, size_t* start, size_t* end, size_t n);
void widthToStart(int* width, size_t* start, size_t n);
size_t leftBound(size_t* values, size_t low, size_t high, size_t query);
int* leftBoundPointer(int* low, int* high, int query);

// bounds.c
SEXP binary_bound(SEXP starts, SEXP stops, SEXP positions);
SEXP binary_bound_by_chr(SEXP nquery, SEXP query_chr_indices, SEXP query_starts, SEXP query_ends, SEXP query_names, SEXP subject_chr_indices, SEXP subject_starts, SEXP subject_ends);

// rangeSummaries.c
SEXP rangeMeans_rle(SEXP Start, SEXP End, SEXP RunValues, SEXP RunLengths, SEXP Na_rm);
SEXP rangeMeans_numeric(SEXP bounds, SEXP x, SEXP Na_rm);
SEXP numCallable_rle(SEXP Start, SEXP End, SEXP RunValues, SEXP RunLengths, SEXP Min);

#endif
