#ifndef GENOSET_H
#define GENOSET_H

// utils.c
void isNA(SEXP vec, char* na);
int numNA(SEXP vec, char* na); // Would rather overload (C++ only) with int isNA(SEXP vec, char* na);
void widthToStartEnd(int* width, size_t* start, size_t* end, size_t n);
void widthToStart(int* width, size_t* start, size_t n);
size_t leftBound(size_t* values, size_t low, size_t high, size_t query);
int* leftBoundPointer(int* low, int* high, int query);

#endif

