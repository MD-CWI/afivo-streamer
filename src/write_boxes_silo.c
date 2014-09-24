/*
Compile:
gcc write_boxes_silo.c -o write_boxes_silo -I../silo/include -L../silo/lib -lsiloh5 -lm -Wall -std=c99
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "silo.h"

const int NN      = 9;
const int ndims   = 2;
const int nmeshes = 125;

void get_coords(double *coords[], const int *ix, const double dr) {
  double r_min[2];
  r_min[0] = ix[0] * (NN-1) * dr;
  r_min[1] = ix[1] * (NN-1) * dr;

  for(int d=0; d<ndims; d++) {
    for(int i=0; i<NN; i++) {
      coords[d][i] = r_min[d] + i * dr;
    }
  }
}

int main(void) {
  DBfile    *file = NULL;
  double    *coords[ndims];
  int        dims[ndims];
  int        ierr;
  int        ix[2];
  int        meshtype;
  int        dummy;
  char       gridname[10];
  DBoptlist *optlist;
  double    *extents;
  int       *zonecounts;
  int       *extzones;

  extents = (double*) malloc(2*nmeshes*nmeshes*ndims*sizeof(double));
  zonecounts = (int*) malloc(nmeshes*nmeshes*sizeof(int));
  extzones = (int*) malloc(nmeshes*nmeshes*sizeof(int));

  for(int d=0; d<ndims; d++) {
    dims[d] = NN;
    coords[d] = (double*) malloc(NN * sizeof(double));
  }

  file = DBCreate("test.silo", DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5);
  DBMkDir(file, "foo");
  DBSetDir(file,"foo");

  optlist = DBMakeOptlist(10);
  dummy = 1;
  DBAddOption(optlist, DBOPT_HIDE_FROM_GUI, &dummy);

  for(int i=0; i<nmeshes; i++) {
    for(int j=0; j<nmeshes; j++) {
      ix[0] = i;
      ix[1] = j;
      get_coords(coords, ix, 0.1);
      zonecounts[i*nmeshes+j] = (NN-1)*(NN-1);
      extzones[i*nmeshes+j] = 0;
      dummy = 4*(i*nmeshes+j);
      extents[dummy] = coords[0][0];
      extents[dummy+1] = coords[1][0];
      extents[dummy+2] = coords[0][NN-1];
      extents[dummy+3] = coords[1][NN-1];
      sprintf(gridname, "grid%d", i*nmeshes+j);
      ierr = DBPutQuadmesh(file, gridname, NULL, coords, dims, ndims,
                           DB_DOUBLE, DB_COLLINEAR, optlist);
    }
  }
  DBFreeOptlist(optlist);

  DBSetDir(file, "..");

  optlist = DBMakeOptlist(10);
  meshtype = DB_QUADMESH;
  DBAddOption(optlist, DBOPT_MB_BLOCK_NS, "|foo/grid%d|n");
  DBAddOption(optlist, DBOPT_MB_BLOCK_TYPE, &meshtype);
  dummy = 4;
  DBAddOption(optlist, DBOPT_EXTENTS_SIZE, &dummy);
  DBAddOption(optlist, DBOPT_EXTENTS, extents);
  DBAddOption(optlist, DBOPT_ZONECOUNTS, zonecounts);
  DBAddOption(optlist, DBOPT_HAS_EXTERNAL_ZONES, extzones);
  ierr = DBPutMultimesh(file, "bigmesh", nmeshes*nmeshes, NULL, NULL, optlist);
  DBFreeOptlist(optlist);

  DBClose(file);
  printf("written test.silo, ierr = %d\n", ierr);
  printf("number of cells %d\n", nmeshes * nmeshes * (NN-1) * (NN-1));
  return 0;
}
