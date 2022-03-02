#include <stdio.h>
#include <stdlib.h>
#include "silo.h"

int main(int argc, char* argv[]){
  if (argc != 4) {
    printf("Usage: %s silo_file variable output_file\n", argv[0]);
    exit(1);
  }

  /* Command line arguments */
  const char *silo_file = argv[1];
  const char *variable_name = argv[2];
  const char *output_file = argv[3];

  DBfile* dbfile = DBOpen(silo_file, DB_UNKNOWN, DB_READ);
  DBmultivar *mvar = DBGetMultivar(dbfile, variable_name);

  FILE *output = fopen(output_file, "wb");
  fwrite(&mvar->nvars, sizeof(int), 1, output);

  for (int i = 0; i < mvar->nvars; i++) {
    DBquadvar* dbqv = DBGetQuadvar(dbfile, mvar->varnames[i]);
    DBquadmesh* dbqm = DBGetQuadmesh(dbfile, dbqv->meshname);

    if (dbqm->datatype != DB_DOUBLE || dbqv->datatype != DB_DOUBLE) {
      fprintf(stderr, "Expected datatype DB_DOUBLE");
      exit(1);
    }

    fwrite(&dbqm->ndims, sizeof(int), 1, output);
    fwrite(dbqm->dims, sizeof(int), dbqm->ndims, output);
    fwrite(dbqm->min_index, sizeof(int), dbqm->ndims, output);
    fwrite(dbqm->max_index, sizeof(int), dbqm->ndims, output);

    for (int d=0; d < dbqm->ndims; d++) {
      fwrite(dbqm->coords[d], sizeof(double), dbqm->dims[d], output);
    }

    /* Write scalar data */
    fwrite(dbqv->vals[0], sizeof(double), dbqv->nels, output);

    DBFreeQuadvar(dbqv);
    DBFreeQuadmesh(dbqm);
  }

  DBFreeMultivar(mvar);
  DBClose(dbfile);
  fclose(output);

  return 0;
}
