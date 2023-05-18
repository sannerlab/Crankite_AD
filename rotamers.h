//
// May 2023
// Author Michel Sanner
//
// code for dynamically loading side chain rotamers into ADCP
// avoiding having to recompile to C code to add new types
// the rotamer description is read from rotamers.lib
//
// the files rotamers.[ch] replace canonicalAA.[ch] where the
// _AASCRotTable rable was hardwired with a capacity of 100 entries
//
// structure for full atoms rotameric side chains
#ifndef ROTAMERDEF
#define ROTAMERDEF

#include "params.h"

struct _AASCRot {
  int nbRot;
  int nbAtoms;
  int *atypes;        /* list of nbAtoms integers describing atom types */
  char *name;         /* residue name */
  char *coarse_type;  /* name of coarse potential entry */
  double *charges;    /* list of nbAtoms doubles holding atomic charges */
  double *rotProbas;  /* list of probabilities for rotamers, adds up to 1.0 */
  char *atnames;
  double **coords;   /* [nbRot][natoms x 3] double holding atomic coordinates for rotamers */
};

extern int getSideChainTemplateIndexFromName(char *str);
extern int getSideChainTemplateIndexFromIDchar(char id);
extern int countRotamers(char *filename);
extern int initialize_AASCRotTable_from_file(char *filename, int lastIndex);
extern int nbCanAA;

extern struct _AASCRot *_AASCRotTable;
#endif
