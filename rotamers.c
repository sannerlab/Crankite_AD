//
// May 2023
// Author Michel Sanner
//
// code for dynamically loading side chain rotamers into ADCP
// avoiding having to recompile to ccode to add new types
// the rotamer description is read from rotamers.lib
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "params.h"
#include "rotamers.h"
#include "params.h"
#include "vector.h"
#include "rotation.h"
#include "peptide.h"
#include "energy.h"

// dynamically allocated table holding pointers to AASCRot
// describing full atom AA side chain rotamers
struct _AASCRot *_AASCRotTable;

// number of entries in the table _AASCRotTable
int nbCanAA = 0;

// get index of table entry for a given rotamer name
int getSideChainTemplateIndexFromName(char *str) {
  for (int i=0; i<nbCanAA; i++) {
    //printf("FUGU %s %d\n", str, i);
    if (strcmp(str, _AASCRotTable[i].name)==0) {
      //printf("FUGU %s %d\n", str, i);
      return i;
    }
  }
  return -1;
}

// 1 letter code for std AA
char IDchars[18] = "RNDCQEHILKMFPSTWYV";

char * aanames[] = {
  "ARG","ASN","ASP","CYS","GLN","GLU","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
};

  
int getSideChainTemplateIndexFromIDchar(char id) {
  char *c = IDchars;
  for (int i=0; i<18; i++) {
    if (c[i]==id) {
      //printf("FAGA %c %s %d\n", id, aanames[i]);
      int ind = getSideChainTemplateIndexFromName(aanames[i]);
      //printf("FAGA %c %s %d\n", id, aanames[i], ind);
      return ind;
    }
  }
  return -1;
}

// count rotamers
int countRotamers(char *filename) {
  FILE *in_file = fopen(filename, "r");
  char line[100];
  char *retval;
  int nbAA=0, firstNonSpace=0;
  
  // test for files not existing.
  if (in_file == NULL)
    {  
      printf("Error! countRotamers: Could not open file %s\n",filename);
      exit(-1);
    }
  while ( 1 ) 
    {
      retval = fgets(line, 100, in_file );
      if (retval == NULL ) break;
      for (firstNonSpace=0; line[firstNonSpace]==' '; firstNonSpace++) {};
      if (line[firstNonSpace]=='/') continue;
      if (strlen(&line[firstNonSpace])==0) continue;
      if (strncmp(&line[firstNonSpace], "rotamer", 7)==0) nbAA++;
    }
  fclose(in_file);
  return nbAA;
}

// read one rotamer description entry from the file and populate the
// structure
int getRotamerForAA(FILE *in_file, struct _AASCRot *aarot, char *aaname)
{
  char line[512], coarsePot[20];
  int firstNonSpace=0;
  int read;
  char *fgetret;
  
  fgetret = fgets(line, sizeof(line), in_file);
  while ((fgetret != NULL)) {
    /* printf("Retrieved line %s:\n", line); */
    /* printf("%s", line); */
    if (strlen(line)==0) {
      /* printf("failed read %s\n", line); */
      return 0; // failed to read
    }
    for (firstNonSpace=0; line[firstNonSpace]==' '; firstNonSpace++) {};
    if (strncmp(&line[firstNonSpace], "rotamer", 7)==0) break;
    fgetret = fgets(line, sizeof(line), in_file);
  }
  // read rotamer declaration line
  read = sscanf(&line[firstNonSpace], "rotamer %s %s %d %d", aaname, coarsePot, &aarot->nbRot, &aarot->nbAtoms);
  // allocate memory for rotamer name, atom names, types, charges and coordinates pointer for
  // rotamers
  aarot->name = (char *)malloc((strlen(aaname)+1)*sizeof(char));
  aarot->rotProbas = (double *)malloc(aarot->nbRot*sizeof(double));
  aarot->coarse_type = (char *)malloc((strlen(coarsePot)+1)*sizeof(char));
  aarot->atypes = (int *)malloc(aarot->nbAtoms*sizeof(int));
  aarot->charges = (double *)malloc(aarot->nbAtoms*sizeof(double));
  aarot->atnames = (char *)malloc(1+aarot->nbAtoms*5*sizeof(char));
  aarot->coords = (double **)malloc(aarot->nbRot*sizeof(double*));

  // fill out the structure
  strcpy(aarot->name, aaname); // rotamer name
  strcpy(aarot->coarse_type, coarsePot); // rotamer name

  // read the rotamer probabilities
  for (int i=0; i<aarot->nbRot; i++) {
    read = fscanf(in_file, "%lf", &aarot->rotProbas[i]);
    if (read!=1) return 0;
  }
  fscanf(in_file, "%c", line); while (line[0]!='\n') fscanf(in_file, "%c", line); // get rid of trailing spaces and NL

  // read the atom types
  fgetret = fgets(line, sizeof(line), in_file);
  char *p = NULL;
  p = strtok(line, " ");
  if (p != NULL) p[strcspn(p, "\r\n")] = 0; // replace LF, CR, CRLF, LFCR, ... with 0
  for (int i=0; i<aarot->nbAtoms && p != NULL; i++) {
    // printf("atomtype str=[%s]\n", p);
    for (int j=0; j<MAX_ATOM_TYPES; j++)
      {
	if (strcmp(p, atypes[j])==0) {
	  aarot->atypes[i] = j;
	  break;
	}
      }
    p = strtok(NULL," ");
    if (p != NULL) p[strcspn(p, "\r\n")] = 0; // replace LF, CR, CRLF, LFCR, ... with 0
  }

  // read the atom types
  // for (int i=0; i<aarot->nbAtoms; i++) {
  //  read = fscanf(in_file, "%d", &aarot->atypes[i]);
  //  if (read!=1) return 0;
  // }

  // read the partial charges
  for (int i=0; i<aarot->nbAtoms; i++) {
    read = fscanf(in_file, "%lf", &aarot->charges[i]);
    if (read!=1) return 0;
  }
  fscanf(in_file, "%c", line); while (line[0]!='\n') fscanf(in_file, "%c", line); // get rid of trailing spaces and NL

  // read the atom names
  for (int i=0; i<aarot->nbAtoms; i++) {
    read = fscanf(in_file, "%5c", &aarot->atnames[i*5]);
    //printf("read %ld [%s]\n", read, &aarot->atnames[i*5]);
  }
  fscanf(in_file, "%c", line); while (line[0]!='\n') fscanf(in_file, "%c", line); // get rid of trailing spaces and NL

  // read the atomic coordinates
  double *cFromFile = malloc(aarot->nbRot*aarot->nbAtoms*3*sizeof(double));
  for (int i=0; i< aarot->nbRot*aarot->nbAtoms*3; i++)
    {
      read = fscanf(in_file, "%lf", &cFromFile[i]);
      //printf("read %ld [%lf]\n", read, cFromFile[i]);    
    }
  for (int i=0; i<aarot->nbRot; i++) {
    double *c = malloc(aarot->nbAtoms*3*sizeof(double));
    int n=0;
    for (int j=0; j<aarot->nbAtoms; j++) {
      for (int k=0; k<3; k++) {
	c[n] = cFromFile[i*aarot->nbAtoms*3 + j*3 + k];
	n+=1;
      }
    }
    aarot->coords[i] = c;
  }

  /* printf("found rotamer %s %d %d\n", aaname, aarot->nbRot,  aarot->nbAtoms); */
  /* printf("    "); */
  /* for (int i=0; i<aarot->nbRot; i++) printf(" %f ", aarot->rotProbas[i]); */
  /* printf("\n    "); */
  /* for (int i=0; i<aarot->nbAtoms; i++) printf(" %d ", aarot->atypes[i]); */
  /* printf("\n    "); */
  /* for (int i=0; i<aarot->nbAtoms; i++) printf("%lf ", aarot->charges[i]); */
  /* printf("\n    %s\n", aarot->atnames); */
  /* printf("    %lf ... %lf\n", aarot->coords[0][0],  aarot->coords[aarot->nbRot-1][aarot->nbAtoms*3-1]); */
  return 1;
}

// create the table from rotamers.lib
int initialize_AASCRotTable_from_file(char *filename, int lastIndex)
{
  /* create _AASCRotTable for 20 standard amino acids */
  int i=0, retval;
  char aaname[] = "UNK                                                ";

  int nbRot = countRotamers(filename);

  // allocate table of rotamers structures
  //_AASCRotTable = malloc(nbCanAA*sizeof(struct _AASCRot));

  // load the rotamers
  FILE *in_file = fopen(filename, "r");
  for (i=0; i<nbRot; i+=1) {
    retval = getRotamerForAA(in_file, &_AASCRotTable[lastIndex+i], &aaname[0]);
    if (retval==1) {
      printf("loaded rotamer %6s at index %3d from %s\n", aaname, lastIndex+i, filename);
    } else {
      printf("Error: failed to load rotamer %s (index %d) from file %s\n", aaname, i, filename);
      exit(-1);
    }
  }
  //if (nbCanAA != i) printf("ERROR: intialize_AASCRotTable error while loading rotamers\n");
  //else printf("loaded rotamers for %d amino acid sidechains\n", nbCanAA);
  fclose(in_file);
  return lastIndex+nbRot;
}
