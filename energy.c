/*
** This is an implementation of model interactions between two amino acids
** as well within a single amino acid. This is a rather simple force-field.
**
** Copyright (c) 2004 - 2010 Alexei Podtelezhnikov
** Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
*/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<float.h>
#include<math.h>

//#include"canonicalAA.h"
#include"rotamers.h"
#include"error.h"
#include"params.h"
#include"vector.h"
#include"rotation.h"
#include"aadict.h"
#include"peptide.h"
#include"vdw.h"
#include"energy.h"

double centerX, centerY, centerZ, spacing;
int NX, NY, NZ;
double targetBest, currTargetEnergy;
const char atypes[MAX_ATOM_TYPES][3]= {
	  "C", "N", "OA", "HD", "SA", "A", "NA", "H", "HS", "NS",
	  "NX","OS","OX", "F",  "Mg", "MG", "P", "S", "SX", "Cl",
	  "CL","Ca","Mn","MN","Fe","FE","Zn","ZN","Br","BR","I"};

double *gridmapvalues[MAX_ATOM_TYPES];
int hasType[MAX_ATOM_TYPES];
double *emapvalues;
double *dmapvalues;

int transPtsCount;
double *Xpts;
double *Ypts;
double *Zpts;
double *ramaprob, *alaprob, *glyprob;

#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

//#define NaN 0.0 / 0.0
#define Erg(I,J)     erg[(I) * chain->NAA + (J)]
//#define Ergt(I,J)   ergt[(I - start) * chain->NAA + (J)]
#define Ergt(I,J)   ergt[(I) * chain->NAA + (J)]

/* contact map and associated stream */
#define Distb(I,J)     biasmap->distb[(I) * biasmap->NAA + (J)]


/***********************************************************/
/****                     CONSTANTS                     ****/
/***********************************************************/


/* CA-CA distance cutoff for Hbond interactions */
const double hbond_cutoff = 49.;

/* extern struct _ILE ILE; */
/* extern struct _LEU LEU; */
/* extern struct _PRO PRO; */
/* extern struct _VAL VAL; */
/* extern struct _PHE PHE; */
/* extern struct _TRP TRP; */
/* extern struct _TYR TYR; */
/* extern struct _ASP ASP; */
/* extern struct _GLU GLU; */
/* extern struct _ARG ARG; */
/* extern struct _HIS HIS; */
/* extern struct _LYS LYS; */
/* extern struct _SER SER; */
/* extern struct _THR THR; */
/* extern struct _CYS CYS; */
/* extern struct _MET MET; */
/* extern struct _ASN ASN; */
/* extern struct _GLN GLN; */

/***********************************************************/
/****       ENERGY MATRIX AND BIASMAP  OPERATIONS       ****/
/***********************************************************/


/* Initialize the energy matrix of a chain, that is,
   fill the already allocated matrix with the amino acid
   interaction energies.  In the diagonal the intraresidual
   energies go, and into (0,0), the global energy */
void energy_matrix_calculate(Chain *chain, Biasmap *biasmap, model_params *mod_params) {
	int i, j;



	//fprintf(stderr,"first row %g %g", chain->Erg(0, 0), chain->Erg(1, 0));
	/* (0,*) and (*,0) */

	for (i = 1; i < chain->NAA; i++){
		chain->Erg(0, i) = chain->Erg(i, 0) = 0.;
		fprintf(stderr,"%g ",chain->Erg(0,i));
	}
	//fprintf(stderr,"\n");

	/* (0,0) */
	chain->Erg(0, 0) = 0.0;
	if (mod_params->external_potential_type == 5){
		double ADenergies[chain->NAA-1];

		ADenergyNoClash(ADenergies, 1, chain->NAA-1, chain, NULL, mod_params, 0);

		for (i = 1; i < chain->NAA; i++) {
			chain->Erg(0, i) = ADenergies[i-1];
			fprintf(stderr," aaa %d %g \n",i, chain->Erg(0,i));
			//chain->Erg(0, i) = ADenergy(chain->aa + i, mod_params);
			chain->Erg(0, 0) += chain->Erg(0, i);
		}
		//free(ADenergies);
		//chain->Erg(0, 0) = global_energy(0,0,chain, NULL,biasmap, mod_params);

	}
	fprintf(stderr,"SS Energy ");
	chain->Erg(chain->NAA - 1, 0) = global_energy(0, 0,chain, NULL,biasmap, mod_params);

	if (mod_params->external_potential_type2 == 4)	chain->Erg(1, 0) = cyclic_energy((chain->aa) + 1, (chain->aa) + chain->NAA - 1, 0);
	/* diagonal */
	fprintf(stderr,"diag ");
	//fprintf(stderr,"ENERGY1 START\n");
	for (i = 1; i < chain->NAA; i++){
		chain->Erg(i, i) = energy1((chain->aa) + i, mod_params);
	//	fprintf(stderr,"%g ",chain->Erg(i,i));
	}
	//fprintf(stderr,"\n");
	//fprintf(stderr,"ENERGY1 END\n");
	for (i = 2; i < chain->NAA - 1; i++){
		chain->Erg(i, i) += ramabias((chain->aa)+i-1, (chain->aa) + i, (chain->aa) +i+1);
	}


	/* offdiagonal */
	fprintf(stderr,"offdiag ");
	for (i = 1; i < chain->NAA; i++){
		for (j = 1; j < i; j++){
			chain->Erg(i, j) = chain->Erg(j, i) = energy2(biasmap,(chain->aa) + i, (chain->aa) + j, mod_params);
	//	fprintf(stderr,"%g ",chain->Erg(i,j));
            
        }
	//fprintf(stderr,"\n");
    }
	if (mod_params->external_potential_type2 == 4)
		chain->Erg(1, chain->NAA-1) = chain->Erg(chain->NAA-1, 1) = energy2cyclic(biasmap,(chain->aa) + 1, (chain->aa) + chain->NAA-1, mod_params);
}

/* Calculate the total energy by adding up the energy matrix. */
double totenergy(Chain *chain)
{
	int i, j;
	double toten = 0.0;

	//fprintf(stderr, "tote... %g\n", chain->Erg(0,0));
	for (i = 0; i < chain->NAA; i++)
		for (j = 0; j <= i; j++) {
			toten += chain->Erg(i, j);
	//		fprintf(stderr, "%g ", chain->Erg(i,j));
		}
	//	fprintf(stderr,"\n");

	return toten;
}

/* Calculate the local energy by adding up the diagonal of the energy matrix. */
double locenergy(Chain *chain)
{
	int i;
	double toten = 0.;

	for (i = 1; i < chain->NAA; i++)
		toten += chain->Erg(i, i);

	return toten;
}

/* Return the external energy. */
double extenergy(Chain *chain)
{
	return chain->Erg(0, 0) + chain->Erg(1, 0);
}

/* Return the energy between first and last energy. */
double firstlastenergy(Chain *chain)
{
	return chain->Erg(1, chain->NAA - 1);
	//return chain->Erg(0, 0);
}


/* Print the energy matrix of a chain */
void energy_matrix_print(Chain *chain, Biasmap *biasmap, model_params *mod_params) {
	int i, j;
    for (i = 0; i < chain->NAA; i++)
		fprintf(stderr,"%g ",chain->Erg(0, i));
	fprintf(stderr,"\n");
    for (i = 0; i < chain->NAA; i++){
		for (j = 0; j <= i; j++){
			fprintf(stderr,"%g ",chain->Erg(i, j));
        }
	fprintf(stderr,"\n");
    }

}
/* Initialize matrix of Go-type contacts from a contact map file, unless
   NULL is given as the contact map filename (return a zero contact map).
   CAUTION!  The diagonal contains the secondary structure information
   of the amino acids.  If the amino acid one letter codes are in the
   diagonal, the contact map is meaningless. */
void biasmap_initialise(Chain *chain, Biasmap *biasmap, model_params *mod_params)
{
	int i, j, k;
	double val;

	FILE *fin = NULL;
	//fprintf(stderr,"Opening contact map file %s\n",mod_params->contact_map_file);
	fin = fopen(mod_params->contact_map_file,"r");
	int abort = 0;
	if (fin == NULL) { // no file opened beforehand
	    if (mod_params->contact_map_file==NULL) {
		fprintf(stderr, "ERROR: Invalid Go-type bias: %s.", mod_params->contact_map_file);
		fprintf(stderr, "  No contact map file specified.\n");
		abort=1;
	    } else if ((fin = fopen(mod_params->contact_map_file, "r")) == NULL && strcmp(mod_params->contact_map_file,"NULL")!=0) { // nonexisting file and not NULL given
		fprintf(stderr, "ERROR: Invalid Go-type bias: %s.", mod_params->contact_map_file);
		fprintf(stderr, "  Could not open file (missing?).\n");
		abort = 1;
	    }
	}
	if (abort) {
		stop("Contact map file has to be either an existing file or 'NULL' for no bias (see -p _B=... option).");
	}

	//if the biasmap has been already initialised (e.g. nested sampling), do not reread it
	if((biasmap)->distb != NULL) {
		if (fin) {
			fclose(fin);
			//fprintf(stderr,"Closing contact map file %s\n",mod_params->contact_map_file);
		}
		return; //when nested sampling - this is only to be done once
	}

	// allocate memory
	(biasmap)->distb = (double *) realloc((biasmap)->distb, chain->NAA * chain->NAA * sizeof(double));
	(biasmap)->NAA = chain->NAA;
	if ((biasmap)->distb == NULL) {
		stop("biasmap_initialise: Insufficient memory");
	}

	// if no biasmap, zero it
	if (strcmp(mod_params->contact_map_file,"NULL")==0) {
	    fprintf(stderr,"WARNING! Contact map file 'NULL' given, so no bias interactions will be calculated.\n");
	    for (i = 1; i < chain->NAA; i++) {
		for (j = 1; j < chain->NAA; j++) {
		    (biasmap)->distb[i * (biasmap)->NAA + j] = 0.0;
		}
	    }
	    if (fin) {
		fclose(fin);
		//fprintf(stderr,"Closing contact map file %s\n",mod_params->contact_map_file);
	    }
	    return;
	}

	// otherwise read in contact map
	for (i = 1; i < chain->NAA; i++)
		for (j = 1; j < chain->NAA; j++) {
			/* the fscanf consumes whitespaces, numbers,
			   decimal points, signs, but not most letters
			   so XOUZ is a good alphabet for symbolic maps */
			k = fscanf(fin, "%lf", &val); /* read in actual number 0, 1 or -1 */
			if (k == 0) { /* symbols used, not numbers */
				switch (fgetc(fin) & 0x3) { /* character code % 4 */
				case 0:	/* X */
					val = 1.0;
					break;
				case 1:	/* U */
					val = 1.0;
					break;
				case 2:	/* Z */
					val = -1.0;
					break;
				case 3:	/* O */
					val = 0.0;
					break;
				}
			} else if (k == EOF)
				goto out;

			/* use 0s in the diagonal, if in the diagonal of the contact map the amino acid letters are printed,
			otherwise they take a value according to the 1 letter code of the amino acid */
			//if (i==j) val=0.0;

			(biasmap)->distb[i * (biasmap)->NAA + j] = val;
		}
      out:
	if (fin) {
		fclose(fin);
		//fprintf(stderr,"Closing contact map file %s\n",mod_params->contact_map_file);
	}
	fprintf(stderr, "Go-type bias: %dx%d\n", i - 1, chain->NAA - 1);

	/* use 0s in the diagonal, if in the diagonal of the contact map the amino acid letters are printed,
	otherwise they take a value according to the 1 letter code of the amino acid */
	//fprintf(stderr,"WARNING:  If the one-letter amino acid codes are in the diagonal, the bias potential's energy contribution will be meaningless!  Check your contact map.\n");

	/* special treatment for glycine */
	if (mod_params->prt != 0.0)
		for (i = 1; i < chain->NAA; i++)
			if (chain->aa[i].id == 'G') {
				//fprintf(stderr,"Gly %d (%d)",i,biasmap->NAA);
				for (j = 1; j < chain->NAA; j++)
					if (abs(i - j) > 1) {
						//fprintf(stderr," %d",j);
						(biasmap)->distb[i * (biasmap)->NAA + j] = (biasmap)->distb[j * (biasmap)->NAA + i] = 0.0;
					}
				//fprintf(stderr,"\n");
			}

	/* check if the biasing matrix is symmetric */
	/* possible reasons for this are:
		the contact map is lopsided
		the contact contains spaces, which are skipped by fscanf
		the one-letter amino acid codes are used in the diagonal, and there was a problem parsing N or I (on francesca) */
	int lopsided = 0;
	for (i = 1; i < chain->NAA; i++) {
		for (j = 1; j < i; j++) {
			if ((biasmap)->distb[i * (biasmap)->NAA + j] != (biasmap)->distb[j * (biasmap)->NAA + i]) {
				fprintf(stderr, "Lopsided bias: %d %d valued %g %g\n", i, j, (biasmap)->distb[i * (biasmap)->NAA + j], (biasmap)->distb[j * (biasmap)->NAA + i]);
				lopsided = 1;
			}
		}
	}
	if (lopsided) stop("Lopsided bias map!");

}

/* finalize matrix of Go-type contacts from contact map file
   close input file  */
void biasmap_finalise(Biasmap *biasmap){

    if(biasmap){   
	if (biasmap->distb) free(biasmap->distb);
	free(biasmap); 
    }

}


/*make energy grid map smoother*/
double lower_gridenergy(double E) {
	//return E;
	if (E > 2.718) {
		//return log10f(E) + 9;
		return log(E) + 1.718;
	}
	//if (E > 10.) {
		//return log10f(E) + 9;
	//	return log(E-9) + 10.;
	//}
	return E;
}

/*init the size and center and spacing of AD gridbox*/
void gridbox_initialise(simulation_params *sim_params) {
	FILE *gridmap = NULL;
	char line[256];
	sprintf(line, "%srigidReceptor.C.map", sim_params->target_folder);
	gridmap = fopen(line, "r");
	int i = 0, j = 0;
	while (fgets(line, sizeof(line), gridmap)) {
		if (i < 3) {
			//printf("%s", line);
			i++;
			continue;
		}
		char * pch;
		pch = strtok(line, " ");
		j = 0;
		while (pch != NULL)
		{
			if (i == 3 && j == 1) {
				spacing = atof(pch);
			}
			else if (i == 4 && j == 1) {
				NX = atoi(pch) + 1;
			}
			else if (i == 4 && j == 2) {
				NY = atoi(pch) + 1;
			}
			else if (i == 4 && j == 3) {
				NZ = atoi(pch) + 1;
			}
			else if (i == 5 && j == 1) {
				centerX = atof(pch);
			}
			else if (i == 5 && j == 2) {
				centerY = atof(pch);
			}
			else if (i == 5 && j == 3) {
				centerZ = atof(pch);
			}
			pch = strtok(NULL, " ");
			j++;
		}
		i++;
		if (i > 5) break;
	}
	printf("grid box initialise success %i %i %i \n", NX, NY, NZ);
	fclose(gridmap);
}


/* initialise the grip maps from the file, elements are 0:C, 1:N, 2:O, 3:H, 4:S, 5:CA, 6:NA*/
void gridmap_initialise(char *filename, int atype) {
	FILE *gridmap_file = NULL;
	gridmap_file = fopen(filename, "r");
	if (gridmap_file == NULL) {
		stop("Missing gridmap_file.map file.");
	}
	char line[256];
	int i = 0;
	double *curr_gridmap_values = malloc(NX*NY*NZ * sizeof(double));
	while (fgets(line, sizeof(line), gridmap_file)) {
		if (i < 6) {
			i++;
			continue;
		}
		curr_gridmap_values[i - 6] = lower_gridenergy(atof(line));
		i++;
	}
	fclose(gridmap_file);
	if (atype==-1) {
	  emapvalues = curr_gridmap_values;
	} else {
	  if (atype==-2) {
	    dmapvalues = curr_gridmap_values;
	  } else {
	    gridmapvalues[atype] = curr_gridmap_values;
	  }
	}
}

/* initialise the tranpoints from the file, if no transpoints found, add the box center */
void transpts_initialise(simulation_params *sim_params) {
	FILE *transpts_file = NULL;
	char filename[255];
	sprintf(filename, "%stranspoints", sim_params->target_folder);
	
	transpts_file = fopen(filename, "r");
	if (transpts_file == NULL) {
	  char msg[512];
	  sprintf(msg, "could not open file %s\n", filename);
	  stop(msg);
	}
	char line[256];
	int i = 0, j = 0;	
	fgets(line, sizeof(line), transpts_file);
	transPtsCount = atoi(line);
	if (transPtsCount == 0 || transpts_file == NULL) {
		printf("no transpoints found \n");

		Xpts = malloc(1 * sizeof(double));
		Ypts = malloc(1 * sizeof(double));
		Zpts = malloc(1 * sizeof(double));
		Xpts[0] = centerX;
		Ypts[0] = centerY;
		Zpts[0] = centerZ;
		transPtsCount = 1;
		fclose(transpts_file);
		return;
	}
	Xpts = malloc(transPtsCount * sizeof(double));
	Ypts = malloc(transPtsCount * sizeof(double));
	Zpts = malloc(transPtsCount * sizeof(double));


	while (fgets(line, sizeof(line), transpts_file)) {
		char * pch;
		pch = strtok(line, " ");
		j = 0;
		while (pch != NULL)
		{
			if (j == 0) {
				Xpts[i] = atof(pch);
			}
			else if (j == 1) {
				Ypts[i] = atof(pch);
			}
			else if (j == 2) {
				Zpts[i] = atof(pch);
			}
			pch = strtok(NULL, " ");
			j++;
		}
		i++;
	}

	printf("transpoints initialise success with %i transpoints \n", transPtsCount);
	fclose(transpts_file);

}

/*initialise the ramachandra probability from ramaprob.data file*/
void ramaprob_initialise(char *folder) {
	FILE *ramaprob_file = NULL;
	char buffer[254];
	strcpy(buffer, folder);
	strcat(buffer, "ramaprob.data");
	ramaprob_file = fopen(buffer, "r");
	if (ramaprob_file == NULL) {
		stop("Missing ramaprob.data file.");
	}
	char line[256];
	int i = 0, j = 0;
	ramaprob = malloc(32400 * sizeof(double));
	alaprob = malloc(32400 * sizeof(double));	
	glyprob = malloc(32400 * sizeof(double));

	while (fgets(line, sizeof(line), ramaprob_file)) {
		char * pch;
		pch = strtok(line, " ");
		j = 0;
		while (pch != NULL)
		{
			if (j == 1) {
				ramaprob[i] = atof(pch);
			}
			else if (j == 2) {
				alaprob[i] = atof(pch);
			}
			else if (j == 3) {
				glyprob[i] = atof(pch);
			}
			pch = strtok(NULL, " ");
			j++;
		}
		i++;
	}

	printf("ramaprob initialise success \n");
	fclose(ramaprob_file);

}

/***********************************************************/
/****               ENERGY  CONTRIBUTIONS               ****/
/***********************************************************/


/* Low level routine to calculate the intensity of a linearly decaying function
   from a maximum value to zero beyond a cutoff function:

   strength_____.
                .\
   	        . \
   0 . . . .  . . .\._________
	        |<->|-decay width
		|
		|-cutoff

 */
inline double linear_decay(double distance,       /* distance of the 2 atoms */
		           double contact_cutoff, /* upto which sum of their contact radii */
			   //double strength, /* the maximum value taken below the cutoff */
			   double decay_width /* the width of the linear decay */ ) {

	if (distance > contact_cutoff + decay_width) return 0.0;
	if (distance < contact_cutoff) return 1.0;
	return 1.0 - (distance - contact_cutoff) / decay_width;
}


/* Energy contribution of proline phi dihedral angles
   E = 30.0 (RT) * (phi - phi_0)^2
   phi_0: equilibrium value, -PI/3 */
/* proline B restraints its phi angle with preceding A */
double proline(AA *a, AA *b)
{
	double phi;

	if (b->id != 'P')
		return 0.0;

	phi = dihedral_4(a->c, b->n, b->ca, b->c) + M_PI_3;

	return 30. * phi * phi;	/* RMSD by Ho et al. (2004) */
}


/* Bending energy contribution (internal stress)
   E = 150 (RT) * ( angle(n,ca,c) - arccos(-1/3) )^2 */
double stress(AA *a, model_params *mod_params)
{
	vector nca, cac;	/* N-Ca and Ca-C bonds */
	double beta, erg = 0.0;

	subtract(nca, a->ca, a->n);
	subtract(cac, a->c, a->ca);
	
	/* ground-state angle is 180 - 111 = 69 (EH2001) */
	beta = angle(nca, cac) - mod_params->stress_angle; //1.20427718387608740808;
	erg += mod_params->stress_k * beta * beta;	/* softer than Engh and Huber (2001) */
	//if (a->id == 1) erg = erg * 200;
	return erg;
}


/*translate rama probability into energy*/
double ramabias(AA *prevaa, AA *a, AA *nextaa)
{
	//return 0.0;
	double phi = 0.0, psi = 0.0;
	double energy = 0.0;
	int segphi = 0;
	int segpsi = 0;
	
	phi = dihedral_rama(prevaa->c, a->n, a->ca, a->c, 1.46);
	psi = dihedral_rama(a->n, a->ca, a->c, nextaa->n, 1.53);

	// MS added inverting phi psi angles to get probabilities for D-amino acids
	// doi: 10.1093/protein/gzy016
	if ((a->etc & LEV)==0) { // D amino acid 
	  phi = -phi;
	  psi = -psi;
	}

	//double anglephi = phi * M_180_PI;
	//double anglepsi = psi * M_180_PI;
	//int sig = (anglephi>0)? 1:-1;
	//int segphi = 2 * int(anglephi/2)+sig;
	//sig = (anglepsi>0)? 1:-1;
	//int segpsi = 2 * int(anglepsi/2)+sig;

	segphi = (int) ((phi + M_PI)/(2*M_PI_180));
	segpsi = (int) ((psi + M_PI)/(2*M_PI_180));

	int ind = segphi * 180 + segpsi;

	switch (a->id) {
		case 'A':
			energy = alaprob[ind];
			break;
		case 'G':
			energy = glyprob[ind];
			break;
		default:
			energy = ramaprob[ind];
			break;
	}
	
	//fprintf(stderr,"aaa num %d id %c ind %d energy %g segphi %d segpsi %d phi %g psi %g \n",a->num, a->id, ind, -energy, segphi, segpsi, anglephi, anglepsi);
	//energy = energy < -3.91 ? energy + 3.91 : 0.0;
	return -energy;	// RT * 0.59219 = kcal/mol
}


/* Csilla: energy contribution of Go-type potential
   E = kappa * C_ij * r_ij^2   for |i-j|>1
   E = eta * cos(gamma_ij)     for i-j=1
   kappa: force constant (different for alpha-helix and beta-sheet)
   C_ij:
   r_ij: C_beta distances
   eta: force constant (different for alpha-helix(positive) and beta-sheet(negative))
   gamma_ij: C_beta,i-C_alpha,i-C_alpha,j-C_beta,j dihedral angle
   a, b: i-th and j-th amino acids
*/
/* Go-type biasing potentials that stabilize alpha-helices and beta-sheets */
double bias(Biasmap *biasmap, AA *a, AA *b, model_params *mod_params)
{
	
	int i = a->num, j = b->num;
	double bb, dst2;
	double r = 0;

    double rs = 2.15; /* rs */
    double bs = -0.25; /*kappa s */

	if (biasmap->distb == NULL) {
		stop("Contact matrix is not initialised.  Perhaps missing contact map?");
	}

	/* works for multi-chain proteins */
	switch ( max(abs(i - j), 1000 * abs(a->chainid - b->chainid)) ) {
	case 0:
		return 0.0;
	case 1:
		if (Distb(i, j) < 0.) {
			bb = mod_params->bias_eta_beta;	/* beta-strand twist */   /* eta_beta */
		} else {
			bb = mod_params->bias_eta_alpha;	/* alpha-helix twist */   /* eta_alpha */
		}
		break;
	case 3:
		if(Distb(i, j) < 0){
          /*Cys bond */
		  bb = bs;
		  r = rs;
		}
		else if (Distb(i, i) > 0. && Distb(j, j) > 0.){
			bb = mod_params->bias_kappa_alpha_3;	/* alpha-helix elasticity */   /* kappa_alpha */
		    r = mod_params->bias_r_alpha;
	        }
		else{
			if (rand() < Distb(i, j)*RAND_MAX) {
				bb = mod_params->bias_kappa_beta;   /* kappa_beta */
			} else {
				bb = 0;
			}
		    r = mod_params->bias_r_beta; 
		    }
		break;
	case 4:
		if(Distb(i, j) < 0){
          /*Cys bond */
		  bb = bs;
		  r = rs;
		}
		else if (Distb(i, i) > 0. && Distb(j, j) > 0.){
			bb = mod_params->bias_kappa_alpha_4;	/* alpha-helix elasticity */   /* kappa_alpha */
			r = mod_params->bias_r_alpha;
		}
		else{
			if (rand() < Distb(i, j)*RAND_MAX) {
				bb = mod_params->bias_kappa_beta;   /* kappa_beta */
			} else {
				bb = 0;
			}
			r = mod_params->bias_r_beta;
	    }
		break;
	default:
		if(Distb(i, j) < 0){
        /*Cys bond */
		  bb = bs;
		  r = rs;
		}
		else{
		  bb = mod_params->bias_kappa_beta;	/* beta-sheet elasticity */  /* kappa_beta */
		  r = mod_params->bias_r_beta;
	    }
	}


	/* works for multi-chain proteins */
	if ( max(abs(i - j), 1000 * abs(a->chainid - b->chainid)) > 1 ) {
		vector x, y;
		lincomb(x, 1. - mod_params->prt, a->ca, mod_params->prt, a->cb);
		lincomb(y, 1. - mod_params->prt, b->ca, mod_params->prt, b->cb);
		dst2 = distance(x, y);
		dst2 = dst2 - 2*sqrt(dst2)*r + r*r;
		
	} else { /* neighbouring amino acids */
		vector x, y, z;
		if (i > j) {	/* reorder */
			AA *temp;
			temp = a;
			a = b;
			b = temp;
		}

		/* either baab or naac pseudo-dihedral is suitable */
		subtract(x, a->ca, a->n);
		subtract(y, b->ca, a->ca);
		subtract(z, b->c, b->ca);
		//set helix eta alpha phaes shift
		if (Distb(i, j) > 0.) dst2 = -phasindihedral(x,y,z, 0.13917, 0.99); 
		else dst2 = -cosdihedral(x, y, z); 
		//dst2 = -phasindihedral(x,y,z, 0.985,0.174); //strand
	}

	return Distb(i, j) * bb * dst2;
}

/* Energy contribution of hydrogen bonds
  E = num_hb * hbs
  num_hb: number of H-bonds between the two amino acids
  hbs: H-bond strength (default: 2.2 (RT))
*/
double hbond(Biasmap *biasmap, AA *a, AA *b, model_params *mod_params)
{

	//return -hbs * (hdonor(a, b) + hdonor(b, a));
	
	double fact = 1.0;
	int i = a->num; int j = b->num;
	if (abs(i-j) == 4)
		fact = 1.0;
	else if (abs(i-j) == 3)
		fact = 1.0;
	else
		fact = 1.0;
	// Gary hack to put more internal Hbond with cyclic CYS-CYS
	if(mod_params->Sbond_strength != 0) fact = 1.;

//	if( Distb( i, i ) * Distb( j, j ) == 0) fact /= 3.0;	
	//fprintf(stderr,"hbond %d %c",a->num,a->id);
	//fprintf(stderr," N: %g",a->n[0]);
	//fprintf(stderr," %g",a->n[1]);
	//fprintf(stderr," %g",a->n[2]);
	//fprintf(stderr," H: %g",a->h[0]);
	//fprintf(stderr," %g",a->h[1]);
	//fprintf(stderr," %g",a->h[2]);
	//fprintf(stderr," , %d %c",b->num,b->id);
	//fprintf(stderr," O: %g",b->o[0]);
	//fprintf(stderr," %g",b->o[1]);
	//fprintf(stderr," %g",b->o[2]);
	//fprintf(stderr," C: %g",b->c[0]);
	//fprintf(stderr," %g",b->c[1]);
	//fprintf(stderr," %g\n",b->c[2]);

	if (a->id == 'P') {
		if (b->id == 'P') // no a->H, no b->H
			return 0.0;
		else // no a->H
			return fact*-mod_params->hbs * hstrength(b->n, b->h, a->o, a->c, mod_params);
	}
	else {
		if (b->id == 'P') // no b->H
			return fact*-mod_params->hbs * hstrength(a->n, a->h, b->o, b->c, mod_params);
		else
			return fact*-mod_params->hbs * (hstrength(a->n, a->h, b->o, b->c, mod_params) + hstrength(b->n, b->h, a->o, a->c, mod_params));
	}
}

/* The scaling factor of hydrophobic interaction between a and b
   R. Srinivasan, ProtSFG 22, 81--99 (1995)
   R. Srinivasan et al., PNAS 96(25), 14258--14263 (1999)
   2 if both amino acids are hydrophobic,
   1 if one amino acid is hydrophobic, and the other one is amphipathic,
   0 otherwise */
inline int hydrophobic_interaction_intensity(AA *a, AA *b, model_params *mod_params) {

	/* works for multi-chain proteins */
	if ( ( abs(a->num - b->num) < mod_params->hydrophobic_min_separation ) && ( a->chainid == b->chainid ) ) return 0; 

	/* hydrophobic -- hydrophobic */
	if ( (a->etc & HYDROPHOBIC) && (b->etc & HYDROPHOBIC) ) return 2;
	if ( (a->etc & HYDROPHOBIC) && (b->etc & AMPHIPATHIC) ) return 1;
	if ( (a->etc & AMPHIPATHIC) && (b->etc & HYDROPHOBIC) ) return 1;

	/* hydrophobic -- polar */
	/* this accounts for the fewer polar -- polar interactions */
	if (a->etc & HYDROPHOBIC) return -2;
	if (b->etc & HYDROPHOBIC) return -2;
	
	return 0;
}

/* Low level routine to calculate 2 atoms hydrophobic energy contribution
   E = -k_h * f(d(a,b)) * intensity,     if a and b are in contact
   k_h is the hydrophobicity parameter (in RT)
   d(a,b) is the distance between the respective atoms of a and b
   f(d(a,b)) = d(a,b)^-1 */
/* inline double hydrophobic_low_recip(double distance,
	       double contact_cutoff,
	       model_params *mod_params) {

	if (distance > mod_params->hydrophobic_max_cutoff) return 0.0;
	if (distance < mod_params->hydrophobic_min_cutoff) distance = mod_params->hydrophobic_min_cutoff;
//fprintf(stderr,"   calc %g\n",(1.0/distance) - mod_params->hydrophobic_max_Eshift);
	return (1.0/distance) - mod_params->hydrophobic_max_Eshift;
} */

/* Low level routine to calculate 2 atoms hydrophobic energy contribution using a spline potential
   E = -k_h * f(d(a,b)) * intensity,     if a and b are in contact
   k_h is the hydrophobicity parameter (in RT)
   d(a,b) is the distance between the respective atoms of a and b
   f(d(a,b)) is a beta spline going from 1 to 0 between min_cutoff and max_cutoff */
/* inline double hydrophobic_low_spline(double distance,
	       double contact_cutoff,
	       model_params *mod_params) {

	if (distance > mod_params->hydrophobic_r + mod_params->hydrophobic_half_delta) return 0.0;
	if (distance < mod_params->hydrophobic_r - mod_params->hydrophobic_half_delta) distance = mod_params->hydrophobic_r - mod_params->hydrophobic_half_delta;
	return (distance - mod_params->hydrophobic_r - mod_params->hydrophobic_half_delta)
		  * (distance - mod_params->hydrophobic_r - mod_params->hydrophobic_half_delta)
		  * (mod_params->hydrophobic_r - 2.0 * mod_params->hydrophobic_half_delta - distance)
		  / (4.0 * mod_params->hydrophobic_half_delta * mod_params->hydrophobic_half_delta * mod_params->hydrophobic_half_delta) ;
} */

/* Low level routine to calculate 2 atoms hydrophobic energy contribution */
/* R. Srinivasan et al., PNAS 96(25), 14258--14263 (1999) */
inline double hydrophobic_low(double distance /* distance of the 2 atoms */,
		       double contact_cutoff /* sum of their contact radii */,
		       model_params *mod_params) {

	if (distance > contact_cutoff + mod_params->hydrophobic_cutoff_range) return 0.0;
	if (distance < contact_cutoff) return 1.0;
	return (distance - contact_cutoff) / mod_params->hydrophobic_cutoff_range;
}

/* Energy contribution of hydrophobic interaction between a and b
   E = -k_h * f(d(a,b)) * intensity,     if a and b are in contact
   k_h is the hydrophobicity parameter (in RT)
   d(a,b) is the distance between the respective atoms of a and b
   intensity depends on the hydrophobicity of the amino acids (see hydrophobic_intensity) */
/* Default: f(d(a,b)) = 1 */
double hydrophobic(Biasmap *biasmap, AA *a, AA *b, model_params *mod_params) {

	/* return 0 if 0 */
	if (mod_params->kauzmann_param == 0.0) return 0.0;

	int intensity;
	if ((intensity = hydrophobic_interaction_intensity(a,b,mod_params)) == 0) return 0.0;

	/* calc hydrophobic contact radii */
	double r_cb_a=0, r_g_a=0, r_g2_a=0;
	double r_cb_b=0, r_g_b=0, r_g2_b=0;
	if (a->etc & CB_) r_cb_a = hydrophobic_contact_radius(a->id, CB_, mod_params->sidechain_properties);
	if (a->etc & G__) r_g_a  = hydrophobic_contact_radius(a->id, G__, mod_params->sidechain_properties);
	if (a->etc & G2_) r_g2_a = hydrophobic_contact_radius(a->id, G2_, mod_params->sidechain_properties);
	if (b->etc & CB_) r_cb_b = hydrophobic_contact_radius(b->id, CB_, mod_params->sidechain_properties);
	if (b->etc & G__) r_g_b  = hydrophobic_contact_radius(b->id, G__, mod_params->sidechain_properties);
	if (b->etc & G2_) r_g2_b = hydrophobic_contact_radius(b->id, G2_, mod_params->sidechain_properties);

	double energy = 0;
	if (mod_params->use_gamma_atoms != NO_GAMMA) {
		/* all side chain contributions */
		/* peptide should have been fixed by now, so no missing coordinates */
		if ( (a->etc & CB_) && (b->etc & CB_) && r_cb_a > 0. && r_cb_b > 0. ) { /* CB -- CB */
			energy += hydrophobic_low(sqrt(distance(a->cb, b->cb)), r_cb_a + r_cb_b, mod_params);
		}
		if ( (a->etc & CB_) && (b->etc & G__) && r_cb_a > 0. && r_g_b > 0. ) { /* CB -- G1 */
			energy += hydrophobic_low(sqrt(distance(a->cb, b->g )), r_cb_a + r_g_b , mod_params);
		}
		if ( (a->etc & CB_) && (b->etc & G2_) && r_cb_a > 0. && r_g2_b > 0. ) { /* CB -- G2 */
			energy += hydrophobic_low(sqrt(distance(a->cb, b->g2)), r_cb_a + r_g2_b, mod_params);
		}
		if ( (a->etc & G__) && (b->etc & CB_) && r_g_a > 0. && r_cb_b > 0. ) { /* G1 -- CB */
			energy += hydrophobic_low(sqrt(distance(a->g , b->cb)), r_g_a  + r_cb_b, mod_params);
		}
		if ( (a->etc & G__) && (b->etc & G__) && r_g_a > 0. && r_g_b > 0. ) { /* G1 -- G1 */
			energy += hydrophobic_low(sqrt(distance(a->g , b->g )), r_g_a  + r_g_b , mod_params);
		}
		if ( (a->etc & G__) && (b->etc & G2_) && r_g_a > 0. && r_g2_b > 0. ) { /* G1 -- G2 */
			energy += hydrophobic_low(sqrt(distance(a->g , b->g2)), r_g_a  + r_g2_b, mod_params);
		}
		if ( (a->etc & G2_) && (b->etc & CB_) && r_g2_a > 0. && r_cb_b > 0. ) { /* G2 -- CB */
			energy += hydrophobic_low(sqrt(distance(a->g2, b->cb)), r_g2_a + r_cb_b, mod_params);
		}
		if ( (a->etc & G2_) && (b->etc & G__) && r_g2_a > 0. && r_g_b > 0. ) { /* G2 -- G1 */
			energy += hydrophobic_low(sqrt(distance(a->g2, b->g )), r_g2_a + r_g_b , mod_params);
		}
		if ( (a->etc & G2_) && (b->etc & G2_) && r_g2_a > 0. && r_g2_b > 0. ) { /* G2 -- G2 */
			energy += hydrophobic_low(sqrt(distance(a->g2, b->g2)), r_g2_a + r_g2_b, mod_params);
		}
		//fprintf(stderr,"%d %d %f ",a->num,b->num, energy);
		return -mod_params->kauzmann_param * energy * (double) intensity;
		
	} else {
		/* OLD VERSION */
		if (contact(a,b, mod_params)) {
			return (-mod_params->kauzmann_param * ((double) intensity));
		} else {
			return 0.0;
		}
	}
	
}


double sidechain_hbond(Biasmap *biasmap, AA *a, AA *b, model_params *mod_params)
{

	double intensity;
	double erg = 0.0;
	double hbond_distance;
	double cos_goc_angle;

	/* only from i, i+X */
	/* works for multi-chain proteins */
	if ( ( abs(a->num - b->num) < mod_params->sidechain_hbond_min_separation ) && ( a->chainid == b->chainid ) ) return 0.0; 

	/* side chain donor - backbone acceptor */
	if (a->etc &G__ && b->etc &C__ && hbond_donor(a->id,G__, mod_params->sidechain_properties)) {
	   // check G--C distance
	   hbond_distance = sqrt(distance(a->g,b->c));
	   if (hbond_distance < sidechain_hbond_donor_radius(a->id,mod_params->sidechain_properties) + BACKBONE_ACCEPTOR_RADIUS + mod_params->sidechain_hbond_decay_width) {
//fprintf(stderr,"%c%d %c%d G-C distance %g ,",a->id,a->num,b->id,b->num,hbond_distance);
	      // check C--O--G angle
	      cos_goc_angle = cosangle(a->g,b->o,b->c);
//fprintf(stderr,"cos_goc_angle %g\n",cos_goc_angle);
	      if (cos_goc_angle < mod_params->sidechain_hbond_angle_cutoff) {
		 if ( (intensity = linear_decay(hbond_distance, sidechain_hbond_donor_radius(a->id,mod_params->sidechain_properties) + BACKBONE_ACCEPTOR_RADIUS, mod_params->sidechain_hbond_decay_width )) > 0.0 ) {
//		      fprintf(stderr,"found hbond %c %d G__ >> %c %d C__ %g\n",a->id,a->num,b->id,b->num,distance(a->g,b->c));
		    erg += -mod_params->sidechain_hbond_strength_s2b * intensity;
		 }
	      }
	   }
	}
	if (b->etc &G__ && a->etc &C__ && hbond_donor(b->id,G__, mod_params->sidechain_properties)) {
	   // check G--C distance
	   hbond_distance = sqrt(distance(b->g,a->c));
	   if (hbond_distance < sidechain_hbond_donor_radius(b->id,mod_params->sidechain_properties) + BACKBONE_ACCEPTOR_RADIUS + mod_params->sidechain_hbond_decay_width) {
//fprintf(stderr,"%c%d %c%d G-C distance %g ,",a->id,a->num,b->id,b->num,hbond_distance);
	      // check C--O--G angle
	      cos_goc_angle = cosangle(b->g,a->o,a->c);
//fprintf(stderr,"cos_goc_angle %g\n",cos_goc_angle);
	      if (cos_goc_angle < mod_params->sidechain_hbond_angle_cutoff) {
		 if ( (intensity = linear_decay(hbond_distance, sidechain_hbond_donor_radius(b->id,mod_params->sidechain_properties) + BACKBONE_ACCEPTOR_RADIUS, mod_params->sidechain_hbond_decay_width )) > 0.0 ) {
//		      fprintf(stderr,"found hbond %c %d G__ >> %c %d C__ %g\n",a->id,a->num,b->id,b->num,distance(a->g,b->c));
		    erg += -mod_params->sidechain_hbond_strength_s2b * intensity;
		 }
	      }
	   }
	}

	/* side chain acceptor - backbone donor */
	if (a->etc &G__ && b->etc &N__ && b->etc &H__ && hbond_acceptor(a->id,G__, mod_params->sidechain_properties)) {
	   // check G--N distance
	   hbond_distance = sqrt(distance(a->g,b->n));
//fprintf(stderr,"check %c %d %c %d G-N distance %g (%g) %g \n",a->id,a->num,b->id,b->num,hbond_distance,sidechain_hbond_acceptor_radius(a->id,mod_params->sidechain_properties) + BACKBONE_DONOR_RADIUS + sidechain_hbond_decay_width,sqrt(distance(a->g,b->h)));
	   if (hbond_distance < sidechain_hbond_acceptor_radius(a->id,mod_params->sidechain_properties) + BACKBONE_DONOR_RADIUS + mod_params->sidechain_hbond_decay_width) {
//fprintf(stderr,"found %c %d %c %d G-N distance %g (%g) %g \n",a->id,a->num,b->id,b->num,hbond_distance,sidechain_hbond_acceptor_radius(a->id,mod_params->sidechain_properties) + BACKBONE_DONOR_RADIUS + sidechain_hbond_decay_width,sqrt(distance(a->g,b->h)));
	      // check C--H--N angle
	      cos_goc_angle = cosangle(a->g,b->h,b->n);
//fprintf(stderr,"cos_ghn_angle %g\n",cos_goc_angle);
	      if (cos_goc_angle < mod_params->sidechain_hbond_angle_cutoff) {
		 if ( (intensity = linear_decay(hbond_distance, sidechain_hbond_acceptor_radius(a->id,mod_params->sidechain_properties) + BACKBONE_DONOR_RADIUS, mod_params->sidechain_hbond_decay_width )) > 0.0 ) {
//		      fprintf(stderr,"found hbond %c %d G__ >> %c %d N__ %g\n",a->id,a->num,b->id,b->num,distance(a->g,b->n));
		    erg += -mod_params->sidechain_hbond_strength_b2s * intensity;
		 }
	      }
	   }
	}
	if (b->etc &G__ && a->etc &N__ && a->etc &H__ && hbond_acceptor(b->id,G__, mod_params->sidechain_properties)) {
	   // check G--N distance
	   hbond_distance = sqrt(distance(b->g,a->n));
//fprintf(stderr,"check %c %d %c %d G-N distance %g (%g) %g\n",b->id,b->num,a->id,a->num,hbond_distance,sidechain_hbond_acceptor_radius(b->id,mod_params->sidechain_properties) + BACKBONE_DONOR_RADIUS + sidechain_hbond_decay_width,sqrt(distance(b->g,a->h)));
	   if (hbond_distance < sidechain_hbond_acceptor_radius(b->id,mod_params->sidechain_properties) + BACKBONE_DONOR_RADIUS + mod_params->sidechain_hbond_decay_width) {
//fprintf(stderr,"found %c %d %c %d G-N distance %g (%g) %g\n",b->id,b->num,a->id,a->num,hbond_distance,sidechain_hbond_acceptor_radius(b->id,mod_params->sidechain_properties) + BACKBONE_DONOR_RADIUS + sidechain_hbond_decay_width,sqrt(distance(b->g,a->h)));
	      // check C--H--N angle
	      cos_goc_angle = cosangle(b->g,a->h,a->n);
//fprintf(stderr,"cos_ghn_angle %g\n",cos_goc_angle);
	      if (cos_goc_angle < mod_params->sidechain_hbond_angle_cutoff) {
		 if ( (intensity = linear_decay(hbond_distance, sidechain_hbond_acceptor_radius(b->id,mod_params->sidechain_properties) + BACKBONE_DONOR_RADIUS, mod_params->sidechain_hbond_decay_width )) > 0.0 ) {
//		      fprintf(stderr,"found hbond %c %d G__ >> %c %d N__ %g\n",a->id,a->num,b->id,b->num,distance(a->g,b->n));
		    erg += -mod_params->sidechain_hbond_strength_b2s * intensity;
		 }
	      }
	   }
	}

	/* side chain donor - side chain acceptor */
	if (a->etc &G__ && b->etc &G__) {
	   if (hbond_donor(a->id,G__, mod_params->sidechain_properties) && hbond_acceptor(b->id,G__, mod_params->sidechain_properties)) {
	      // check G-G' distance
	      hbond_distance = sqrt(distance(a->g,b->g));
	      if (hbond_distance < sidechain_hbond_donor_radius(a->id,mod_params->sidechain_properties) + sidechain_hbond_acceptor_radius(b->id,mod_params->sidechain_properties) + mod_params->sidechain_hbond_decay_width) {
//fprintf(stderr,"%c%d %c%d G-G distance %g ,",a->id,a->num,b->id,b->num,hbond_distance);
		 // check B-G,G'-B' angle
		 vector x, z;
		 subtract(x, a->g, a->cb);
		 subtract(z, b->g, b->cb);
//fprintf(stderr,"cos_gbbg = %g, gb.bg = %g\n",cosine(x,z),dotprod(x,z));
//		 if (cosine(x,z)<mod_params->sidechain_hbond_angle_cutoff) { // && dotprod(x,z)<0) {
		    if ( (intensity = linear_decay(hbond_distance, sidechain_hbond_donor_radius(a->id,mod_params->sidechain_properties) + sidechain_hbond_acceptor_radius(b->id,mod_params->sidechain_properties), mod_params->sidechain_hbond_decay_width )) > 0.0 ) {
//	         fprintf(stderr,"found hbond %c %d G__ >> %c %d G__ %g\n",a->id,a->num,b->id,b->num,distance(a->g,b->g));
		       erg += -mod_params->sidechain_hbond_strength_s2s * intensity;
		    }
//		 }
	      }
	   }
	   if (hbond_donor(b->id,G__, mod_params->sidechain_properties) && hbond_acceptor(a->id,G__, mod_params->sidechain_properties)) {
	      // check G-G' distance
	      hbond_distance = sqrt(distance(a->g,b->g));
	      if (hbond_distance < sidechain_hbond_donor_radius(b->id,mod_params->sidechain_properties) + sidechain_hbond_acceptor_radius(a->id,mod_params->sidechain_properties) + mod_params->sidechain_hbond_decay_width) {
//fprintf(stderr,"%c%d %c%d G-G distance %g ,",a->id,a->num,b->id,b->num,hbond_distance);
		 // check B-G,G'-B' angle
		 vector x, z;
		 subtract(x, a->g, a->cb);
		 subtract(z, b->g, b->cb);
//fprintf(stderr,"cos_gbbg = %g, gb.bg = %g\n",cosine(x,z),dotprod(x,z));
//		 if (cosine(x,z)<mod_params->sidechain_hbond_angle_cutoff) { // && dotprod(x,z)<0) {
		    if ( (intensity = linear_decay(hbond_distance, sidechain_hbond_donor_radius(b->id,mod_params->sidechain_properties) + sidechain_hbond_acceptor_radius(a->id,mod_params->sidechain_properties), mod_params->sidechain_hbond_decay_width )) > 0.0 ) {
//	         fprintf(stderr,"found hbond %c %d G__ >> %c %d G__ %g\n",a->id,a->num,b->id,b->num,distance(a->g,b->g));
		       erg += -mod_params->sidechain_hbond_strength_s2s * intensity;
		    }
//		 }
	      }
	   }
	}
	return erg;
}


/* electrostatic energy contribution between side chain gamma atoms (G__)
   E = \frac{q_1 q_2}{\epsilon d_{12}} \exp(-d_{12}/\lambda)
   where q_1 and q_2 are the charges on atoms,
         \epsilon is the dielectric permittivity,
         \lambda is the Debye length of the electrostatic screening and
         d_{12} is the distance between the charged atoms
   with \lambda = \infty
   without a cutoff distance */
double electrostatic(Biasmap*biasmap, AA *a, AA *b, model_params *mod_params) {

	double q1, q2;
	double d2, dist;

	/* only for gamma atoms */
	if (mod_params->use_gamma_atoms != NO_GAMMA) {
	
	   /* return 0 if 0 or negative */
	   if (mod_params->recip_dielectric_param <= 0.0) return 0.0;

	   /* ignore if they are too close in the sequence */
	   /* works for multi-chain proteins */
	   if ( ( abs(a->num - b->num) < mod_params->electrostatic_min_separation ) && ( a->chainid == b->chainid ) ) return 0.0; 

	   /* ignore if either is not charged */
	   if ( a->etc & ELECTROSTATIC && b->etc & ELECTROSTATIC ) {
	      q1 = charge(a->id, mod_params->sidechain_properties);
	      q2 = charge(b->id, mod_params->sidechain_properties);
	   } else {
	      return 0.0;
	   }

	   
	   /*Nik's unphysical force */
	   /*if(q1 * q2 != 0){
		 d2 = distance(a->g,b->g);
	     dist = sqrt(d2);   
	     if(dist < 4.5) return -mod_params->recip_dielectric_param;
	     //if(dist > 5.5) return 0;
	     return -mod_params->recip_dielectric_param * q1 * q2 * (5.5 - dist); 
	   }*/
	   
	   
	   /* only count if both are charged */
	   if (q1 * q2 != 0) {
	      d2 = distance(a->g,b->g);
	      dist = sqrt(d2);
	      if(dist < 1.0) dist = 1.0;
	      if (mod_params->debye_length_param > 0.0) /* screening */ {
	         return (q1*q2)/dist * mod_params->recip_dielectric_param * exp(-dist/mod_params->debye_length_param) ;
	      } else /* no screening */ {
	         return (q1*q2)/dist * mod_params->recip_dielectric_param;
	      }
	   }
	}

	return 0.0;
}

/*
	//no charges
	if(q1 == q2 && q1 == 0) return 0.0;
	if(a->id == 'G' || a->id == 'A' || b->id == 'G' || b->id == 'A') return 0.0;
	d2 = distance(a->g,b->g);
	dist = sqrt(d2);

	//both charged	
	if(q1 * q2 != 0) return(q1*q2)/(dist*dielectric_param);
	
	//1 charged, 1 not	
	double d2_2 = 0, dist2_2 = 0, etemp = 0;
	if(a->id == 'V' || a->id  == 'T' || a->id == 'I'){
	  d2_2 = distance(a->g2,b->g);
	  dist2_2 = sqrt(d2_2);
	  etemp = 0; 
	}
	if(b->id == 'V' || b->id  == 'T' || b->id == 'I'){
	  d2_2 = distance(a->g,b->g2);
	  dist2_2 = sqrt(d2_2);
	  etemp = 0;  
	}
	
	return etemp + 0; 
*/
	
	/*double d2,q1,q2,dist;
	
	// return 0 if 0 or negative 
	if (dielectric_param <= 0.0) return 0.0;

	q1 = charge(a->id);
	q2 = charge(b->id);
	if(q1 * q2 == 0) return 0.0;
	
	d2 = distance(a->g,b->g);
	dist = sqrt(d2);

	return (q1 * q2)/(dielectric_param*dist);// *exp(-dist/debye_param);
	*/
/*}*/

double lowlevel_sbond(AA *a, AA *b, model_params *mod_params){
  double dis = sqrt(distance(a->g,b->g));
  double specific_strength;
  specific_strength = linear_decay(dis,mod_params->Sbond_distance,mod_params->Sbond_cutoff);
  if(specific_strength == 0.0) return 0.0;
  vector x, y, z;
  subtract(x, a->g, a->cb);
  subtract(y, b->g, a->g);
  subtract(z, b->cb, b->g);
  
  //double ang1 = cosine(x,y);
  //if(dis < 3.0 && (ang1 > 0.5 || ang1 < 0)) return 0.0001;
  //ang1 = cosine(y,z);
  //if(dis < 3.0 && (ang1 > 0.5 || ang1 < 0)) return 0.0001;
  
  double chi3 = -cosdihedral(x, y, z);
  if (fabs(chi3) > mod_params->Sbond_dihedral_cutoff && dis < 3.2) {
  	specific_strength = specific_strength - (fabs(chi3) - mod_params->Sbond_dihedral_cutoff);
   	//return 0.0001;
  }
  if (dis < 3.2) specific_strength *= 1.2;
	//specific_strength = 0.2 + specific_strength - 2 * (fabs(chi3) - mod_params->Sbond_dihedral_cutoff);
  // we have an S-S bond
  // compensate for CB(a)-CG(b) interactions
  //double energy_comp = vdw(a->g, b->cb, mod_params->rs + mod_params->rcb) +
  //		       vdw(a->cb, b->g, mod_params->rs + mod_params->rcb);
  //return - specific_strength * mod_params->Sbond_strength - energy_comp;
  return - specific_strength * mod_params->Sbond_strength;

}

///Gary modified version, can only handle up to 2 disulfide...might be only suitable for ADCP
double sbond_energy(int start, int end, Chain *chain,  Chaint *chaint, Biasmap *biasmap,model_params *mod_params){
  if(mod_params->Sbond_strength == 0) return 0;
  int *cyslist = NULL;
  int number_of_cys=0;
  double ans = 0;
  int i, j;
  AA *a, *b;
  if(cyslist == NULL){
    for(i = 1; i < chain->NAA; i++){
	  if(chain->aa[i].id == 'C'){
		number_of_cys++;
		cyslist = (int *) realloc(cyslist, number_of_cys *sizeof(int));  
	    cyslist[number_of_cys-1] = i;
	  }
    }
    if(number_of_cys < 2) {
		free(cyslist);
		mod_params->Sbond_strength = 0.0; 
		return 0.0;  	  
    }
  }
  

  vector *cyspos = (vector*)malloc(number_of_cys*sizeof(vector));
  for(i = 0; i < number_of_cys; i++){
	if(cyslist[i] <= end && cyslist[i] >= start){
	  a = chaint->aat + cyslist[i];  
	}
	else{
	  a = chain->aa + cyslist[i];  
	}
	cyspos[i][0] = a->g[0];
	cyspos[i][1] = a->g[1];
	cyspos[i][2] = a->g[2];
  }	  
  

  double shortestdist = 10000000;
  int shorti = 0;
  int shortj = 0;
  double temp = 10000000;
  double *cysdist = (double*)malloc(number_of_cys*number_of_cys*sizeof(double));
  for(i = 0; i < number_of_cys; i++){
	  for(j = i+1; j < number_of_cys; j++){
	  	if (cyslist[j]-cyslist[i]!=1) {
			temp = distance(cyspos[i],cyspos[j]); 
		}
		cysdist[i*number_of_cys+j] = cysdist[j*number_of_cys+i] = temp;
		if (temp < shortestdist) {
			shortestdist=temp;
			shorti = i;
			shortj = j;
		}
	  }
  }

  if (shorti+shortj == 0) {
	free(cyspos);
	free(cysdist);
    free(cyslist);
    return 0.0;
  }

  
  
  if(cyslist[shorti] <= end && cyslist[shorti] >= start){
	a = chaint->aat + cyslist[shorti];  
  }
  else{
    a = chain->aa + cyslist[shorti];  
  }

  if(cyslist[shortj] <= end && cyslist[shortj] >= start){
	b = chaint->aat + cyslist[shortj];  
  }
  else{
    b = chain->aa + cyslist[shortj];  
  }
  
  temp = lowlevel_sbond(a,b,mod_params);
  ans += temp;

  if (number_of_cys<4) {
  	free(cyspos);
	free(cysdist);
    free(cyslist);
    return ans;
  }
  temp = 0;
  shortestdist = 100000;
  int shortii=0;
  int shortjj=0;
  for(int i = 0; i < number_of_cys - 1; i++){
    for(int j = i+1; j < number_of_cys; j++) {
      if (i==shorti || j==shortj || i==shortj || j==shorti) {
	cysdist[i*number_of_cys+j] = 10000000;
	cysdist[i*number_of_cys+j] = 10000000;
      } else if (cysdist[i*number_of_cys+j]<shortestdist && cyslist[j]-cyslist[i]!=1) {
	shortestdist=cysdist[i*number_of_cys+j];
	shortii = i;
	shortjj = j;
      }
    }
  }
    
  if (shortii+shortjj == 0) {
	free(cyspos);
	free(cysdist);
    free(cyslist);
    return ans;
  }
    
  if(cyslist[shortii] <= end && cyslist[shortii] >= start){
	a = chaint->aat + cyslist[shortii];  
  }
  else{
    a = chain->aa + cyslist[shortii];  
  }

  if(cyslist[shortjj] <= end && cyslist[shortjj] >= start){
	b = chaint->aat + cyslist[shortjj];
  }
  else{
    b = chain->aa + cyslist[shortjj];  
  }
  
  temp = lowlevel_sbond(a,b,mod_params);
  ans += temp;

  
  free(cyspos);
  free(cysdist);
  free(cyslist);

  return ans;
  
  
  
  /*
  double ans = 0;
  int i, j;
  struct AA *a, *b;
  for(i = 1; i < NAA; i++){
	if(aa[i].id != 'C') continue;
	if( i <= end && i >= start){
		 a = aat+i;
    }
    else a = aa+i;
    
	int number_of_Sbond_this_cys = 0;
	for(j = 1; j < NAA; j++){
	  if(aa[j].id != 'C' || (j == i )) continue; 	
	    if( j <= end && j >= start){
		 b = aat+j;
        }
        else b = aa+j;
	    double temp = lowlevel_sbond(a,b,mod_params);
	    if(temp != 0.0){
			number_of_Sbond_this_cys++;
	        ans += temp;
	    }
	}
	if(number_of_Sbond_this_cys > 1) ans+= 1000.0;
  } 
  
  
  return ans/2;*/

	
}


/* Calculate the radius of gyration of the secondary structure elements,
   representing each secondary structure element by the centre of mass of its alpha carbons */
double secondary_radius_of_gyration(int start, int end, Chain *chain, Chaint *chaint, Biasmap *biasmap, model_params *mod_params, int which_atom, int only_hydrophobic){


  //if(mod_params->srgy_param == 0.0) return 0;

  if (which_atom != CA_ && which_atom != CB_) {
	fprintf(stderr,"WARNING! Unknown request for which_atom %x, must be one of %x or %x.  Usin CA_ as default.\n",which_atom,CA_,CB_);
	which_atom = CA_;
  }
  
  double w = 0;
  vector ri, rc = { 0.0, 0.0, 0.0 };
  vector d = { 0.0, 0.0, 0.0 };
  
  
  vector* s_com = malloc(sizeof(vector)*chain->NAA);
  int* weights = malloc(sizeof(int)*chain->NAA); int totalweight = 0;
  int state = 0;
  int count = -1;
  int state_count = 0;
  int i;
  
  //extract the C-O-Ms of the helices and strands 
  
  for(i = 1; i < chain->NAA; i++){

	if(  Distb(i, i) == state && state != 0 ){ /* still in the same secondary structure element */

	   if ((only_hydrophobic && chain->aa[i].etc & HYDROPHOBIC) || !only_hydrophobic) {
	      /* Add this amino acid */
	      state_count++;
	      if( i <= end && i >= start){ /* changed peptide section */
	         if (which_atom == CA_) {
		   add(s_com[count],s_com[count],chaint->aat[i].ca);
	         } else if (which_atom == CB_) {
		   add(s_com[count],s_com[count],chaint->aat[i].cb);
	         }
	      } else{ /* original peptide section */
	         if (which_atom == CA_) {
		   add(s_com[count],s_com[count],chain->aa[i].ca);
	         } else if (which_atom == CB_) {
		   add(s_com[count],s_com[count],chain->aa[i].cb);
	         }
	      //fprintf(stderr,"ca %f %f %f\n",aa[i].ca[0],aa[i].ca[1],aa[i].ca[2]);
              }
	   }

	} else if(Distb(i , i) != state){ /* change in the secondary structure */

	   if(state != 0){ /* end previous (nonzero) secondary structure */
		scale(s_com[count], 1.0/state_count, s_com[count]);		
		//fprintf(stderr,"c-o-m %f %f %f\n",s_com[count][0],s_com[count][1],s_com[count][2]);
		//weights[count] = state_count; totalweight += state_count;  
	        weights[count] = 1; totalweight += 1;  
	   }

	  if(Distb(i,i) != 0){ /* beginning of new secondary structure */
		state_count = 1;  
		count++;
		if( i <= end && i >= start)
		   castvec(s_com[count],chaint->aat[i].ca);	
		else
		   castvec(s_com[count],chain->aa[i].ca);
	  }

	  /* reset state */
	  state = Distb(i,i);

	}
  }
  /* end last nonzero state */
  if(state != 0){
	scale(s_com[count], 1.0/state_count, s_com[count]);
	//weights[count] = state_count; totalweight += state_count;  
	weights[count] = 1; totalweight += 1;
  }
  
  if(count < 1) {
	free(s_com);
    free(weights);
    return 0;
  }
  
  
  // calculate C-O-M of the C-O-Ms  
  for(i = 0; i <= count; i++){
	//fprintf(stderr,"Vecs: %f %f %f weight %d\n",s_com[i][0],s_com[i][1],s_com[i][2],weights[i]);
	fling(rc, rc, weights[i],s_com[i]);   
  }
  scale(rc, 1.0/totalweight, rc); 
  //fprintf(stderr,"Total weight %d Count %d\n",totalweight,count);
  //fprintf(stderr,"C-O-M: %f %f %f\n",rc[0],rc[1],rc[2]);
  
  
  //calculate the radius of gyration of the C-O-Ms 
  for (i = 0; i <= count; i++) {
	 subtract(ri, s_com[i], rc);
	 d[0] += weights[i] * ri[0] * ri[0];
	 d[1] += weights[i] * ri[1] * ri[1];
	 d[2] += weights[i] * ri[2] * ri[2];
  }
  
  scale(d, 1.0/totalweight, d);
  w = sqrt(d[0] + d[1] + d[2]);
  
  //fprintf(stderr,"D: %f %f %f\n,RGY: %f\n",d[0],d[1],d[2],w); 
   
  free(s_com);
  free(weights);

  return w;

}

int getindex(int x, int y, int z) {
	return (z * NX * NY + y * NX + x);
}


void vectorProduct(float *a, float *b, float *c) {
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

void normalizedVector(float *a, float *b, float *v) {
	int i;
	float n;

	for (i = 0; i<3; i++) {
		v[i] = b[i] - a[i];
	}
	n = 1. / sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	for (i = 0; i<3; i++) v[i] = v[i] * n;
}


//score side chain and also set gamma position
float scoreSideChain(int nbRot, int nbAtoms, double *charges, int *atypes,  double coords[nbRot][nbAtoms][3], AA *a,  int numRand)
{
	int i, j;
	float n; /* used to normalized vectors */
	float sideChainCenter[3], bestSideChainCenter[3];
	float N[3] = { a->n[0], a->n[1], a->n[2] }; /* coordiantes from 1crn.pdb:TYR29:N */
	float CA[3] = { a->ca[0], a->ca[1], a->ca[2] }; /* coordiantes from 1crn.pdb:TYR29:CA */
	float CB[3] = { a->cb[0], a->cb[1], a->cb[2] }; /* coordiantes from 1crn.pdb:TYR29:CB */
	float v1[3], v2[3], v3[3], mat[3][4]; /* used to compute xform matrix to align canonical rotamer to amino acid */
	float tc[nbRot][nbAtoms][3]; /* list of transformed coordinates */
					   /*
					   printf("VAL, %d atoms %d rotamers\n", VAL.nbAtoms, VAL.nbRot);
					   for (i=0; i<VAL.nbRot; i++) {
					   printf("Rotamer %d\n", i);
					   for (j=0; j<VAL.nbAtoms; j++) {
					   printf("  %c %f %f %f\n",VAL.atypes[j], VAL.coords[i][j][0],
					   VAL.coords[i][j][1], VAL.coords[i][j][2]);
					   }
					   }
					   */
					   /* compute matrix to align canonical TYR rotamer to 1crn:TYR29 backbone */
	float score = 0.0;
	float bestScore = 99999.0;
	float randx, randy, randz;

	int nbHeavyAtoms = 0;

	for (int nn =0; nn < nbAtoms; nn++){
		if (atypes[nn]!=3) nbHeavyAtoms++;
	}
	/*scan a little bit more space, number of random trials*/
	for (int pertInd=0; pertInd < numRand; pertInd++){
		if (pertInd!=0){
			randx = 0.5 * rand() / (double) RAND_MAX;
			randy = 0.5 * rand() / (double) RAND_MAX;
			randz = 0.5 * rand() / (double) RAND_MAX;
			N[0] = a->n[0] + randx - 0.25;
			N[1] = a->n[1] + randy - 0.25;
			N[2] = a->n[2] + randz - 0.25;
		}

		//N[3] = { a->n[0] + randx - 0.5, a->n[1] + randy - 0.5, a->n[2] + randz - 0.5 }
		normalizedVector(N, CA, v1); /* X vector */
								 /* printf("V1 %f %f %f %f\n", v1[0], v1[1], v1[2], v1[0]*v1[0]+ v1[1]*v1[1]+ v1[2]*v1[2]); */
		normalizedVector(CA, CB, v3);
		/* printf("V3 %f %f %f %f\n", v3[0], v3[1], v3[2], v3[0]*v3[0]+ v3[1]*v3[1]+ v3[2]*v3[2]); */
		vectorProduct(v3, v1, v2); /* Y vector*/
								   /* printf("V2 %f %f %f %f\n", v2[0], v2[1], v2[2], v2[0]*v2[0]+ v2[1]*v2[1]+ v2[2]*v2[2]); */
		n = 1. / sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
		for (i = 0; i < 3; i++) v2[i] = v2[i] * n;
		vectorProduct(v1, v2, v3); /* Z vector*/
		n = 1. / sqrt(v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2]);
		for (i = 0; i < 3; i++) v3[i] = v3[i] * n;
		/* printf("V3 %f %f %f %f\n", v3[0], v3[1], v3[2], v3[0]*v3[0]+ v3[1]*v3[1]+ v3[2]*v3[2]); */
		/* xform matrix */
		for (i = 0; i < 3; i++) {
			mat[i][0] = v1[i];
			mat[i][1] = v2[i];
			mat[i][2] = v3[i];
			mat[i][3] = CA[i];
		}
		/*
		printf("%8.3f %8.3f %8.3f %8.3f\n",mat[0][0],mat[0][1],mat[0][2],mat[0][3]);
		printf("%8.3f %8.3f %8.3f %8.3f\n",mat[1][0],mat[1][1],mat[1][2],mat[1][3]);
		printf("%8.3f %8.3f %8.3f %8.3f\n",mat[2][0],mat[2][1],mat[2][2],mat[2][3]);
		*/
		/* apply transformation to canonical all rot side chains coordinates */
		score = 0.0;

		for (i = 0; i < nbRot; i++) {
			score = 0.0;
			sideChainCenter[0] = 0.0;
			sideChainCenter[1] = 0.0;
			sideChainCenter[2] = 0.0;
			for (j = 0; j < nbAtoms; j++) {
				//fprintf(stderr, "test type %i score %g \n", atypes[i], score);
				tc[i][j][0] = mat[0][0] * coords[i][j][0] + mat[0][1] * coords[i][j][1] + mat[0][2] * coords[i][j][2] + mat[0][3];
				tc[i][j][1] = mat[1][0] * coords[i][j][0] + mat[1][1] * coords[i][j][1] + mat[1][2] * coords[i][j][2] + mat[1][3];
				tc[i][j][2] = mat[2][0] * coords[i][j][0] + mat[2][1] * coords[i][j][1] + mat[2][2] * coords[i][j][2] + mat[2][3];
				//fprintf(stderr, "test type %i\n", atypes[i]);
				if (atypes[j]!=3) {
					sideChainCenter[0] += tc[i][j][0];
					sideChainCenter[1] += tc[i][j][1];
					sideChainCenter[2] += tc[i][j][2];
				}

				score += gridenergy(tc[i][j][0], tc[i][j][1], tc[i][j][2], atypes[j], charges[j]);
				//fprintf(stderr, "test nbROT %i type %i score %g \n", i, atypes[j], score);
			}
			if (score < bestScore) {
				bestScore = score;
				a->SCRot = i;
				bestSideChainCenter[0] = sideChainCenter[0]/nbHeavyAtoms;
				bestSideChainCenter[1] = sideChainCenter[1]/nbHeavyAtoms;
				bestSideChainCenter[2] = sideChainCenter[2]/nbHeavyAtoms;

			}
		}
	}
	
	switch (a->id)
	{
		case 'I':
			a->g[0] = (tc[a->SCRot][0][0]+tc[a->SCRot][1][0])/2;
			a->g[1] = (tc[a->SCRot][0][1]+tc[a->SCRot][1][1])/2;
			a->g[2] = (tc[a->SCRot][0][2]+tc[a->SCRot][1][2])/2;
			a->g2[0] = tc[a->SCRot][2][0];
			a->g2[1] = tc[a->SCRot][2][1];
			a->g2[2] = tc[a->SCRot][2][2];
			break;
		case 'S':
			a->g[0] = tc[a->SCRot][0][0];
			a->g[1] = tc[a->SCRot][0][1];
			a->g[2] = tc[a->SCRot][0][2];
			break;
		case 'P':
			a->g[0] = tc[a->SCRot][0][0];
			a->g[1] = tc[a->SCRot][0][1];
			a->g[2] = tc[a->SCRot][0][2];
			break;
		default:
			a->g[0] = bestSideChainCenter[0];
			a->g[1] = bestSideChainCenter[1];
			a->g[2] = bestSideChainCenter[2];
			break;

	}
		

	//fprintf(stderr, "score %g \n", bestScore);
	//free(tc),free(v1),free(v2),free(v3),free(mat);

	return bestScore;

}

int checkClash(double x, double y, double z, double *setCoords, int ind){
	double dist = 100;

	//if (sizeof(setCoords)/sizeof(double)<ind-1) {
	//	fprintf(stderr,"error %g %d \n",sizeof(setCoords)/sizeof(double),ind); 
	//}
	//fprintf(stderr,"error %g %d \n",sizeof(setCoords)/sizeof(double),ind); 
	for (int i = 0; i < ind; i = i + 3){
		//if (setCoords)
		//fprintf(stderr,"error %g %g \n",setCoords[i], setCoords[i+1]); 
		dist = (x - setCoords[i]) * (x - setCoords[i]) + (y - setCoords[i+1]) * (y - setCoords[i+1]) + (z - setCoords[i+2]) * (z - setCoords[i+2]);
		if (dist < 4.84) {
			return 1;
		}
		
	}
	return 0;
}


double scoreSideChainNoClash(int nbRot, int nbAtoms, double charges[nbAtoms], int atypes[nbAtoms],  double coords[nbRot][nbAtoms][3], AA *a, double* setCoords, int ind, int numRand)
{
	int i, j;
	double n; /* used to normalized vectors */
	double sideChainCenter[3], bestSideChainCenter[3];
	float N[3] = { a->n[0], a->n[1], a->n[2] }; /* coordinates from 1crn.pdb:TYR29:N */
	float CA[3] = { a->ca[0], a->ca[1], a->ca[2] }; /* coordinates from 1crn.pdb:TYR29:CA */
	float CB[3] = { a->cb[0], a->cb[1], a->cb[2] }; /* coordinates from 1crn.pdb:TYR29:CB */
	float v1[3], v2[3], v3[3], mat[3][4]; /* used to compute xform matrix to align canonical rotamer to amino acid */
	float tc[nbRot][nbAtoms][3]; /* list of transformed coordinates */
					   /*
					   printf("VAL, %d atoms %d rotamers\n", VAL.nbAtoms, VAL.nbRot);
					   for (i=0; i<VAL.nbRot; i++) {
					   printf("Rotamer %d\n", i);
					   for (j=0; j<VAL.nbAtoms; j++) {
					   printf("  %c %f %f %f\n",VAL.atypes[j], VAL.coords[i][j][0],
					   VAL.coords[i][j][1], VAL.coords[i][j][2]);
					   }
					   }
					   */
					   /* compute matrix to align canonical TYR rotamer to 1crn:TYR29 backbone */
	double score = 0.0;
	double bestScore = 99999.0;
	double randx, randy, randz;
    int clash = 0;
	int nbHeavyAtoms = 0;

	for (int nn =0; nn < nbAtoms; nn++){
		if (atypes[nn]!=3) nbHeavyAtoms++;
	}
	/*scan a little bit more space, number of random trials*/
	for (int pertInd=0; pertInd < numRand; pertInd++){
		if (pertInd!=0){
			randx = 0.5 * rand() / (double) RAND_MAX;
			randy = 0.5 * rand() / (double) RAND_MAX;
			randz = 0.5 * rand() / (double) RAND_MAX;
			N[0] = a->n[0] + randx - 0.25;
			N[1] = a->n[1] + randy - 0.25;
			N[2] = a->n[2] + randz - 0.25;
		}

		//N[3] = { a->n[0] + randx - 0.5, a->n[1] + randy - 0.5, a->n[2] + randz - 0.5 }
		normalizedVector(N, CA, v1); /* X vector */
								 /* printf("V1 %f %f %f %f\n", v1[0], v1[1], v1[2], v1[0]*v1[0]+ v1[1]*v1[1]+ v1[2]*v1[2]); */
		normalizedVector(CA, CB, v3);
		/* printf("V3 %f %f %f %f\n", v3[0], v3[1], v3[2], v3[0]*v3[0]+ v3[1]*v3[1]+ v3[2]*v3[2]); */
		vectorProduct(v3, v1, v2); /* Y vector*/
								   /* printf("V2 %f %f %f %f\n", v2[0], v2[1], v2[2], v2[0]*v2[0]+ v2[1]*v2[1]+ v2[2]*v2[2]); */
		n = 1. / sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
		for (i = 0; i < 3; i++) v2[i] = v2[i] * n;
		vectorProduct(v1, v2, v3); /* Z vector*/
		n = 1. / sqrt(v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2]);
		for (i = 0; i < 3; i++) v3[i] = v3[i] * n;
		/* printf("V3 %f %f %f %f\n", v3[0], v3[1], v3[2], v3[0]*v3[0]+ v3[1]*v3[1]+ v3[2]*v3[2]); */
		/* xform matrix */
		for (i = 0; i < 3; i++) {
			mat[i][0] = v1[i];
			mat[i][1] = v2[i];
			mat[i][2] = v3[i];
			mat[i][3] = CA[i];
		}
		/*
		printf("%8.3f %8.3f %8.3f %8.3f\n",mat[0][0],mat[0][1],mat[0][2],mat[0][3]);
		printf("%8.3f %8.3f %8.3f %8.3f\n",mat[1][0],mat[1][1],mat[1][2],mat[1][3]);
		printf("%8.3f %8.3f %8.3f %8.3f\n",mat[2][0],mat[2][1],mat[2][2],mat[2][3]);
		*/
		/* apply transformation to canonical all rot side chains coordinates */
                int rotCount = 0;
		for (i = 0; i < nbRot; i++) {
			if (nbRot > 20) {
				i = rand()%nbRot;
				rotCount ++;
				if (rotCount > 20) break;
			}
			score = 0.0;
			clash = 0;
			sideChainCenter[0] = 0.0;
			sideChainCenter[1] = 0.0;
			sideChainCenter[2] = 0.0;
			for (j = 0; j < nbAtoms; j++) {
				//fprintf(stderr, "test type %i score %g \n", atypes[i], score);
				tc[i][j][0] = mat[0][0] * coords[i][j][0] + mat[0][1] * coords[i][j][1] + mat[0][2] * coords[i][j][2] + mat[0][3];
				tc[i][j][1] = mat[1][0] * coords[i][j][0] + mat[1][1] * coords[i][j][1] + mat[1][2] * coords[i][j][2] + mat[1][3];
				tc[i][j][2] = mat[2][0] * coords[i][j][0] + mat[2][1] * coords[i][j][1] + mat[2][2] * coords[i][j][2] + mat[2][3];
				//fprintf(stderr, "test type %i\n", atypes[i]);
				if (atypes[j]!=3) { // ignore H atoms when calculating center of side chain
					sideChainCenter[0] += tc[i][j][0];
					sideChainCenter[1] += tc[i][j][1];
					sideChainCenter[2] += tc[i][j][2];
					clash = checkClash(tc[i][j][0], tc[i][j][1], tc[i][j][2], setCoords, ind);
					if (clash) score += 6.5;
				}

				score += gridenergy(tc[i][j][0], tc[i][j][1], tc[i][j][2], atypes[j], charges[j]);
			}
			//fprintf(stderr, "num %d id %c test nbROT %i type %i score %g \n",a->num,a->id, i, atypes[j], score);
			//if (clash) break;
			if (score < bestScore) {
				bestScore = score;
				a->SCRot = i;
				bestSideChainCenter[0] = sideChainCenter[0]/nbHeavyAtoms;
				bestSideChainCenter[1] = sideChainCenter[1]/nbHeavyAtoms;
				bestSideChainCenter[2] = sideChainCenter[2]/nbHeavyAtoms;
			}
		}
	}

	if (bestScore>90000) return 1000.0;

	switch (a->id)
	{
		case 'I':
			a->g[0] = (tc[a->SCRot][0][0]+tc[a->SCRot][1][0])/2;
			a->g[1] = (tc[a->SCRot][0][1]+tc[a->SCRot][1][1])/2;
			a->g[2] = (tc[a->SCRot][0][2]+tc[a->SCRot][1][2])/2;
			a->g2[0] = tc[a->SCRot][2][0];
			a->g2[1] = tc[a->SCRot][2][1];
			a->g2[2] = tc[a->SCRot][2][2];
			break;
		case 'S':
			a->g[0] = tc[a->SCRot][0][0];
			a->g[1] = tc[a->SCRot][0][1];
			a->g[2] = tc[a->SCRot][0][2];
			break;
		case 'P':
			a->g[0] = tc[a->SCRot][0][0];
			a->g[1] = tc[a->SCRot][0][1];
			a->g[2] = tc[a->SCRot][0][2];
			break;
		case 'T':
			a->g[0] = tc[a->SCRot][1][0];
			a->g[1] = tc[a->SCRot][1][1];
			a->g[2] = tc[a->SCRot][1][2];
			a->g2[0] = tc[a->SCRot][0][0];
			a->g2[1] = tc[a->SCRot][0][1];
			a->g2[2] = tc[a->SCRot][0][2];
			break;
		case 'V':
			a->g[0] = tc[a->SCRot][1][0];
			a->g[1] = tc[a->SCRot][1][1];
			a->g[2] = tc[a->SCRot][1][2];
			a->g2[0] = tc[a->SCRot][0][0];
			a->g2[1] = tc[a->SCRot][0][1];
			a->g2[2] = tc[a->SCRot][0][2];
			break;
		default:
                        //a->g[0] = tc[a->SCRot][0][0];
                        //a->g[1] = tc[a->SCRot][0][1];
                        //a->g[2] = tc[a->SCRot][0][2];

			a->g[0] = bestSideChainCenter[0];
			a->g[1] = bestSideChainCenter[1];
			a->g[2] = bestSideChainCenter[2];
			break;

	}
		

	//fprintf(stderr, "score %g \n", bestScore);
	//free(tc),free(v1),free(v2),free(v3),free(mat);

	return bestScore;

}


//score side chain and also set gamma position
float scoreSideChainNew(AA *a,  int numRand)
{
	int nbRot = _AASCRotTable[a->sideChainTemplateIndex].nbRot;
	int nbAtoms = _AASCRotTable[a->sideChainTemplateIndex].nbAtoms;
	double *charges = _AASCRotTable[a->sideChainTemplateIndex].charges;
	int *atypes = _AASCRotTable[a->sideChainTemplateIndex].atypes;
	double **coords = _AASCRotTable[a->sideChainTemplateIndex].coords;

	int i, j;
	float n; /* used to normalized vectors */
	float sideChainCenter[3], bestSideChainCenter[3];
	float N[3] = { a->n[0], a->n[1], a->n[2] }; /* coordiantes from 1crn.pdb:TYR29:N */
	float CA[3] = { a->ca[0], a->ca[1], a->ca[2] }; /* coordiantes from 1crn.pdb:TYR29:CA */
	float CB[3] = { a->cb[0], a->cb[1], a->cb[2] }; /* coordiantes from 1crn.pdb:TYR29:CB */
	float v1[3], v2[3], v3[3], mat[3][4]; /* used to compute xform matrix to align canonical rotamer to amino acid */
	float tc[nbRot][nbAtoms][3]; /* list of transformed coordinates */
					   /*
					   printf("VAL, %d atoms %d rotamers\n", VAL.nbAtoms, VAL.nbRot);
					   for (i=0; i<VAL.nbRot; i++) {
					   printf("Rotamer %d\n", i);
					   for (j=0; j<VAL.nbAtoms; j++) {
					   printf("  %c %f %f %f\n",VAL.atypes[j], VAL.coords[i][j][0],
					   VAL.coords[i][j][1], VAL.coords[i][j][2]);
					   }
					   }
					   */
					   /* compute matrix to align canonical TYR rotamer to 1crn:TYR29 backbone */
	float score = 0.0;
	float bestScore = 99999.0;
	float randx, randy, randz;

	int nbHeavyAtoms = 0;

	for (int nn =0; nn < nbAtoms; nn++){
		if (atypes[nn]!=3) nbHeavyAtoms++;
	}
	/*scan a little bit more space, number of random trials*/
	for (int pertInd=0; pertInd < numRand; pertInd++){
		if (pertInd!=0){
			randx = 0.5 * rand() / (double) RAND_MAX;
			randy = 0.5 * rand() / (double) RAND_MAX;
			randz = 0.5 * rand() / (double) RAND_MAX;
			N[0] = a->n[0] + randx - 0.25;
			N[1] = a->n[1] + randy - 0.25;
			N[2] = a->n[2] + randz - 0.25;
		}

		//N[3] = { a->n[0] + randx - 0.5, a->n[1] + randy - 0.5, a->n[2] + randz - 0.5 }
		normalizedVector(N, CA, v1); /* X vector */
								 /* printf("V1 %f %f %f %f\n", v1[0], v1[1], v1[2], v1[0]*v1[0]+ v1[1]*v1[1]+ v1[2]*v1[2]); */
		normalizedVector(CA, CB, v3);
		/* printf("V3 %f %f %f %f\n", v3[0], v3[1], v3[2], v3[0]*v3[0]+ v3[1]*v3[1]+ v3[2]*v3[2]); */
		vectorProduct(v3, v1, v2); /* Y vector*/
								   /* printf("V2 %f %f %f %f\n", v2[0], v2[1], v2[2], v2[0]*v2[0]+ v2[1]*v2[1]+ v2[2]*v2[2]); */
		n = 1. / sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
		for (i = 0; i < 3; i++) v2[i] = v2[i] * n;
		vectorProduct(v1, v2, v3); /* Z vector*/
		n = 1. / sqrt(v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2]);
		for (i = 0; i < 3; i++) v3[i] = v3[i] * n;
		/* printf("V3 %f %f %f %f\n", v3[0], v3[1], v3[2], v3[0]*v3[0]+ v3[1]*v3[1]+ v3[2]*v3[2]); */
		/* xform matrix */
		for (i = 0; i < 3; i++) {
			mat[i][0] = v1[i];
			mat[i][1] = v2[i];
			mat[i][2] = v3[i];
			mat[i][3] = CA[i];
		}
		/*
		printf("%8.3f %8.3f %8.3f %8.3f\n",mat[0][0],mat[0][1],mat[0][2],mat[0][3]);
		printf("%8.3f %8.3f %8.3f %8.3f\n",mat[1][0],mat[1][1],mat[1][2],mat[1][3]);
		printf("%8.3f %8.3f %8.3f %8.3f\n",mat[2][0],mat[2][1],mat[2][2],mat[2][3]);
		*/
		/* apply transformation to canonical all rot side chains coordinates */
		score = 0.0;

		for (i = 0; i < nbRot; i++) {
			score = 0.0;
			sideChainCenter[0] = 0.0;
			sideChainCenter[1] = 0.0;
			sideChainCenter[2] = 0.0;
			for (j = 0; j < nbAtoms; j++) {
				//fprintf(stderr, "test type %i score %g \n", atypes[i], score);
			        double x = coords[i][3*j];
				double y = coords[i][3*j+1];
				double z = coords[i][3*j+2];
				tc[i][j][0] = mat[0][0] * x + mat[0][1] * y + mat[0][2] * z + mat[0][3];
				tc[i][j][1] = mat[1][0] * x + mat[1][1] * y + mat[1][2] * z + mat[1][3];
				tc[i][j][2] = mat[2][0] * x + mat[2][1] * y + mat[2][2] * z + mat[2][3];
				//fprintf(stderr, "test type %i\n", atypes[i]);
				if (atypes[j]!=3) {
					sideChainCenter[0] += tc[i][j][0];
					sideChainCenter[1] += tc[i][j][1];
					sideChainCenter[2] += tc[i][j][2];
				}

				score += gridenergy(tc[i][j][0], tc[i][j][1], tc[i][j][2], atypes[j], charges[j]);
				//fprintf(stderr, "test nbROT %i type %i score %g \n", i, atypes[j], score);
			}
			if (score < bestScore) {
				bestScore = score;
				a->SCRot = i;
				bestSideChainCenter[0] = sideChainCenter[0]/nbHeavyAtoms;
				bestSideChainCenter[1] = sideChainCenter[1]/nbHeavyAtoms;
				bestSideChainCenter[2] = sideChainCenter[2]/nbHeavyAtoms;

			}
		}
	}

	//FIXME MS should not rely on AA type to set a->g and a->g2
	// although it seems this function is only called fro PRO
	switch (a->id)
	{
		case 'I':
			a->g[0] = (tc[a->SCRot][0][0]+tc[a->SCRot][1][0])/2;
			a->g[1] = (tc[a->SCRot][0][1]+tc[a->SCRot][1][1])/2;
			a->g[2] = (tc[a->SCRot][0][2]+tc[a->SCRot][1][2])/2;
			a->g2[0] = tc[a->SCRot][2][0];
			a->g2[1] = tc[a->SCRot][2][1];
			a->g2[2] = tc[a->SCRot][2][2];
			break;
		case 'S':
			a->g[0] = tc[a->SCRot][0][0];
			a->g[1] = tc[a->SCRot][0][1];
			a->g[2] = tc[a->SCRot][0][2];
			break;
		case 'P':
		  //MS this caused proline righ (D) to look bad and CD to be too close to C==O
		  //    a->g[0] = tc[a->SCRot][0][0];
		  //	a->g[1] = tc[a->SCRot][0][1];
		  //	a->g[2] = tc[a->SCRot][0][2];
			break;
		default:
			a->g[0] = bestSideChainCenter[0];
			a->g[1] = bestSideChainCenter[1];
			a->g[2] = bestSideChainCenter[2];
			break;

	}
		

	//fprintf(stderr, "score %g \n", bestScore);
	//free(tc),free(v1),free(v2),free(v3),free(mat);

	return bestScore;

}

void buildSideChain(AA *a, float *coords)
{
  // build the side chain atom coordinates for residue 'a' using rotamer rotNum
  // the coordinates are stored in 'coords' which managed by the caller
  
  int i, j;
  int nbAtoms = _AASCRotTable[a->sideChainTemplateIndex].nbAtoms;
  int rotNum = a->SCRot; // rotamer number
  float N[3] = { a->n[0], a->n[1], a->n[2] }; /* coordinates from 1crn.pdb:TYR29:N */
  float CA[3] = { a->ca[0], a->ca[1], a->ca[2] }; /* coordinates from 1crn.pdb:TYR29:CA */
  float CB[3] = { a->cb[0], a->cb[1], a->cb[2] }; /* coordinates from 1crn.pdb:TYR29:CB */
  float v1[3], v2[3], v3[3], mat[3][4]; /* used to compute xform matrix to align canonical rotamer to amino acid */
  double **can_coords = _AASCRotTable[a->sideChainTemplateIndex].coords;
  double n, x, y, z; /* used to normalized vectors */

  // build transformation matrix to transform caninocal rotamer to current residue
  normalizedVector(N, CA, v1); /* X vector */
  normalizedVector(CA, CB, v3);
  vectorProduct(v3, v1, v2); /* Y vector*/
  n = 1. / sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
  for (i = 0; i < 3; i++) v2[i] = v2[i] * n;
  vectorProduct(v1, v2, v3); /* Z vector*/
  n = 1. / sqrt(v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2]);
  for (i = 0; i < 3; i++) v3[i] = v3[i] * n;

  /* xform matrix */
  for (i = 0; i < 3; i++) {
    mat[i][0] = v1[i];
    mat[i][1] = v2[i];
    mat[i][2] = v3[i];
    mat[i][3] = CA[i];
  }
  for (j = 0; j < nbAtoms; j++) {
    x = can_coords[rotNum][3*j];
    y = can_coords[rotNum][3*j+1];
    z = can_coords[rotNum][3*j+2];
    coords[j*3]   = mat[0][0] * x + mat[0][1] * y + mat[0][2] * z + mat[0][3];
    coords[j*3+1] = mat[1][0] * x + mat[1][1] * y + mat[1][2] * z + mat[1][3];
    coords[j*3+2] = mat[2][0] * x + mat[2][1] * y + mat[2][2] * z + mat[2][3];
  }
}
	
double scoreSideChainNoClashNew(AA *a, double* setCoords, int ind, int numRand)
{
  // for histidine we consider rotamers for HIE and HID protonation states
  int nbRot = _AASCRotTable[a->sideChainTemplateIndex].nbRot;
  int nbAtoms = _AASCRotTable[a->sideChainTemplateIndex].nbAtoms;

  double *charges = _AASCRotTable[a->sideChainTemplateIndex].charges;
  int *atypes = _AASCRotTable[a->sideChainTemplateIndex].atypes;
  double **coords = _AASCRotTable[a->sideChainTemplateIndex].coords;

  double *charges2 = _AASCRotTable[a->sideChainTemplateIndex].charges;
  int *atypes2 = _AASCRotTable[a->sideChainTemplateIndex].atypes;
  double **coords2 = _AASCRotTable[a->sideChainTemplateIndex].coords;

  // MS this is brittle as is depend on the order of the rotamers in rotamers.lib
  //int isHis = a->sideChainTemplateIndex==6 || a->sideChainTemplateIndex==7 || a->sideChainTemplateIndex==8;

  // MS maybe we shoudl implement if input was 'h' we try HIE or HID, if input was h<HIE> we only try HIE
  //  but for now the best one overwrites a->sideChainTemplateIndex
  int isHis = strcmp(_AASCRotTable[a->sideChainTemplateIndex].name, "HIS")==0 || strcmp(_AASCRotTable[a->sideChainTemplateIndex].name, "HIE")==0 || strcmp(_AASCRotTable[a->sideChainTemplateIndex].name, "HID")==0;
  int indHIE = getSideChainTemplateIndexFromName("HIE");
  int indHID = getSideChainTemplateIndexFromName("HID");
  if (isHis) { // HIS or HIE or HID
    charges = _AASCRotTable[indHIE].charges; // HIE
    atypes = _AASCRotTable[indHIE].atypes;
    coords = _AASCRotTable[indHIE].coords;

    charges2 = _AASCRotTable[indHID].charges; // HID
    atypes2 = _AASCRotTable[indHID].atypes;
    coords2 = _AASCRotTable[indHID].coords;
  } else {
    charges = _AASCRotTable[a->sideChainTemplateIndex].charges;
    atypes = _AASCRotTable[a->sideChainTemplateIndex].atypes;
    coords = _AASCRotTable[a->sideChainTemplateIndex].coords;
  }
    
  int i, j;
  double n; /* used to normalized vectors */
  double sideChainCenter[3], bestSideChainCenter[3];
  float N[3] = { a->n[0], a->n[1], a->n[2] }; /* coordinates from 1crn.pdb:TYR29:N */
  float CA[3] = { a->ca[0], a->ca[1], a->ca[2] }; /* coordinates from 1crn.pdb:TYR29:CA */
  float CB[3] = { a->cb[0], a->cb[1], a->cb[2] }; /* coordinates from 1crn.pdb:TYR29:CB */
  float v1[3], v2[3], v3[3], mat[3][4]; /* used to compute xform matrix to align canonical rotamer to amino acid */
  float tc[nbRot][nbAtoms][3]; /* list of transformed coordinates */
                               /*
				 printf("VAL, %d atoms %d rotamers\n", VAL.nbAtoms, VAL.nbRot);
				 for (i=0; i<VAL.nbRot; i++) {
				 printf("Rotamer %d\n", i);
				 for (j=0; j<VAL.nbAtoms; j++) {
				 printf("  %c %f %f %f\n",VAL.atypes[j], VAL.coords[i][j][0],
				 VAL.coords[i][j][1], VAL.coords[i][j][2]);
				 }
				 }
			       */

  double score = 0.0;  // score for a given rotameric side chain
  double score2 = 0.0;  // score for a given rotameric side chain for HID
  double bestScore = 99999.0; //best score across rotamers
  double randx, randy, randz;
  int clash = 0;
  int nbHeavyAtoms = 0;

  // count heavy atoms
  //
  for (int nn =0; nn < nbAtoms; nn++){
    if (atypes[nn]!=3) nbHeavyAtoms++; // type 3 is 'HD'
  }

  /*scan a little bit more space, number of random trials*/
  /* FIXME MS if randnum is > 1 this will not work because tc get overwritten */
  for (int pertInd=0; pertInd < numRand; pertInd++){
    if (pertInd!=0){
      randx = 0.5 * rand() / (double) RAND_MAX;
      randy = 0.5 * rand() / (double) RAND_MAX;
      randz = 0.5 * rand() / (double) RAND_MAX;
      N[0] = a->n[0] + randx - 0.25;
      N[1] = a->n[1] + randy - 0.25;
      N[2] = a->n[2] + randz - 0.25;
    }
    
    // compute xform matrix to align canonical rotamer to amino acid
    //
    //N[3] = { a->n[0] + randx - 0.5, a->n[1] + randy - 0.5, a->n[2] + randz - 0.5 }
    normalizedVector(N, CA, v1); /* X vector */
    /* printf("V1 %f %f %f %f\n", v1[0], v1[1], v1[2], v1[0]*v1[0]+ v1[1]*v1[1]+ v1[2]*v1[2]); */
    normalizedVector(CA, CB, v3);
    /* printf("V3 %f %f %f %f\n", v3[0], v3[1], v3[2], v3[0]*v3[0]+ v3[1]*v3[1]+ v3[2]*v3[2]); */
    vectorProduct(v3, v1, v2); /* Y vector*/
    /* printf("V2 %f %f %f %f\n", v2[0], v2[1], v2[2], v2[0]*v2[0]+ v2[1]*v2[1]+ v2[2]*v2[2]); */
    n = 1. / sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
    for (i = 0; i < 3; i++) v2[i] = v2[i] * n;
    vectorProduct(v1, v2, v3); /* Z vector*/
    n = 1. / sqrt(v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2]);
    for (i = 0; i < 3; i++) v3[i] = v3[i] * n;
    /* printf("V3 %f %f %f %f\n", v3[0], v3[1], v3[2], v3[0]*v3[0]+ v3[1]*v3[1]+ v3[2]*v3[2]); */
    /* xform matrix to align canonical rotamer to amino acid in peptide */
    for (i = 0; i < 3; i++) {
      mat[i][0] = v1[i];
      mat[i][1] = v2[i];
      mat[i][2] = v3[i];
      mat[i][3] = CA[i];
    }
    /*
      printf("%8.3f %8.3f %8.3f %8.3f\n",mat[0][0],mat[0][1],mat[0][2],mat[0][3]);
      printf("%8.3f %8.3f %8.3f %8.3f\n",mat[1][0],mat[1][1],mat[1][2],mat[1][3]);
      printf("%8.3f %8.3f %8.3f %8.3f\n",mat[2][0],mat[2][1],mat[2][2],mat[2][3]);
    */

    // apply transformation to canonical all rot side chains coordinates and calculate grid score
    //
    int rotCount = 0;
    for (i = 0; i < nbRot; i++) {
      if (nbRot > 20) {
	i = rand()%nbRot;
	rotCount ++;
	if (rotCount > 20) break;
      }
      score = 0.0;
      clash = 0;
      sideChainCenter[0] = 0.0;
      sideChainCenter[1] = 0.0;
      sideChainCenter[2] = 0.0;
      for (j = 0; j < nbAtoms; j++) {
	//fprintf(stderr, "test type %i score %g \n", atypes[i], score);
	double x = coords[i][3*j];
	double y = coords[i][3*j+1];
	double z = coords[i][3*j+2];
	tc[i][j][0] = mat[0][0] * x + mat[0][1] * y + mat[0][2] * z + mat[0][3];
	tc[i][j][1] = mat[1][0] * x + mat[1][1] * y + mat[1][2] * z + mat[1][3];
	tc[i][j][2] = mat[2][0] * x + mat[2][1] * y + mat[2][2] * z + mat[2][3];
	//fprintf(stderr, "test type %i\n", atypes[i]);
	if (atypes[j]!=3) { // ignore H atoms when calculating center of side chain
	  sideChainCenter[0] += tc[i][j][0];
	  sideChainCenter[1] += tc[i][j][1];
	  sideChainCenter[2] += tc[i][j][2];
	  clash = checkClash(tc[i][j][0], tc[i][j][1], tc[i][j][2], setCoords, ind);
	  if (clash) score += 6.5;
	}
	score += gridenergy(tc[i][j][0], tc[i][j][1], tc[i][j][2], atypes[j], charges[j]);
      }

      if (isHis) { // HIS or HIE or HID the previous loop used HIE rotamers
	for (j = 0; j < nbAtoms; j++) {   // here we use the HID rotamers 
	  //fprintf(stderr, "test type %i score %g \n", atypes[i], score);
	  double x = coords2[i][3*j];
	  double y = coords2[i][3*j+1];
	  double z = coords2[i][3*j+2];
	  tc[i][j][0] = mat[0][0] * x + mat[0][1] * y + mat[0][2] * z + mat[0][3];
	  tc[i][j][1] = mat[1][0] * x + mat[1][1] * y + mat[1][2] * z + mat[1][3];
	  tc[i][j][2] = mat[2][0] * x + mat[2][1] * y + mat[2][2] * z + mat[2][3];
	  //fprintf(stderr, "test type %i\n", atypes[i]);
	  if (atypes2[j]!=3) { // ignore H atoms when calculating center of side chain
	    sideChainCenter[0] += tc[i][j][0];
	    sideChainCenter[1] += tc[i][j][1];
	    sideChainCenter[2] += tc[i][j][2];
	    clash = checkClash(tc[i][j][0], tc[i][j][1], tc[i][j][2], setCoords, ind);
	    if (clash) score += 6.5;
	  }
	  score2 += gridenergy(tc[i][j][0], tc[i][j][1], tc[i][j][2], atypes2[j], charges2[j]);
	}
      }
      //fprintf(stderr, "num %d id %c test nbROT %i type %i score %g \n",a->num,a->id, i, atypes[j], score);
      //if (clash) break;
      if (score < bestScore) {
	bestScore = score;
	a->SCRot = i;
	bestSideChainCenter[0] = sideChainCenter[0]/nbHeavyAtoms;
	bestSideChainCenter[1] = sideChainCenter[1]/nbHeavyAtoms;
	bestSideChainCenter[2] = sideChainCenter[2]/nbHeavyAtoms;
	if (isHis) a->sideChainTemplateIndex=indHIE;
      }
      if (isHis && score2 < bestScore) {
	bestScore = score2;
	a->SCRot = i;
	bestSideChainCenter[0] = sideChainCenter[0]/nbHeavyAtoms;
	bestSideChainCenter[1] = sideChainCenter[1]/nbHeavyAtoms;
	bestSideChainCenter[2] = sideChainCenter[2]/nbHeavyAtoms;
	if (isHis) a->sideChainTemplateIndex=indHID;
      }
    }
  }

  if (bestScore>90000) return 1000.0;

  switch (a->id)
    {
    case 'I':
      a->g[0] = (tc[a->SCRot][0][0]+tc[a->SCRot][1][0])/2;
      a->g[1] = (tc[a->SCRot][0][1]+tc[a->SCRot][1][1])/2;
      a->g[2] = (tc[a->SCRot][0][2]+tc[a->SCRot][1][2])/2;
      a->g2[0] = tc[a->SCRot][2][0];
      a->g2[1] = tc[a->SCRot][2][1];
      a->g2[2] = tc[a->SCRot][2][2];
      break;
    case 'S':
      a->g[0] = tc[a->SCRot][0][0];
      a->g[1] = tc[a->SCRot][0][1];
      a->g[2] = tc[a->SCRot][0][2];
      break;
    case 'P':
      a->g[0] = tc[a->SCRot][0][0];
      a->g[1] = tc[a->SCRot][0][1];
      a->g[2] = tc[a->SCRot][0][2];
      break;
    case 'T':
      a->g[0] = tc[a->SCRot][1][0];
      a->g[1] = tc[a->SCRot][1][1];
      a->g[2] = tc[a->SCRot][1][2];
      a->g2[0] = tc[a->SCRot][0][0];
      a->g2[1] = tc[a->SCRot][0][1];
      a->g2[2] = tc[a->SCRot][0][2];
      break;
    case 'V':
      a->g[0] = tc[a->SCRot][1][0];
      a->g[1] = tc[a->SCRot][1][1];
      a->g[2] = tc[a->SCRot][1][2];
      a->g2[0] = tc[a->SCRot][0][0];
      a->g2[1] = tc[a->SCRot][0][1];
      a->g2[2] = tc[a->SCRot][0][2];
      break;
    default:
      //a->g[0] = tc[a->SCRot][0][0];
      //a->g[1] = tc[a->SCRot][0][1];
      //a->g[2] = tc[a->SCRot][0][2];

      a->g[0] = bestSideChainCenter[0];
      a->g[1] = bestSideChainCenter[1];
      a->g[2] = bestSideChainCenter[2];
      break;
      
    }
		

  //fprintf(stderr, "score %g \n", bestScore);
  //free(tc),free(v1),free(v2),free(v3),free(mat);

  return bestScore;
}

double gridenergy(double X, double Y, double Z, int i, double charge) {
	//fprintf(stderr, "X %g Y %g Z %g charge \n", X, Y, Z, i);
	double erg = 0.0;
	
	double exactGridX = (X - centerX) / spacing + (NX - 1) / 2;
	double exactGridY = (Y - centerY) / spacing + (NY - 1) / 2;
	double exactGridZ = (Z - centerZ) / spacing + (NZ - 1) / 2;
	double perAtomtype = 0.0, deSolv = 0.0, eStatic = 0.0;
	//fprintf(stderr, "type %i charge %g \n", i, charge);
	double *mapvalue = gridmapvalues[i];
	//double *emapvalue = gridmapvalues[8];
	//double *dmapvalue = gridmapvalues[9];
	//fprintf(stderr, "type %i charge %g \n", i, charge);
	/* elements are 0:C, 1:N, 2:O, 3:H, 4:S, 5:CA, 6:NA ,7:f, 8:elec, 9:desolv     */


	double abscharge = (charge >= 0. ? charge : -charge);
    //charge = 0.;
	int lowGridX = (int)exactGridX,
		lowGridY = (int)exactGridY,
		lowGridZ = (int)exactGridZ;
	double highFracX = exactGridX - lowGridX,
		highFracY = exactGridY - lowGridY,
		highFracZ = exactGridZ - lowGridZ;
	double lowFracX = 1. - highFracX,
		lowFracY = 1. - highFracY,
		lowFracZ = 1. - highFracZ;
	double lowLowLowFrac = lowFracX * lowFracY * lowFracZ,
		lowLowHighFrac = lowFracX * lowFracY * highFracZ,
		lowHighLowFrac = lowFracX * highFracY * lowFracZ,
		lowHighHighFrac = lowFracX * highFracY * highFracZ,
		highLowLowFrac = highFracX * lowFracY * lowFracZ,
		highLowHighFrac = highFracX * lowFracY * highFracZ,
		highHighLowFrac = highFracX * highFracY * lowFracZ,
		highHighHighFrac = highFracX * highFracY * highFracZ;

	int lowLowLowIndex = getindex(exactGridX, exactGridY, exactGridZ);
	//if (lowLowLowIndex < 0 || lowLowLowIndex > NX*NY*NZ) return 0;
	int	lowLowHighIndex = lowLowLowIndex + NX * NY,
		lowHighLowIndex = lowLowLowIndex + NX,
		lowHighHighIndex = lowLowHighIndex + NX,
		highLowLowIndex = lowLowLowIndex + 1,
		highLowHighIndex = lowLowHighIndex + 1,
		highHighLowIndex = lowHighLowIndex + 1,
		highHighHighIndex = lowHighHighIndex + 1;

	int outofBox = 0;
	double outofBoxPen = 0.0;
        double outD = 0.0;
	if (exactGridX <= 0 || exactGridX >= NX - 1) {
                outD = exactGridX <= 0 ?  (-exactGridX) : (exactGridX - NX + 1);
                outofBoxPen = outD*outD;
		//outofBoxPen = ((exactGridX - (NX - 1) / 2)*(exactGridX - (NX - 1) / 2)) / 20;
		//outofBoxPen = ((exactGridX - (NX - 1) / 2)*(exactGridX - (NX - 1) / 2)) - (NX - 1) * (NX - 1) / 4;
		outofBox = 1;
		if (outofBoxPen > 1000000000) {
			fprintf(stderr, "xX %g Y %g Z %g Erg %g \n", exactGridX, exactGridY, exactGridZ, outofBoxPen);
			return erg + 10;
		}
		erg += outofBoxPen;
                //return 2.;
	}
	if (exactGridY <= 0 || exactGridY >= NY - 1) {
                outD = exactGridY <= 0 ?  (-exactGridY) : (exactGridY - NY + 1);
                outofBoxPen = outD*outD;
		//outofBoxPen = ((exactGridY - (NY - 1) / 2)*(exactGridY - (NY - 1) / 2)) / 20;
		//outofBoxPen = ((exactGridY - (NY - 1) / 2)*(exactGridY - (NY - 1) / 2)) - (NY - 1) * (NY - 1) / 4;
		outofBox = 1;
		if (outofBoxPen > 1000000000) {
			fprintf(stderr, "X %g yY %g Z %g Erg %g \n", exactGridX, exactGridY, exactGridZ, outofBoxPen);
			return erg + 10;
		}
		erg += outofBoxPen;
                //return 2.;
	}
	if (exactGridZ <= 0 || exactGridZ >= NZ - 1) {
                outD = exactGridZ <= 0 ?  (-exactGridZ) : (exactGridZ - NZ + 1);
                outofBoxPen = outD*outD;

		//outofBoxPen = ((exactGridZ - (NZ - 1) / 2)*(exactGridZ - (NZ - 1) / 2)) / 20;
		//outofBoxPen = ((exactGridZ - (NZ - 1) / 2)*(exactGridZ - (NZ - 1) / 2)) - (NZ - 1) * (NZ - 1) / 4;
		outofBox = 1;
		if (outofBoxPen > 1000000000) {
			fprintf(stderr, "X %g Y %g zZ %g Erg %g \n", exactGridX, exactGridY, exactGridZ, outofBoxPen);
			return erg + 10;
		}
		erg += outofBoxPen;
                //return 2.;
	}
	if (outofBox)
                return erg>10000?10000.:(erg/100); 
		//fprintf(stderr, "X %g Y %g Z %g Erg %g \n", exactGridX, exactGridY, exactGridZ, outofBoxPen);
	if (!outofBox)	{
		perAtomtype = lowLowLowFrac * mapvalue[lowLowLowIndex] +
			lowLowHighFrac * mapvalue[lowLowHighIndex] +
			lowHighLowFrac * mapvalue[lowHighLowIndex] +
			lowHighHighFrac * mapvalue[lowHighHighIndex] +
			highLowLowFrac * mapvalue[highLowLowIndex] +
			highLowHighFrac * mapvalue[highLowHighIndex] +
			highHighLowFrac * mapvalue[highHighLowIndex] +
			highHighHighFrac * mapvalue[highHighHighIndex];
		eStatic = charge * (lowLowLowFrac * emapvalues[lowLowLowIndex] +
			lowLowHighFrac * emapvalues[lowLowHighIndex] +
			lowHighLowFrac * emapvalues[lowHighLowIndex] +
			lowHighHighFrac * emapvalues[lowHighHighIndex] +
			highLowLowFrac * emapvalues[highLowLowIndex] +
			highLowHighFrac * emapvalues[highLowHighIndex] +
			highHighLowFrac * emapvalues[highHighLowIndex] +
			highHighHighFrac * emapvalues[highHighHighIndex]);
		deSolv = abscharge * (lowLowLowFrac * dmapvalues[lowLowLowIndex] +
			lowLowHighFrac * dmapvalues[lowLowHighIndex] +
			lowHighLowFrac * dmapvalues[lowHighLowIndex] +
			lowHighHighFrac * dmapvalues[lowHighHighIndex] +
			highLowLowFrac * dmapvalues[highLowLowIndex] +
			highLowHighFrac * dmapvalues[highLowHighIndex] +
			highHighLowFrac * dmapvalues[highHighLowIndex] +
			highHighHighFrac * dmapvalues[highHighHighIndex]);
		//if (deSolv < 0.001)
		//	deSolv = -0.1 * abscharge;
		erg = perAtomtype + deSolv + eStatic;
	}

	//fprintf(stderr, "index %i exenergy %g atom %g estatic %g des %g \n", highHighHighIndex, erg, perAtomtype, deSolv, eStatic);
	if (erg>1000000|| erg<-1000000){
		//0 / 0;
		//fprintf(stderr, "index %i exenergy %g atom %g estatic %g des %g \n", i, erg, perAtomtype, deSolv, eStatic);
		if (perAtomtype != 0) {
			fprintf(stderr, "exenergy %g atom %g estatic %g des %g %g %g %g %g \n",
				mapvalue[lowLowLowIndex], mapvalue[lowLowHighIndex],
				mapvalue[lowHighLowIndex], mapvalue[lowHighHighIndex],
				mapvalue[highLowLowIndex], mapvalue[highLowHighIndex],
				mapvalue[highHighLowIndex], mapvalue[highHighHighIndex]);
			fprintf(stderr, "X %g Y %g Z %g \n", exactGridX, exactGridY, exactGridZ);
			fprintf(stderr, "X %g Y %g Z %g \n", X, Y, Z);
			stop("baddd");
		}

		//fprintf(stderr, "X %g Y %g Z %g \n", exactGridX, exactGridY, exactGridZ);
		//fprintf(stderr, "X %g Y %g Z %g \n", X, Y, Z);
		//stop("baddd");
		//0/0;
		
	}
	//fprintf(stderr, "index %i exenergy %g atom %g estatic %g des %g \n", highHighHighIndex, erg, perAtomtype, deSolv, eStatic);
	return erg;
}

static int indMoved(int ind, int start, int end){
	if (start>=end){
		if (ind>end && ind<start)
			return 0;
		else
			return 1;
	} else {
		if (ind>=start && ind<=end)
			return 1;
		else
			return 0;
	}
}


void ADenergyNoClash(double* ADEnergies, int start, int end, Chain *chain, Chaint *chaint, model_params *mod_params, int mod)
{
  //MS what does mod do in ADenergyNoClash
  /* only calculate for constrained amino acids */
	/* TODO: add constraint type other than 1 */
	//if ((mod_params->external_potential_type != 1 && mod_params->external_potential_type != 3) || !(a->etc & CONSTRAINED)) return 0.0;
	/* C-O-M or n, ca, c */
	//gridmap_initialise();

        // MS that limits the size of peptide to 30*NAA atoms
	double coordsSet[30 * chain->NAA];
	//double *currgridmapvalues = malloc(NX*NY*NZ * sizeof(double));

	//double *coordsSet = malloc(21 * chain->NAA * sizeof(double));

	int ind = 0;
	int numRand = 1;
	int numDir = 1;
	if (mod == 1) {
		numRand = 3;
		numDir = 2;
	}

	AA* a;
	int i = 0; int j = 0; int m = 0;
	for (i = 0; i <= 30 * chain->NAA -1; i++){
		coordsSet[i] = 9999.;
	}

	int linked;
	if (end > chain->NAA-1) linked = 1;
	else linked = 0;

	for (i = 1; i <= chain->NAA-1; i++){
		if (chaint!=NULL && indMoved(i, start, (end-1)%(chain->NAA-1)+1 )) {
			a = chaint->aat + i;
		}
		else {
			a = chain->aa + i;
		}
		if (!indMoved(i, start, (end-1)%(chain->NAA-1)+1 )) {
			if (a->id != 'G'){
				coordsSet[ind]=a->cb[0]; ind++;
				coordsSet[ind]=a->cb[1]; ind++;
				coordsSet[ind]=a->cb[2]; ind++;
			}
			if (a->etc & G__){
				coordsSet[ind]=a->g[0]; ind++;
				coordsSet[ind]=a->g[1]; ind++;
				coordsSet[ind]=a->g[2]; ind++;
			}
			if (a->etc & G2_){
				coordsSet[ind]=a->g2[0]; ind++;
				coordsSet[ind]=a->g2[1]; ind++;
				coordsSet[ind]=a->g2[2]; ind++;
			}
		}
		coordsSet[ind]=a->ca[0]; ind++;
		coordsSet[ind]=a->ca[1]; ind++;
		coordsSet[ind]=a->ca[2]; ind++;
		coordsSet[ind]=a->c[0]; ind++;
		coordsSet[ind]=a->c[1]; ind++;
		coordsSet[ind]=a->c[2]; ind++;
		coordsSet[ind]=a->n[0]; ind++;
		coordsSet[ind]=a->n[1]; ind++;
		coordsSet[ind]=a->n[2]; ind++;
		coordsSet[ind]=a->o[0]; ind++;
		coordsSet[ind]=a->o[1]; ind++;
		coordsSet[ind]=a->o[2]; ind++;
	}
	//fprintf(stderr, "count %d %d %d\n", ind,start,end);
	int notmovedind = ind;
	double sideChainEnergy = 0.0;
	double erg = 0.0;
	double exC = 0.0, exCa = 0.0, exN = 0.0, exO = 0.0, exCb = 0.0, exH = 0.0;
	double CCharge = 0.241, CaCharge = 0.186, NCharge = -0.346, OCharge = -0.271, CbCharge = 0.050, HCharge = 0.163;
	//for (int i =0; i< ind; i++) fprintf(stderr, "count C %g \n", coordsSet[i]);
	//double *energiesforward = malloc((end-start+1) * sizeof(double));
	//double *energiesbackward = malloc((end-start+1) * sizeof(double));
	double energiesforward[end-start+1];
	double energiesbackward[end-start+1];
	for(int m=start; m<=end; m++){
		energiesforward[m-start] = 99999.0;
		energiesbackward[m-start] = 99999.0;
		//fprintf(stderr, "bb Energy %g %g\n", energiesforward[m-start],energiesbackward[m-start]);
	}
	int direction = 1;
	if (mod == 0) direction = rand()%100<50 ? 1 : 0;



	for (m=0; m<numDir; m++) {
		for (int ii = notmovedind; ii <= 30 * chain->NAA -1; ii++){
			coordsSet[ii] = 9999.;
		}
		
		ind = notmovedind;
		for (j = start; j <= end; j++) {
			CCharge = 0.241, CaCharge = 0.186, NCharge = -0.346, OCharge = -0.271, CbCharge = 0.050, HCharge = 0.163;
			if ((mod == 1 && m == 0) || direction == 0) 
				i = j;
			else
				i = end - j + start;

			if (chaint!=NULL)
				a = chaint->aat + (1 + (i-1)%(chain->NAA-1));
			else
				a = chain->aa + (1 + (i-1)%(chain->NAA-1));
			//if (chaint!=NULL)
			//	a = chaint->aat + i;
			//else
			//	a = chain->aa + i;

			/* element types are 0:C, 1:N, 2:O, 3:H, 4:S, 5:CA, 6:NA           */
			
			//fprintf(stderr, "energies C %g CA %g N %g O %g \n", a->c[0], a->c[1], a->c[2], exO);

			if (a->id == 'G') {
				CaCharge = 0.218;
			}
			else if (a->id == 'S') {
				CaCharge = 0.219;
				CbCharge = 0.199;
			}
			else if (a->id == 'P') {
				CaCharge = 0.165;
				CbCharge = 0.034;
				NCharge = -0.29;
			}
			else if (a->id == 'C') {
				CbCharge = 0.120;
			}
			else if (a->id == 'T' || a->id == 'D' || a->id == 'N') {
				CbCharge = 0.146;
			}

			if (a->num == 1) {
				NCharge = -0.06;
				HCharge = 0.275;
				CCharge = 0.21;
				//CaCharge = 0.28;
			}

			if (a->num == chain->NAA - 1) {
				//NCharge = -0.06;
				OCharge = -0.65;
				//CCharge = 0.484;
				CCharge = 0.21;
			}

			sideChainEnergy = 0.0;
			erg = 0.0;
			//exC = 0.0, exCa = 0.0, exN = 0.0, exO = 0.0, exCb = 0.0, exH = 0.0;
			if (a->id != 'P') {
				exH = gridenergy(a->h[0], a->h[1], a->h[2], 3, HCharge);
			}
			exC = gridenergy(a->c[0], a->c[1], a->c[2], 0, CCharge);
			//fprintf(stderr, "energies C %g CA %g N %g O %g \n", exC, exCa, exN, exO);
			exCa = gridenergy(a->ca[0], a->ca[1], a->ca[2], 0, CaCharge);
			//fprintf(stderr, "energies C %g CA %g N %g O %g \n", exC, exCa, exN, exO);
			//exH = gridenergy(a->h[0], a->h[1], a->h[2], 3, HCharge);
			if (a->id != 'G') {
				exCb = gridenergy(a->cb[0], a->cb[1], a->cb[2], 0, CbCharge);
			}
			//fprintf(stderr, "energies C %g CA %g N %g O %g \n", exC, exCa, exN, exO);
			exN = gridenergy(a->n[0], a->n[1], a->n[2], 1, NCharge);
			//fprintf(stderr, "energies C %g CA %g N %g O %g \n", exC, exCa, exN, exO);
			exO = gridenergy(a->o[0], a->o[1], a->o[2], 2, OCharge);
			//fprintf(stderr, "energies C %g CA %g N %g O %g \n", exC, exCa, exN, exO);
			erg = (exC + exCa + exH + exN + exO + exCb);
			
			if (erg > 10000000 || erg < -10000000) {
				fprintf(stderr, "energies C %g CA %g N %g O %g Cb %g H %g", exC, exCa, exN, exO, exCb, exH);
				stop("Grid energy exceeds limits, something wrong!");
			}

				
			/*Here is a hack, external_r0[0] term 1.x indicate to reconstruct full-atom sidechain score grid energy */
			if ((int) mod_params->external_r0[0] == 1) {
			  /* sideChainEnergy = scoreSideChainNoClashNew(a, coordsSet, ind, numRand); */
			    if (a->id=='P') {
			        sideChainEnergy = scoreSideChainNew(a, 1);
			    } else if (a->sideChainTemplateIndex==-1) {
			        sideChainEnergy = 0.0; // A or G
			    } else {
			        sideChainEnergy = scoreSideChainNoClashNew(a, coordsSet, ind, numRand);
			    }
			  /* switch (a->id) */
			  /* 	{ */
			  /* 	case 'I': */
			  /* 		//sideChainEnergy = gridenergy(a->g2[0], a->g2[1], a->g2[2], 0, 0.012) + gridenergy(a->g[0], a->g[1], a->g[2], 0, 0.012); */
			  /* 		sideChainEnergy = scoreSideChainNoClash(ILE.nbRot, ILE.nbAtoms, ILE.charges, ILE.atypes, ILE.coords, a, coordsSet, ind, numRand); */
			  /* 		break; */
			  /* 	case 'L': */
			  /* 		sideChainEnergy = scoreSideChainNoClash(LEU.nbRot, LEU.nbAtoms, LEU.charges, LEU.atypes, LEU.coords, a, coordsSet, ind, numRand); */
			  /* 		break; */
			  /* 	case 'P': */
			  /* 		sideChainEnergy = scoreSideChain(PRO.nbRot, PRO.nbAtoms, PRO.charges, PRO.atypes, PRO.coords, a, 1); */
			  /* 		break; */
			  /* 	case 'V': */
			  /* 		sideChainEnergy = scoreSideChainNoClash(VAL.nbRot, VAL.nbAtoms, VAL.charges, VAL.atypes, VAL.coords, a, coordsSet, ind, numRand); */
			  /* 		//sideChainEnergy = gridenergy(a->g2[0], a->g2[1], a->g2[2], 0, 0.012) + gridenergy(a->g[0], a->g[1], a->g[2], 0, 0.012); */
			  /* 		break; */
			  /* 	case 'F': */
			  /* 		sideChainEnergy = scoreSideChainNoClash(PHE.nbRot, PHE.nbAtoms, PHE.charges, PHE.atypes, PHE.coords, a, coordsSet, ind, numRand); */
			  /* 		break; */
			  /* 	case 'W': */
			  /* 		sideChainEnergy = scoreSideChainNoClash(TRP.nbRot, TRP.nbAtoms, TRP.charges, TRP.atypes, TRP.coords, a, coordsSet, ind, numRand); */
			  /* 		break; */
			  /* 	case 'Y': */
			  /* 		sideChainEnergy = scoreSideChainNoClash(TYR.nbRot, TYR.nbAtoms, TYR.charges, TYR.atypes, TYR.coords, a, coordsSet, ind, numRand); */
			  /* 		break; */
			  /* 	case 'D': */
			  /* 		sideChainEnergy = scoreSideChainNoClash(ASP.nbRot, ASP.nbAtoms, ASP.charges, ASP.atypes, ASP.coords, a, coordsSet, ind, numRand); */
			  /* 		break; */
			  /* 	case 'E': */
			  /* 		sideChainEnergy = scoreSideChainNoClash(GLU.nbRot, GLU.nbAtoms, GLU.charges, GLU.atypes, GLU.coords, a, coordsSet, ind, numRand); */
			  /* 		break; */
			  /* 	case 'R': */
			  /* 		sideChainEnergy = scoreSideChainNoClash(ARG.nbRot, ARG.nbAtoms, ARG.charges, ARG.atypes, ARG.coords, a, coordsSet, ind, numRand); */
			  /* 		break; */
			  /* 	case 'H': */
			  /* 		sideChainEnergy = scoreSideChainNoClash(HIS.nbRot, HIS.nbAtoms, HIS.charges, HIS.atypes, HIS.coords, a, coordsSet, ind, numRand); */
			  /* 		break; */
			  /* 	case 'K': */
			  /* 		sideChainEnergy = scoreSideChainNoClash(LYS.nbRot, LYS.nbAtoms, LYS.charges, LYS.atypes, LYS.coords, a, coordsSet, ind, numRand); */
			  /* 		break; */
			  /* 	case 'S': */
			  /* 		sideChainEnergy = scoreSideChainNoClash(SER.nbRot, SER.nbAtoms, SER.charges, SER.atypes, SER.coords, a, coordsSet, ind, numRand); */
			  /* 		//sideChainEnergy = gridenergy(a->g[0], a->g[1], a->g[2], 2, -0.398); */
			  /* 		break; */
			  /* 	case 'T': */
			  /* 		sideChainEnergy = scoreSideChainNoClash(THR.nbRot, THR.nbAtoms, THR.charges, THR.atypes, THR.coords, a, coordsSet, ind, numRand); */
			  /* 		//sideChainEnergy = gridenergy(a->g2[0], a->g2[1], a->g2[2], 2, -0.393) +  gridenergy(a->g[0], a->g[1], a->g[2], 0, 0.042); */
			  /* 		break; */
			  /* 	case 'C': */
			  /* 		sideChainEnergy = scoreSideChainNoClash(CYS.nbRot, CYS.nbAtoms, CYS.charges, CYS.atypes, CYS.coords, a, coordsSet, ind, numRand); */
			  /* 		//sideChainEnergy = gridenergy(a->g[0], a->g[1], a->g[2], 4, -0.095); */
			  /* 		break; */
			  /* 	case 'M': */
			  /* 		sideChainEnergy = scoreSideChainNoClash(MET.nbRot, MET.nbAtoms, MET.charges, MET.atypes, MET.coords, a, coordsSet, ind, numRand); */
			  /* 		break; */
			  /* 	case 'N': */
			  /* 		sideChainEnergy = scoreSideChainNoClash(ASN.nbRot, ASN.nbAtoms, ASN.charges, ASN.atypes, ASN.coords, a, coordsSet, ind, numRand); */
			  /* 		break; */
			  /* 	case 'Q': */
			  /* 		sideChainEnergy = scoreSideChainNoClash(GLN.nbRot, GLN.nbAtoms, GLN.charges, GLN.atypes, GLN.coords, a, coordsSet, ind, numRand); */
			  /* 		break; */
			  /* 	default: */
			  /* 		break; */
			  /* 	} */
			}
			erg += sideChainEnergy;
			// AD energy is in kcal/mol, scale down kcal/mol to RT!
			erg = erg / 0.59219;
			//if (chaint!=NULL)
			//	chaint->Ergt(0, i) = erg;
			//else
			//	chain->Erg(0, i) = erg;
			if ((mod == 1 && m == 0) || direction == 0) 
				energiesforward[i-start] = erg;
			else
				energiesbackward[i-start] = erg;
			
			if (a->id != 'G'){
				coordsSet[ind]=a->cb[0]; ind++;
				coordsSet[ind]=a->cb[1]; ind++;
				coordsSet[ind]=a->cb[2]; ind++;
			}
			
			if (a->etc & G__){
				coordsSet[ind]=a->g[0]; ind++;
				coordsSet[ind]=a->g[1]; ind++;
				coordsSet[ind]=a->g[2]; ind++;
			}
			if (a->etc & G2_){
				coordsSet[ind]=a->g2[0]; ind++;
				coordsSet[ind]=a->g2[1]; ind++;
				coordsSet[ind]=a->g2[2]; ind++;
			}

			//if (mod == 1){
			//	fprintf(stderr, "Energy  i %d j %d type %c %g start %d end %d mod %d ind %d\n", i, j, a->id, erg,start,end,mod,ind);
			//}
		}


	}
	double totE1 = 0.;
	double totE2 = 0.;
	for(m=start; m<=end; m++){
		totE1 += energiesforward[m-start];
		totE2 += energiesbackward[m-start];
		//if (mod == 1)
		//	fprintf(stderr, "bb Energy %g %g\n", energiesforward[m-start],energiesbackward[m-start]);
	}
	//if (mod == 1)
	//	fprintf(stderr, "bb Energy %g %g\n", totE1,totE2);
	//free(coordsSet);
	//double* ADenergies;
	if (totE1 < totE2){
		//ADenergies = energiesforward;
		for(m=start; m<=end; m++){
			
			ADEnergies[m-start] = energiesforward[m-start];
			//if (mod == 1){
			//	fprintf(stderr, "Forward Energy  i %d %g \n", m, ADEnergies[m-start]);
			//}
		}
		//free(energiesbackward);
	} else {
		for(m=start; m<=end; m++){
			ADEnergies[m-start] = energiesbackward[m-start];
			//if (mod == 1){
			//	fprintf(stderr, "Backward Energy  i %d %g \n", m, ADEnergies[m-start]);
			//}
		}
		//ADenergies = energiesbackward;
		//free(energiesforward);
	}




	//free(coordsSet);
	//ADenergies = energiesforward;
	//for(int m=start; m<=end; m++){
	//	fprintf(stderr, "bb Energy %g %d %d\n", energiesforward[m-start],start,end);
	//	fprintf(stderr, "bb Energy %g %d %d\n", energiesbackward[m-start],start,end);
	//}
	//return ADenergies;
	//for(int m=start; m<=end; m++){
	//	fprintf(stderr, "bb Energy %g %d %d\n", energiesforward[m-start],start,end);
	//
	//}

	//return energiesforward;
	//free(coordsSet);
}


/* external potential depending on atomic position */
double external(AA *a, model_params *mod_params, vector molcom)
{

	/* only calculate for constrained amino acids */
	/* TODO: add constraint type other than 1 */
	if ((mod_params->external_potential_type != 1 && mod_params->external_potential_type != 3) || !(a->etc & CONSTRAINED)) return 0.0;
	/* C-O-M or n, ca, c */
	//gridmap_initialise();
	vector com;	/* N-Ca and Ca-C bonds */
	double erg = 0.0;
	double dr;

	add(com, a->ca, a->n);
	add(com, com, a->c);
	scale(com,1.0/3.0, com);

	//fprintf(stderr, "calculating constraint on amino acid %d", a->num);
	/* constraining to the z axis using a harmonic potential
		E =	k * ( sqrt(x^2+y^2) - r0 )^2	if x^2+y^2 > r0^2,
			0				otherwise.  */
	if (mod_params->external_potential_type == 1) {

		//fprintf(stderr,"calculating constraint on amino acid %d", a->num);
		//fprintf(stderr, "calculating constraint on amino acid %g %g", com[0], com[1]);
    		//fprintf(stderr," (etc: %x), %x\n",a->etc,(a->etc & CONSTRAINED));
		//Original CRANKITE
		double dr2 = com[0]*com[0] + com[1]*com[1]; //distance^2 from (0,0) in the (x,y) plane
		if (dr2 > mod_params->external_r0[0]*mod_params->external_r0[0]) {
			dr = sqrt(dr2) - mod_params->external_r0[0];
			erg += mod_params->external_k[0] * dr * dr;
		}


	} else if (mod_params->external_potential_type == 2) {
		stop("unimplemented type 2");
		if (mod_params->external_direction[0] == EXTERNAL_POSITIVE || mod_params->external_direction[0] == EXTERNAL_POSNEG ) {
			dr = com[0] - mod_params->external_r0[0];
			if (dr > 0.0) erg += mod_params->external_k[0] * dr * dr;
		}
		if (mod_params->external_direction[0] == EXTERNAL_NEGATIVE || mod_params->external_direction[0] == EXTERNAL_POSNEG ) {
			dr = com[0] + mod_params->external_r0[0];
			if (dr < 0) erg += mod_params->external_k[0] * dr * dr;
		}
		if (mod_params->external_direction[1] == EXTERNAL_POSITIVE || mod_params->external_direction[1] == EXTERNAL_POSNEG ) {
			dr = com[1] - mod_params->external_r0[1];
			if (dr < 0) erg += mod_params->external_k[1] * dr * dr;
		}
		if (mod_params->external_direction[1] == EXTERNAL_NEGATIVE || mod_params->external_direction[1] == EXTERNAL_POSNEG ) {
			dr = com[1] + mod_params->external_r0[1];
			if (dr < 0) erg += mod_params->external_k[1] * dr * dr;
		}
		if (mod_params->external_direction[2] == EXTERNAL_POSITIVE || mod_params->external_direction[2] == EXTERNAL_POSNEG ) {
			dr = com[2] - mod_params->external_r0[2];
			if (dr < 0) erg += mod_params->external_k[2] * dr * dr;
		}
		if (mod_params->external_direction[2] == EXTERNAL_NEGATIVE || mod_params->external_direction[2] == EXTERNAL_POSNEG ) {
			dr = com[2] + mod_params->external_r0[2];
			if (dr < 0) erg += mod_params->external_k[2] * dr * dr;
		} 
	/* conical potential around the z axis.
		E =	k * ( sqrt(x^2+y^2) - r0 * (1 - (z-COM)/z_tip) )^2	if z/z_tip in [0,1] and x^2+y^2 > r0^2,
			0						otherwise.
	COM is the system centre-of-mass.  */
	} else if (mod_params->external_potential_type == 3) {
		//stop("unimplemented type 3");
		double zratio = (com[2] - molcom[2]) / mod_params->external_ztip;
		if (zratio < 0) zratio = 0; /* we are on the other side, cylindrical potential */
		if (zratio > 1) zratio = 1;
		double dr2 = com[0]*com[0] + com[1]*com[1]; //distance^2 from (0,0) in the (x,y) plane
		if (dr2 > (1.0 - zratio) * (1.0 - zratio) * mod_params->external_r0[0]*mod_params->external_r0[0]) {
			dr = sqrt(dr2) - (1.0 - zratio) * mod_params->external_r0[0];
			erg += mod_params->external_k[0] * dr * dr;
			//fprintf(stderr,"%d com=(%g,%g,%g) molcom=(%g,%g,%g) zratio=%g Ee=%g*%g^2=%g\n",
			//	a->num, com[0], com[1], com[2], molcom[0], molcom[1], molcom[2],
			//	zratio, mod_params->external_k[0], dr, mod_params->external_k[0] * dr * dr);
		}
	}
	//fprintf(stderr, "external energy2 %g \n", erg);
	return erg;
}





/* external potential depending on atomic position */
double external2(AA *a, model_params *mod_params, vector molcom)
{
	
	/* only calculate for constrained amino acids */
	/* TODO: add constraint type other than 1 */
	if ((mod_params->external_potential_type2 != 1 && mod_params->external_potential_type2 != 3) || !(a->etc & CONSTRAINED2)) return 0.0;

	/* C-O-M or n, ca, c */
	vector com;	/* N-Ca and Ca-C bonds */
	double erg = 0.0;
	double dr;

	add(com, a->ca, a->n);
	add(com, com, a->c);
	scale(com,1.0/3.0, com);

	if (mod_params->external_potential_type2 == 1) {

		//fprintf(stderr,"calculating constraint on amino acid %d", a->num);
    		//fprintf(stderr," (etc: %x), %x\n",a->etc,(a->etc & CONSTRAINED));

		double dr2 = com[0]*com[0] + com[1]*com[1]; //distance^2 from (0,0) in the (x,y) plane
		if (dr2 > mod_params->external_r0[0]*mod_params->external_r0[0]) {
			dr = sqrt(dr2) - mod_params->external_r0[0];
			erg += mod_params->external_k[0] * dr * dr;
		}
//		//x+
//		dr = com[0] - mod_params->external_r02[0];
//		if (dr > 0.0) erg += mod_params->external_k2[0] * dr * dr;
//		//x-
//		dr = com[0] + mod_params->external_r02[0];
//		if (dr < 0.0) erg += mod_params->external_k2[0] * dr * dr;
//		//y+
//		dr = com[1] - mod_params->external_r02[1];
//		if (dr > 0.0) erg += mod_params->external_k2[1] * dr * dr;
//		//y-
//		dr = com[1] + mod_params->external_r02[1];
//		if (dr < 0.0) erg += mod_params->external_k2[1] * dr * dr;

	} else if (mod_params->external_potential_type2 == 2) {
		stop("unimplemented type 2");
		if (mod_params->external_direction2[0] == EXTERNAL_POSITIVE || mod_params->external_direction2[0] == EXTERNAL_POSNEG ) {
			dr = com[0] - mod_params->external_r02[0];
			if (dr > 0.0) erg += mod_params->external_k2[0] * dr * dr;
		}
		if (mod_params->external_direction2[0] == EXTERNAL_NEGATIVE || mod_params->external_direction2[0] == EXTERNAL_POSNEG ) {
			dr = com[0] + mod_params->external_r02[0];
			if (dr < 0) erg += mod_params->external_k2[0] * dr * dr;
		}
		if (mod_params->external_direction2[1] == EXTERNAL_POSITIVE || mod_params->external_direction2[1] == EXTERNAL_POSNEG ) {
			dr = com[1] - mod_params->external_r02[1];
			if (dr < 0) erg += mod_params->external_k2[1] * dr * dr;
		}
		if (mod_params->external_direction2[1] == EXTERNAL_NEGATIVE || mod_params->external_direction2[1] == EXTERNAL_POSNEG ) {
			dr = com[1] + mod_params->external_r02[1];
			if (dr < 0) erg += mod_params->external_k2[1] * dr * dr;
		}
		if (mod_params->external_direction2[2] == EXTERNAL_POSITIVE || mod_params->external_direction2[2] == EXTERNAL_POSNEG ) {
			dr = com[2] - mod_params->external_r02[2];
			if (dr < 0) erg += mod_params->external_k2[2] * dr * dr;
		}
		if (mod_params->external_direction2[2] == EXTERNAL_NEGATIVE || mod_params->external_direction2[2] == EXTERNAL_POSNEG ) {
			dr = com[2] + mod_params->external_r02[2];
			if (dr < 0) erg += mod_params->external_k2[2] * dr * dr;
		}
	/* conical potential around the z axis.
		E =	k * ( sqrt(x^2+y^2) - r0 * (1 - (z-COM)/z_tip) )^2	if z/z_tip in [0,1] and x^2+y^2 > r0^2,
			0						otherwise.
	COM is the system centre-of-mass.  */
	} else if (mod_params->external_potential_type2 == 3) {
		//stop("unimplemented type 3");
		double zratio = (com[2] - molcom[2]) / mod_params->external_ztip2;
		if (zratio < 0) zratio = 0; /* we are on the other side, cylindrical potential */
		if (zratio > 1) zratio = 1;
		double dr2 = com[0]*com[0] + com[1]*com[1]; //distance^2 from (0,0) in the (x,y) plane
		if (dr2 > (1.0 - zratio) * (1.0 - zratio) * mod_params->external_r02[0]*mod_params->external_r02[0]) {
			dr = sqrt(dr2) - (1.0 - zratio) * mod_params->external_r02[0];
			erg += mod_params->external_k2[0] * dr * dr;
		}
	}
	return erg;
}





/***********************************************************/
/****          ENERGY CONTRIBUTIONS  SUMMED UP          ****/
/***********************************************************/

/* internal amino acid interactions */
/* or external potential depending on amino acid position */
double energy1(AA *a, model_params *mod_params)
{
	double retval = 0.0;

	/* internal potential */
	retval += stress(a, mod_params) + clash(a, mod_params);

	//MOVED TO GLOBAL_ENERGY
	///* external potential */
	//retval += external(a, mod_params);
	//retval += external2(a, mod_params);

	return retval;
}

/* interactions between two amino acids */
double energy2cyclic(Biasmap *biasmap, AA *a,  AA *b, model_params *mod_params)
{
	double retval = 0.0;

	/* Go-type bias potential */
	
	if (biasmap->distb && Distb(a->num, b->num) != 0.0)		
		retval += bias(biasmap, a, b, mod_params);
	//fprintf(stderr,"e21 %g\n",retval);

    
	retval += hydrophobic(biasmap,a,b, mod_params);
	//fprintf(stderr,"e22 %d %d %g\n",a->num,b->num,hydrophobic(biasmap,a,b, mod_params));
	if (mod_params->use_gamma_atoms != NO_GAMMA) {
		retval += electrostatic(biasmap,a,b, mod_params);
		//if (electrostatic(biasmap,a,b, mod_params)>0.02) {
		//	fprintf(stderr,"e23 %d %d %g\n",a->num,b->num,electrostatic(biasmap,a,b, mod_params));
		//}
		retval += sidechain_hbond(biasmap,a,b, mod_params);
		//fprintf(stderr,"e24 %g\n",retval);
	}

	retval += exclude_neighbor(b, a, mod_params) + hbond(biasmap, b, a, mod_params) + proline(b, a);


	return retval;
}

/* interactions between two amino acids */
double energy2(Biasmap *biasmap, AA *a,  AA *b, model_params *mod_params)
{
	double d2, retval = 0.0;

	/* Go-type bias potential */
	
	if (biasmap->distb && Distb(a->num, b->num) != 0.0)		
		retval += bias(biasmap, a, b, mod_params);
	//fprintf(stderr,"e21 %g\n",retval);

    
	retval += hydrophobic(biasmap,a,b, mod_params);
	//fprintf(stderr,"e22 %d %d %g\n",a->num,b->num,hydrophobic(biasmap,a,b, mod_params));
	if (mod_params->use_gamma_atoms != NO_GAMMA) {
		retval += electrostatic(biasmap,a,b, mod_params);
		//if (electrostatic(biasmap,a,b, mod_params)>0.02) {
		//	fprintf(stderr,"e23 %d %d %g\n",a->num,b->num,electrostatic(biasmap,a,b, mod_params));
		//}
		retval += sidechain_hbond(biasmap,a,b, mod_params);
		//fprintf(stderr,"e24 %g\n",retval);
	}
	
	double CaDistance = distance(a->ca, b->ca);
	int seqdist;
	if (a->chainid == b->chainid)
		seqdist = b->num - a->num;
	else
		seqdist = 1000 * abs(b->chainid - a->chainid);

	switch ( seqdist) {
	case 1:
		retval += exclude_neighbor(a, b, mod_params) + hbond(biasmap,a, b, mod_params) + proline(a, b);
		retval += (sqrt(CaDistance) - 3.819)*(sqrt(CaDistance) - 3.819);
		//fprintf(stderr,"e25a %d %d %g\n",a->num,b->num,hbond(biasmap,a, b, mod_params));
		//fprintf(stderr,"e25a %g\n",retval);
		break;
	case -1:
		retval += exclude_neighbor(b, a, mod_params) + hbond(biasmap, b, a, mod_params) + proline(b, a);
                retval += (sqrt(CaDistance) - 3.819)*(sqrt(CaDistance) - 3.819);
		//fprintf(stderr,"e25b %d %d %g\n",a->num,b->num,hbond(biasmap,a, b, mod_params));
		//fprintf(stderr,"e25b %g\n",retval);
		break;
	default:
		d2 = distance(a->ca, b->ca);
		if (d2 < mod_params->vdw_extended_cutoff) {
			retval += exclude(a, b, d2, mod_params);
			if (d2 < hbond_cutoff) {
				retval += hbond(biasmap,a, b, mod_params);
				//fprintf(stderr,"e25c %d %d %g\n",a->num,b->num,hbond(biasmap,a, b, mod_params));
			}
		}
		break;
	}

	return retval;
}


// Gary Hack cyclic peptides type 0: C-N bond, type 1: -S-S- bond to be added if needed
double cyclic_energy(AA *a, AA *b, int type) {
	double ans = 0.;
	if (type == 0) {

		double CaDistance = 0.0;
		double NCDistance = 0.0;
		//double HODistance = 0.0;
		//double NODistance = 0.0;
		//double HCDistance = 0.0;

		//NCDistance = distance(a->n, b->c);
		CaDistance = distance(a->ca, b->ca);
		//HODistance = distance(a->h, b->o);
		//NODistance = distance(a->n, b->o);
		//HCDistance = distance(a->h, b->c);

		//if (1 || CaDistance < 5) ans += 10 * (sqrt(CaDistance) - 3.819);
		if (CaDistance < 5) {
			ans += (sqrt(CaDistance) - 3.819)*(sqrt(CaDistance) - 3.819);
			ans += (sqrt(NCDistance) - 1.345)*(sqrt(NCDistance) - 1.345);
		} else ans += CaDistance;
		//if (NCDistance > 2) ans += (sqrt(NCDistance) - 1.345)*(sqrt(NCDistance) - 1.345);
		//else ans += NCDistance;
		//if (a->id != 'P') ans += (sqrt(HODistance) - 3.13)*(sqrt(HODistance) - 3.13);
		//if (1 || NODistance > 3.5 || NODistance < 1.2) ans += (sqrt(NODistance) - 2.25)*(sqrt(NODistance) - 2.25);
		//if (a->id != 'P') ans += 5 * (sqrt(HCDistance) - 2.02)*(sqrt(HCDistance) - 2.02);
	}
	return ans;
	//return 0.;
}

/*This calculates energy which depends on more than 2 residues, e.g. srgy 
 * residues i: start <= i <= end will be from aat all others from aa */ 
double global_energy(int start, int end, Chain *chain, Chaint *chaint, Biasmap *biasmap, model_params *mod_params){
  
  /* S-S bonds */
  double ans = sbond_energy(start,end,chain,chaint,biasmap,mod_params);

  /* secondary radius of gyration */
  if(mod_params->srgy_param != 0.0) {
    double r_gyr = secondary_radius_of_gyration(start, end, chain, chaint, biasmap, mod_params, CA_, 0);
    if(r_gyr >= mod_params->srgy_offset) {
    	ans += mod_params->srgy_param*(r_gyr-mod_params->srgy_offset)*(r_gyr-mod_params->srgy_offset);
    }
  }

  /* hydrophobic secondary radius of gyration */
  if(mod_params->hphobic_srgy_param != 0.0) {
    double r_gyr = secondary_radius_of_gyration(start, end, chain, chaint, biasmap, mod_params, CB_, 1);
    if(r_gyr >= mod_params->hphobic_srgy_offset) {
        ans += mod_params->hphobic_srgy_param*(r_gyr-mod_params->hphobic_srgy_offset)*(r_gyr-mod_params->hphobic_srgy_offset);
    }
  }


  //fprintf(stderr, "calcuaaalating constraint on amino acidhaha1");
  /* external potentials (may depend on the molecule centre-of-mass) */
  vector mol_com;
  mol_com[0] = mol_com[1] = mol_com[2] = 0.0;
  for (int i = 1; i < chain->NAA; i++){
	vector com;
	if( i <= end && i >= start){ /* changed peptide section */
		add(com, ((chaint->aat) + i)->ca, ((chaint->aat) + i)->n);
		add(com, com, ((chaint->aat) + i)->c);
	} else {
		add(com, ((chain->aa) + i)->ca, ((chain->aa) + i)->n);
		add(com, com, ((chain->aa) + i)->c);
	}
	scale(com,1.0/3.0, com);
	add(mol_com, mol_com, com);
  }
  scale(mol_com, 1.0/(double)(chain->NAA-1), mol_com);
  double test_external = 0.0;
  for (int i = 1; i < chain->NAA; i++){
	double ee;
	if (!chaint) {
		ee = external((chain->aa) + i, mod_params, mol_com);
		//chain->Erg(0, i) = ee;
		//fprintf(stderr,"external energy3 %g\n", ee);
	}
	else if(i <= end && i >= start){ /* changed peptide section */
		ee = external((chaint->aat) + i, mod_params, mol_com);
		//chaint->Ergt(0, i) = ee;
		//fprintf(stderr,"external energy3 %g\n", ee);
	} else {
		//ee = chain->Erg(0, i);
		//chaint->Ergt(0, i) = ee;
		ee = external((chain->aa) + i, mod_params, mol_com);
		//fprintf(stderr, "external energy4 %g\n", ee);
	}
	
	//ee = 10.0;
	ans += ee;
	test_external += ee;
	ans += external2(((chain->aa) + i), mod_params, mol_com);
  }

  return ans;

}

/* The sum of all vdW interactions, including intraresidual interactions (clash)
   and interresidual ones (exclude_neighbour and exclude). */
double all_vdw(Biasmap *biasmap, Chain *chain, model_params *mod_params) {

	int i, j;
	double val = 0.;
	double d2 = 0.;
#ifdef LJ_HBONDED_HARD
	int hbond_proximity;
#endif

	for (i = 1; i < chain->NAA; i++) {
		val += clash(chain->aa + i, mod_params);
		for (j = 1; j < i; j++) {
			switch (i - j) {
			case 1:
				val += exclude_neighbor(chain->aa + i, chain->aa + j, mod_params);
				break;
			default:

#ifdef LJ_HBONDED_HARD
				/* calc if the residues are in contact due to a H-bond */
				/* to exclude the LJ interactions due to H-bond proximity */
				hbond_proximity = 0;
				for (int i1 = i-1; i1<i+2; i1++) {
				    for (int j1 = j-1; j1<j+2; j1++) {
					if (((i1>=1) || (i1<chain->NAA)) && ((j1>=1) || (j1<chain->NAA))) {
					    if (( hstrength(chain->aa[i1].n,chain->aa[i1].h,chain->aa[j1].o,chain->aa[j1].c, mod_params) != 0 ) ||
						( hstrength(chain->aa[j1].n,chain->aa[j1].h,chain->aa[i1].o,chain->aa[i1].c, mod_params) != 0 )) {
						hbond_proximity = 1;
					    }
					}
				    }
				}
#endif
				d2 = distance((chain->aa + i)->ca, (chain->aa + j)->ca);
				if (d2 < mod_params->vdw_extended_cutoff) {
#ifdef LJ_HBONDED_HARD
					if (hbond_proximity) {
					val += exclude_hard(chain->aa + i, chain->aa + j, d2, mod_params,hbond_proximity);
					} else {
#endif
					val += exclude(chain->aa + i, chain->aa + j, d2, mod_params);
#ifdef LJ_HBONDED_HARD
					}
#endif
				}
				break;
			}
		}
	}

	return val;

}


/***********************************************************/
/****                   ENERGY  TESTS                   ****/
/***********************************************************/

/****       FINITE DERIVATIVES WRT THE PARAMETERS       ****/


/* The sum of a function of two amino acids over all amino acid pairs */
static double sumf(Chain* chain, Biasmap* biasmap, double (*fx) ( Biasmap * biasmap,AA *, AA *, model_params *mod_params), model_params *mod_params)
{
	int i, j;
	double val = 0.;

	for (i = 1; i < chain->NAA; i++)
		for (j = 1; j < i; j++)
			val += (*fx) (biasmap, chain->aa + i, chain->aa + j, mod_params);

	return val;
}

/* The sum of a function of one amino acid only over all amino acids */
static double sumf_diag(Chain* chain, Biasmap* biasmap, double (*fx) ( AA *, model_params *mod_params), model_params *mod_params)
{
	int i;
	double val = 0.;

	for (i = 1; i < chain->NAA; i++)
		val += (*fx) (chain->aa + i, mod_params);

	return val;
}

/* The derivative of sumf(fx) using finite difference
   d sumf(fx)(x) / d x ~= 0.5*[sumf(fx)(x+dx)-sumf(fx)(x-dx)]/dx */
static double dfdx(Chain *chain, Biasmap *biasmap, double (*fx) (Biasmap *, AA *, AA *, model_params *mod_params),
		   double *x, double dx, model_params *mod_params)
{
	double val = 0.;

	*x += dx;
	/* update vdw radii of side chains */
	initialize_sidechain_properties(mod_params);
	vdw_param_calculate(mod_params);
	val += sumf(chain, biasmap, fx, mod_params);

	*x -= 2. * dx;
	initialize_sidechain_properties(mod_params);
	vdw_param_calculate(mod_params);
	val -= sumf(chain, biasmap, fx, mod_params);

	*x += dx;
	initialize_sidechain_properties(mod_params);
	vdw_param_calculate(mod_params);

	return 0.5 * val / dx;
}

/* The derivative of sumf_diag(fx) using finite difference
   d sumf_diag(fx)(x) / d x ~= 0.5*[sumf_diag(fx)(x+dx)-sumf_diag(fx)(x-dx)]/dx */
static double intraresidual_dfdx(Chain *chain, Biasmap *biasmap, double (*fx) (AA *, model_params *mod_params),
		   double *x, double dx, model_params *mod_params)
{
	double val = 0.;

	*x += dx;
	initialize_sidechain_properties(mod_params);
	vdw_param_calculate(mod_params);
	val += sumf_diag(chain, biasmap, fx, mod_params);

	*x -= 2. * dx;
	initialize_sidechain_properties(mod_params);
	vdw_param_calculate(mod_params);
	val -= sumf_diag(chain, biasmap, fx, mod_params);

	*x += dx;
	initialize_sidechain_properties(mod_params);
	vdw_param_calculate(mod_params);

	return 0.5 * val / dx;
}

/* The derivative of a global energy function fx(x) using finite difference
   d fx(x) / d x ~= 0.5*[fx(x+dx)-fx(x-dx)]/dx */
static double global_dfdx(Chain *chain, Biasmap *biasmap, double (*fx) (Biasmap *, Chain *, model_params *mod_params),
		   double *x, double dx, model_params *mod_params)
{
	double val = 0.;

	*x += dx;

	if (mod_params->vdw_uniform_depth) {
		mod_params->vdw_depth_cb = mod_params->vdw_depth_c =
					   mod_params->vdw_depth_n =
					   mod_params->vdw_depth_o =
					   mod_params->vdw_depth_s = mod_params->vdw_depth_ca;
	}
	initialize_sidechain_properties(mod_params);
	vdw_param_calculate(mod_params);
	val += (*fx) (biasmap, chain, mod_params);

	*x -= 2. * dx;

	if (mod_params->vdw_uniform_depth) {
		mod_params->vdw_depth_cb = mod_params->vdw_depth_c =
					   mod_params->vdw_depth_n =
					   mod_params->vdw_depth_o =
					   mod_params->vdw_depth_s = mod_params->vdw_depth_ca;
	}
	initialize_sidechain_properties(mod_params);
	vdw_param_calculate(mod_params);
	val -= (*fx) (biasmap, chain, mod_params);

	*x += dx;

	if (mod_params->vdw_uniform_depth) {
		mod_params->vdw_depth_cb = mod_params->vdw_depth_c =
					   mod_params->vdw_depth_n =
					   mod_params->vdw_depth_o =
					   mod_params->vdw_depth_s = mod_params->vdw_depth_ca;
	}
	initialize_sidechain_properties(mod_params);
	vdw_param_calculate(mod_params);

	return 0.5 * val / dx;
}




/* this is an energy-related probe used in contrastive divergence
   avoid %g end %e here, bc hates exponential notation */
void energy_probe_1(Chain* chain, Biasmap *biasmap, simulation_params *sim_params)
{
	double this[36] = { 0., 0., 0., 0., 0.,
				   0., 0., 0., 0., 0.,
				   0., 0., 0., 0., 0.,
				   0., 0., 0., 0., 0.,
				   0., 0., 0., 0., 0.,
				   0., 0., 0., 0., 0.,
				   0., 0., 0., 0., 0.,
				   0. };


	//do not change the original mod_params
	model_params * mod_params = malloc(sizeof(model_params));
	model_params_copy(mod_params,&(sim_params->protein_model));

	/* save the previous results into last */
	for (int i=0; i<36; i++) {
	    sim_params->energy_probe_1_last[i] = sim_params->energy_probe_1_this[i];
	}

	/* hydrogen bonds */
	if (sim_params->energy_probe_1_calc[0])
	    this[0] = dfdx(chain,biasmap,hbond, &(mod_params->hboh2), 0.07, mod_params);
	if (sim_params->energy_probe_1_calc[1])
	    this[1] = dfdx(chain,biasmap,hbond, &(mod_params->hbohn), 0.02, mod_params);
	if (sim_params->energy_probe_1_calc[2])
	    this[2] = dfdx(chain,biasmap,hbond, &(mod_params->hbcoh), 0.02, mod_params);
	if (sim_params->energy_probe_1_calc[3])
	    this[3] = dfdx(chain,biasmap,hbond, &(mod_params->hbs), 0.01, mod_params);

	/* bias potential */
	if (biasmap->distb == NULL)
		goto out;
	/* gradient of energy derivative (likelihood) */
	if (sim_params->energy_probe_1_calc[4])
	    this[4] = dfdx(chain,biasmap,bias, &(mod_params->bias_eta_beta), 0.01, mod_params);
	if (sim_params->energy_probe_1_calc[5])
	    this[5] = dfdx(chain,biasmap,bias, &(mod_params->bias_eta_alpha), 0.01, mod_params);
	if (sim_params->energy_probe_1_calc[6])
	    this[6] = dfdx(chain,biasmap,bias, &(mod_params->bias_kappa_alpha_3), 0.01, mod_params);
	if (sim_params->energy_probe_1_calc[7])
	    this[7] = dfdx(chain,biasmap,bias, &(mod_params->bias_kappa_alpha_4), 0.01, mod_params);
	if (sim_params->energy_probe_1_calc[8])
	    this[8] = dfdx(chain,biasmap,bias, &(mod_params->bias_kappa_beta), 0.01, mod_params);

	/* hydrophobicity */
	if (sim_params->energy_probe_1_calc[9])
	    this[9] = dfdx(chain,biasmap,hydrophobic, &(mod_params->kauzmann_param), 0.01, mod_params);
	if (sim_params->energy_probe_1_calc[10])
	    this[10] = dfdx(chain,biasmap,hydrophobic, &(mod_params->hydrophobic_cutoff_range), 0.07, mod_params);
    
	/* electrostatics */
	if (sim_params->energy_probe_1_calc[11])
	    this[11] = dfdx(chain,biasmap,electrostatic, &(mod_params->recip_dielectric_param), 0.01, mod_params);

	/* bias potential */
	if (sim_params->energy_probe_1_calc[12])
	    this[12] = dfdx(chain,biasmap,bias,&(mod_params->bias_r_alpha),0.01, mod_params);
	if (sim_params->energy_probe_1_calc[13])
	    this[13] = dfdx(chain,biasmap,bias,&(mod_params->bias_r_beta),0.01, mod_params);

	/* electrostatics */
	if (sim_params->energy_probe_1_calc[14])
	    this[14] = dfdx(chain,biasmap,electrostatic, &(mod_params->debye_length_param), 0.01, mod_params);

	/* side chain hydrogen bonds */
	if (sim_params->energy_probe_1_calc[15])
	    this[15] = dfdx(chain,biasmap,sidechain_hbond,&(mod_params->sidechain_hbond_strength_s2b),0.01, mod_params);
	if (sim_params->energy_probe_1_calc[16])
	    this[16] = dfdx(chain,biasmap,sidechain_hbond,&(mod_params->sidechain_hbond_strength_b2s),0.01, mod_params);
	if (sim_params->energy_probe_1_calc[17])
	    this[17] = dfdx(chain,biasmap,sidechain_hbond,&(mod_params->sidechain_hbond_strength_s2s),0.01, mod_params);
	if (sim_params->energy_probe_1_calc[18])
	    this[18] = dfdx(chain,biasmap,sidechain_hbond,&(mod_params->sidechain_hbond_angle_cutoff),0.05, mod_params);

	/* global energy */
    double val = 0.;
    double dx = 0.09;
	if (sim_params->energy_probe_1_calc[19]) {
	    mod_params->srgy_param += dx;
	    val += global_energy(0,0,chain,NULL,biasmap, mod_params);
	    mod_params->srgy_param -= 2. * dx;
	    val -= global_energy(0,0,chain,NULL,biasmap, mod_params);
	    mod_params->srgy_param += dx;

	    this[19] = 0.5 * val / dx;
	}

	if (sim_params->energy_probe_1_calc[20]) {
	    val = 0;
	    mod_params->srgy_offset += dx;
	    val += global_energy(0,0,chain,NULL, biasmap,mod_params);
	    mod_params->srgy_offset -= 2. * dx;
	    val -= global_energy(0,0,chain,NULL, biasmap,mod_params);
	    mod_params->srgy_offset += dx;
	    this[20] = 0.5 * val / dx;
	}

	/* atomic radii of LJ vdW potential */
	if (sim_params->energy_probe_1_calc[21])
	    this[21] = global_dfdx(chain,biasmap,all_vdw,&(mod_params->rca),0.01, mod_params);
	if (sim_params->energy_probe_1_calc[22])
	    this[22] = global_dfdx(chain,biasmap,all_vdw,&(mod_params->rcb),0.01, mod_params);
	if (sim_params->energy_probe_1_calc[23])
	    this[23] = global_dfdx(chain,biasmap,all_vdw,&(mod_params->rc),0.01, mod_params);
	if (sim_params->energy_probe_1_calc[24])
	    this[24] = global_dfdx(chain,biasmap,all_vdw,&(mod_params->rn),0.01, mod_params);
	if (sim_params->energy_probe_1_calc[25])
	    this[25] = global_dfdx(chain,biasmap,all_vdw,&(mod_params->ro),0.01, mod_params);
	if (sim_params->energy_probe_1_calc[26])
	    this[26] = global_dfdx(chain,biasmap,all_vdw,&(mod_params->rs),0.01, mod_params);

	/* stress */
	if (sim_params->energy_probe_1_calc[27])
	    this[27] = intraresidual_dfdx(chain,biasmap,stress,&(mod_params->stress_k),1.0, mod_params);

	/* depth of LJ vdW potential */
	if (sim_params->energy_probe_1_calc[28])
	    this[28] = global_dfdx(chain,biasmap,all_vdw,&(mod_params->vdw_depth_ca),0.001, mod_params);
	if (mod_params->vdw_uniform_depth) {
	    for (int i=29; i<=33; i++) this[i] = 0;
	} else {
	    if (sim_params->energy_probe_1_calc[29])
		this[29] = global_dfdx(chain,biasmap,all_vdw,&(mod_params->vdw_depth_cb),0.001, mod_params);
	    if (sim_params->energy_probe_1_calc[30])
		this[30] = global_dfdx(chain,biasmap,all_vdw,&(mod_params->vdw_depth_c),0.001, mod_params);
	    if (sim_params->energy_probe_1_calc[31])
		this[31] = global_dfdx(chain,biasmap,all_vdw,&(mod_params->vdw_depth_n),0.001, mod_params);
	    if (sim_params->energy_probe_1_calc[32])
		this[32] = global_dfdx(chain,biasmap,all_vdw,&(mod_params->vdw_depth_o),0.001, mod_params);
	    if (sim_params->energy_probe_1_calc[33])
		this[33] = global_dfdx(chain,biasmap,all_vdw,&(mod_params->vdw_depth_s),0.001, mod_params);
	}

	/* secondary radius of gyration global energy */
	if (sim_params->energy_probe_1_calc[34]) {
	    val = 0.;
	    dx = 0.09;
	    mod_params->hphobic_srgy_param += dx;
	    val += global_energy(0,0,chain,NULL,biasmap, mod_params);
	    mod_params->hphobic_srgy_param -= 2. * dx;
	    val -= global_energy(0,0,chain,NULL,biasmap, mod_params);
	    mod_params->hphobic_srgy_param += dx;

	    this[34] = 0.5 * val / dx;
	}

	if (sim_params->energy_probe_1_calc[35]) {
	    val = 0;
	    dx = 0.09;
	    mod_params->hphobic_srgy_offset += dx;
	    val += global_energy(0,0,chain,NULL, biasmap,mod_params);
	    mod_params->hphobic_srgy_offset -= 2. * dx;
	    val -= global_energy(0,0,chain,NULL, biasmap,mod_params);
	    mod_params->hphobic_srgy_offset += dx;
	    this[35] = 0.5 * val / dx;
	}

	// If the hydrophobic potential form is 1/dist, more parameters need optimising)
	//if (sim_params->energy_probe_1_calc[34])
	//    this[34] = dfdx(chain,biasmap,hydrophobic, &(mod_params->hydrophobic_r), 0.01, mod_params);
	//if (sim_params->energy_probe_1_calc[35])
	//    this[35] = dfdx(chain,biasmap,hydrophobic, &(mod_params->hydrophobic_half_delta), 0.01, mod_params);


    out:

	/* save these results into this and calculate the gradient */
	for (int i=0; i<36; i++) {
	    sim_params->energy_probe_1_this[i] = this[i];
	    sim_params->energy_gradient[i] = sim_params->energy_probe_1_this[i] - sim_params->energy_probe_1_last[i];
	//    fprintf(sim_params->outfile,"%f ", sim_params->energy_gradient[i]);
	}
	//putchar('\n');

	model_param_finalise(mod_params);
	free(mod_params);

}

/****                ENERGY CONTRIBUTIONS               ****/

void exclude_energy_contributions_in_energy_c(Chain * chain,Biasmap *biasmap, double tote, model_params *mod_params, FILE *outfile)
{
	
	int i,j;
	double exclude_energy = 0;
	double d2;

	//FILE *fptr = fopen("g.hydro","a");
	
	for(i = 1; i < chain->NAA; i++) {
	  for(j = i+1; j < chain->NAA; j++){
		switch (j - i) {
		case 1:
			exclude_energy = exclude_neighbor(&(chain->aa[i]),&(chain->aa[j]), mod_params);
			break;
		default:
			d2 = distance(chain->aa[i].ca, chain->aa[j].ca);
			if (d2 < mod_params->vdw_extended_cutoff) {
				exclude_energy = exclude(&(chain->aa[i]), &(chain->aa[j]), d2, mod_params);
/*				if(exclude(&aa[i],&aa[j],d2) > 5){ 
				  fprintf(outfile,"%d %d %c %c %f %f",i,j,aa[i].id,aa[j].id,d2, exclude(&aa[i],&aa[j],d2));
                  fprintf(outfile," %f\n", sqrt(distance(*&aa[i].g,*(&aa[j].g))));  
			} 
*/
			}
			break;
		}
		fprintf(outfile,"%d %d %g\n",i,j,exclude_energy);
	  }
	fprintf(outfile,"\n");
    } 
	fprintf(outfile,"\n");
}


void energy_contributions_in_energy_c(Chain * chain,Biasmap *biasmap, double tote, model_params *mod_params, FILE *outfile) {
	
	int i,j;
	double eenergy = 0, henergy = 0;
	double stress_energy = 0, clash_energy = 0;
	double bias_energy = 0, hbond_energy = 0, exclude_energy = 0, proline_energy = 0;
	double exclude_neighbour_energy = 0;
	double sidechain_hbond_energy = 0;
			double d2;

	//FILE *fptr = fopen("g.hydro","a");
	
	for(i = 1; i < chain->NAA; i++) {
	  stress_energy += stress(&(chain->aa[i]), mod_params);
	  clash_energy += clash(&(chain->aa[i]), mod_params); 
	  //vdwrad(&aa[i],&aa[i]);
	  for(j = i+1; j < chain->NAA; j++){
	    //vdwrad(&aa[i],&aa[j]);	    
	    eenergy += electrostatic(biasmap,&(chain->aa[i]),&(chain->aa[j]), mod_params);
	    double temp = hydrophobic(biasmap,&(chain->aa[i]),&(chain->aa[j]), mod_params);
	    //fprintf(fptr,"%f ",temp);
	    henergy += temp;
	    bias_energy += bias(biasmap,&(chain->aa[i]),&(chain->aa[j]), mod_params);

	    sidechain_hbond_energy += sidechain_hbond(biasmap,&(chain->aa[i]),&(chain->aa[j]), mod_params);
	    int seqdist;
	    if (chain->aa[i].chainid == chain->aa[j].chainid )
		seqdist = chain->aa[i].num - chain->aa[j].num;
	    else
		seqdist = 1000 * abs( chain->aa[j].chainid - chain->aa[i].chainid);

	    switch ( seqdist ) {
		//switch (j - i) {
		case 1:
			exclude_energy += exclude_neighbor(&(chain->aa[i]),&(chain->aa[j]), mod_params);
			exclude_neighbour_energy += exclude_neighbor(&(chain->aa[i]),&(chain->aa[j]), mod_params);
			proline_energy += proline(&(chain->aa[i]),&(chain->aa[j]));
			hbond_energy += hbond(biasmap,&(chain->aa[i]), &(chain->aa[j]), mod_params);
			break;
		default:
			d2 = distance(chain->aa[i].ca, chain->aa[j].ca);
			if (d2 < mod_params->vdw_extended_cutoff) {
				exclude_energy += exclude(&(chain->aa[i]), &(chain->aa[j]), d2, mod_params);
/*				if(exclude(&aa[i],&aa[j],d2) > 5){ 
				  fprintf(outfile,"%d %d %c %c %f %f",i,j,aa[i].id,aa[j].id,d2, exclude(&aa[i],&aa[j],d2));
                  fprintf(outfile," %f\n", sqrt(distance(*&aa[i].g,*(&aa[j].g))));  
			} 
*/
				if (d2 < hbond_cutoff) {
					hbond_energy += hbond(biasmap,&(chain->aa[i]), &(chain->aa[j]), mod_params);
				}
			}
			break;
		}
	  }
    } 
	//fprintf(fptr,"\n"); fclose(fptr); //free(fptr);
	double global = global_energy(0,0,chain,NULL, biasmap,mod_params);
    tote += global;
	fprintf(outfile,"#%20s || %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s\n","total","stress","clash","bias","hbond","exclude","proline","hydrophobic","electrostatic","sidechain hbond","global","exclude_neighbour");
	fprintf(outfile,"%20.6f || %20.6f %20.6f %20.6f %20.6f %20.6f %20.6f %20.6f %20.6f %20.6f %20.6f %20.6f\n",tote,stress_energy,clash_energy,bias_energy,hbond_energy,exclude_energy,proline_energy,henergy,eenergy,sidechain_hbond_energy,global,exclude_neighbour_energy);
	
}
