/*
** Metropolis Monte Carlo sampling procedure for simplified polypeptides.
**
** Copyright (c) 2004 - 2008 Alexei Podtelezhnikov
** Copyright (c) 2007 - 2013 Nikolas Burkoff, Csilla Varnai and David Wild
*/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<time.h>
#include<signal.h>
#include<math.h>
#include <dirent.h>
#include <errno.h>

#ifdef PARALLEL
#include<mpi.h>
#include"random16.h"
#endif

//#include"canonicalAA.h"
#include"rotamers.h"
#include"error.h"
#include"params.h"
#include"aadict.h"
#include"vector.h"
#include"rotation.h"
#include"peptide.h"
#include"vdw.h"
#include"energy.h"
#include"metropolis.h"
#include"probe.h"
#include"nested.h"
#include"checkpoint_io.h"
#include"flex.h"

#define VER "ADCP 0.1, Copyright (c) Yuqi Zhang, Michel Sanner, CCSB Scripps \n\
2004 - 2010 Alexei Podtelezhnikov\n\
Nested Sampling by N.Burkoff\n"
#define USE "Usage: %s [options] [[-f] infile | SEQuenCE] [-o outfile]\n\
Options:\n\
 -f infile            input PDB file with initial conformation\n\
 or SEQuenCE          peptide sequence in ALPHA and beta states\n\
 -o outfile           redirected output file\n\
 -L rotamerLibs       ':' separated list of rotamer files from ADCP/data/rotamers e.g. 'fluo:swiss'\n\
 -l userRotamerLibs   ':' separated list of user-specide rotamer files eg. -l ./myRots.lib'\n\
 -T targetFolder      folder providing the .map and transpoint files\n\
 -a ACCEPTANCE        crankshaft rotation acceptance rate\n\
 -A AMPLITUDE,FIX_AMP crankshaft rotation amplitude and whether the amplitude should be kept fixed (default is no fixing) \n\
 -b BETA1-BETA2:INT   thermodynamic beta schedule\n\
 -p STRING            undocumented custom parameters\n\
 -d vdw_model         Which vdW model to set the defaults for (default: LJ, could also be hard_cutoff and LJ_hard_cutoff)\n\
 -r PACExSTRETCH      test interval x total number\n\
 -s SEED              random seed\n\
 -t MASK,OPTIONS      hexadecimal mask of active tests\n\
 -c TEMP			  temperature (Celcius) to run serial MC simulation\n\
 \n\
 -n                   nested sampling procedure\n\
 -c	TEMP (NS)	  	  minimum temperature (Celcius)\n\
 -m NUM               number of MC moves between NS sample points\n\
 -r OUTPUTxMAX (NS)   output every OUTPUT sample point x max number of NS iterations\n\
 -C NUM,FILENAME      checkpointing: number NS iterations until checkpoint,checkpointfilename\n\
 -R N                 restart checkpointing - need -C, NUM=checkpointfile number\n"


#define M_LOG2E        1.4426950408889634074   


#ifdef PARALLEL
/* parallel job parameters */
int size = 1, rank = 0;

/* set thermodynamic beta for parallel tempering replicas */
void thermoset(simulation_params *sim_params)
{
	int i;
	double beta = 1.0, factor;

	if (size > 1 && sim_params->beta2 > 0.0 && sim_params->beta1 > 0.0)
		factor = pow(sim_params->beta2 / sim_params->beta1, 1.0 / (size - 1));
	else
		factor = 1.0 - 0.5 * M_LOG2E * log(size + 1) / size;

	if (sim_params->beta1 > 0.0)
		beta = sim_params->beta1;

	for (i = 0; i < rank; i++)
		beta *= factor;

	sim_params->thermobeta = beta;
	sim_params->bstp = factor;

    

	fprintf(stderr, "Rank = %d  ThermoBeta = %g\n", rank, sim_params->thermobeta);
}

void thermoswap(Chain* chain, simulation_params *sim_params)
{
	int i, j;
	static int swap = 0;
	double loss, p;
	double send[2], recv[2];
	MPI_Status status;

	swap++;

	/* who's up for a swap and what's the probability? */
	j = rand16();
	i = (j >> 8) % size;
	j = j % size;
	p = rand16();

	if ((rank != i && rank != j) || i == j)
		return;

	send[0] = sim_params->thermobeta;
	send[1] = totenergy(chain);

	if (rank == i) {
		MPI_Recv(recv, 2, MPI_DOUBLE, j, j * size + i, MPI_COMM_WORLD,
			 &status);
		MPI_Send(send, 2, MPI_DOUBLE, j, i * size + j, MPI_COMM_WORLD);
	}

	if (rank == j) {
		MPI_Send(send, 2, MPI_DOUBLE, i, j * size + i, MPI_COMM_WORLD);
		MPI_Recv(recv, 2, MPI_DOUBLE, i, i * size + j, MPI_COMM_WORLD,
			 &status);
	}

	loss = (recv[0] - send[0]) * (recv[1] - send[1]);
	if (loss < 0.0 && exp(loss) * RAND16_MAX < p)
		return;

	sim_params->thermobeta = recv[0];	/* swap */

	fprintf(stderr, "swap %d : rank %d : %g x %g <=> %g x %g\n",
		swap, rank, send[0], send[1], recv[0], recv[1]);
}

#else
void thermoset(simulation_params *sim_params)
{
	if (sim_params->beta2 > 0.0 && sim_params->beta1 > 0.0 && sim_params->pace > 0 && sim_params->stretch > 0 && sim_params->intrvl > 0)
		sim_params->bstp = pow(sim_params->beta2 / sim_params->beta1, sim_params->intrvl / ((double)sim_params->pace * sim_params->stretch));
	else
		sim_params->bstp = 1.0;

	if (sim_params->beta1 > 0.0 && sim_params->lowtemp == 0)
		sim_params->thermobeta = sim_params->beta1;
		
}

void thermoswap(Chain *chain,simulation_params *sim_params)
{
	static int swap = 0;
	//fprintf(stderr, "annealing %d : %g \n ", swap, sim_params->thermobeta);
	if (sim_params->bstp == 1.0)
		return;

	swap++;
	sim_params->thermobeta *= sim_params->bstp;

	fprintf(stderr, "annealing %d : %g\n", swap, sim_params->thermobeta);
}
#endif

static double calculateRMSD(Chain *chain, Chain *chain2)
{
	//return 3.5;
	double RMSD = 0;
	double dist = 0;
	for (int i = 1; i < chain->NAA; i++) {
		dist = (chain->aa[i].ca[0] - chain2->aa[i].ca[0])*(chain->aa[i].ca[0] - chain2->aa[i].ca[0]);
		dist += (chain->aa[i].ca[1] - chain2->aa[i].ca[1])*(chain->aa[i].ca[1] - chain2->aa[i].ca[1]);
		dist += (chain->aa[i].ca[2] - chain2->aa[i].ca[2])*(chain->aa[i].ca[2] - chain2->aa[i].ca[2]);
		RMSD += dist;
	}
	return (RMSD / (chain->NAA-1));	
}

void simulate(Chain * chain, Chaint *chaint, Biasmap* biasmap, simulation_params *sim_params)
{
	unsigned int i, j, k = sim_params->intrvl;
	double temp;
	Chain *chain2 = (Chain *)malloc(sizeof(Chain));
	chain2->aa = NULL; chain2->xaa = NULL; chain2->erg = NULL; chain2->xaa_prev = NULL;
	allocmem_chain(chain2,chain->NAA,chain->Nchains);

	energy_matrix_print(chain, biasmap, &(sim_params->protein_model));
	//stop("I will stop here,\n");
	if (sim_params->protein_model.opt == 1) {
		// targetBest and currTargetEnergy are two global variables
		targetBest = 99999.;
		currTargetEnergy = 99999.;
		//double targetBestTemp = targetBest;
		//double targetBestPrev = targetBest;
		double lastTargetEnergy = 9999.;		
		int lastIndex = 0; //Index for last last good energy 
		int resetIndex = 0; //Index for last annealing reset
		int lastGoodIndex = 0; //Index for last good energy swap in
		int bestIndex = 0; //Index for last best energy found
		int mutateIndex = -999999; //Index for last transmutate
		int moved = 0;
		int stopSignal = 0;
		int swapInd = 0;
		int swapInd2 = 0;
		int swapLength = 10; //size of swapping pool
		int inCache = 0;
		int ind = 0;
		int currIndex = 0;
		int stuckcount = 0;
		//FILE *swapFile = NULL;
		char swapname[12];
		sprintf(swapname, "swap%d.pdb", swapLength);
		double swapEnergy[swapLength + 1];
		Chain* swapChains[swapLength + 1];

		/*optimizing parameters*/
		int noImprovHeatSteps = 1000000;
		int noImprovStopSteps = 30000000;
		int swapBadSteps = 50000;
		int swapMutateSteps = 200000;
		int swapGoodSteps = 100000;
		double goodEnergyDiff = 5; //5kcal=8.33 3kcal=5
		double rmsdCutoff = 4; // swapping clusters rmsd cutoff
		double heatFactor = 0.5; // starting temp while annealing
		double annealFactor = 1.02; // 2 = 1.02^35 = 1.015^47 = 1.012^58 = 2 ,first 35 *10000 to reach room temperature
		int annealSteps = 10000; //temp change step in annealing
		int swapAneal = 1; //swapping while annealing, 0 false, 1 true.
		int swapGoodProb = 5000; //chance of swapping a good pose after swapGoodSteps, out of 10000
		int swapBadProb = 9500; //chance of swapping a bad pose after swapBadSteps, out of 10000
		int initTransMutate = 0; //initialize with randomize translation
		int swapTransMutate = 1; //transMutate during simulation
		int swapFlipChain = 0; //flipChain during simulation

		if (sim_params->infile == NULL) initTransMutate = 1; //default for seq input initialize with randomize translation

		//initialize swapping pool, last element is with the best energy
		for (int i = 0; i < swapLength + 1; i++) {
			swapChains[i] = (Chain *)malloc(sizeof(Chain));
			swapChains[i]->aa = NULL; swapChains[i]->xaa = NULL; swapChains[i]->erg = NULL; swapChains[i]->xaa_prev = NULL;
			allocmem_chain(swapChains[i], chain->NAA, chain->Nchains);
			copybetween(swapChains[i], chain);
			swapEnergy[i] = 9999.;
		}
		//hack, external2_k overwrite external_k
		if (sim_params->protein_model.external_potential_type2 == 4) sim_params->protein_model.external_k[0] = sim_params->protein_model.external_k2[0];
		double external_k = sim_params->protein_model.external_k[0];
		sim_params->protein_model.external_k[0] = external_k * heatFactor;
		if (initTransMutate == 1) {
			fprintf(stderr, "initial transmutate %s \n",sim_params->sequence);
			transmutate(chain, chaint, biasmap, 0, 0, &temp, sim_params);
		}

		fprintf(stderr, "begin run with pace %d stretch %d \n", sim_params->pace, sim_params->stretch);
		for (i = 1; i < sim_params->stretch * sim_params->pace; i++) {
			if (stopSignal) break;
			//targetBestTemp = targetBest;
			if (!sim_params->keep_amplitude_fixed) { // potentially alter amplitude
				if ((i % 1000000 == 1 && i < 10000000) || (i % 10000000 == 1)) {
					energy_matrix_print(swapChains[swapLength], biasmap, &(sim_params->protein_model));
					/* This bit ensures the amplitude of the moves
					is independent of the chain's history*/
					//copybetween(chain2, chain);
					move(chain, chaint, biasmap, 0.0, &temp, -1, sim_params);
					//for (j = 1; (j < sim_params->pace || j < 1024); j++) {
					//	move(chain2, chaint, biasmap, 0.0, &temp, 1, sim_params);
					//}
				}else{
					moved = move(chain, chaint, biasmap, 0, &temp, 0, sim_params);
				}
			}

			if (moved && rand()%1000 < 10)
				transopt(chain, chaint, biasmap, 0, 0, &temp, sim_params, 0);

			inCache = 0;
			currIndex = i;
			currTargetEnergy = sim_params->protein_model.opt_totE_weight*totenergy(chain)
				+ sim_params->protein_model.opt_extE_weight*extenergy(chain)
				+ sim_params->protein_model.opt_firstlastE_weight*locenergy(chain);

			//do a hard minimization if energy is good
			if (currTargetEnergy - targetBest <= goodEnergyDiff) {
				// also do a hard minimization
				transopt(chain, chaint, biasmap, 0, 0, &temp, sim_params, 1);
				// update the energy
				currTargetEnergy = sim_params->protein_model.opt_totE_weight*totenergy(chain)
								+ sim_params->protein_model.opt_extE_weight*extenergy(chain)
								+ sim_params->protein_model.opt_firstlastE_weight*locenergy(chain);
			}

			if (sim_params->protein_model.external_k[0] < external_k && currIndex % annealSteps == 0) {
				sim_params->protein_model.external_k[0] = annealFactor*sim_params->protein_model.external_k[0];
				//fprintf(stderr, "annealing %g \n", sim_params->protein_model.external_k[0]);
			}
			else if (sim_params->protein_model.external_k[0] > external_k) {
				fprintf(stderr, "annealing complete \n");
				sim_params->protein_model.external_k[0] = external_k;
			}

			if (currTargetEnergy - lastTargetEnergy < 0.001 && currTargetEnergy - lastTargetEnergy > -0.001) {
				stuckcount++;
				if (stuckcount >= 1000) {
					swapInd = rand() % (swapLength + 1);
					while (swapEnergy[swapInd] >= currTargetEnergy) swapInd = rand() % (swapLength + 1);
					//swapInd = swapLength;
					fprintf(stderr, "swap out stuck curr %g swap %g best %g\n", currTargetEnergy, swapEnergy[swapInd], targetBest);
					copybetween(chain, swapChains[swapInd]);
					lastTargetEnergy = swapEnergy[swapInd];
					lastIndex = currIndex;
					stuckcount = 0;
				}
				continue;
			}
			else {
				stuckcount = 0;
				lastTargetEnergy = currTargetEnergy;

			}
			//best energy found
			if (currTargetEnergy - targetBest < -0.001) {
				//reset temp;
				if (sim_params->protein_model.external_k[0] != external_k && currIndex > 100000) {
					fprintf(stderr, "best energy found, reset temp\n");
					sim_params->protein_model.external_k[0] = external_k;
				}

				//record energy and reset indices;
				resetIndex = currIndex;
				bestIndex = currIndex;
				lastIndex = currIndex;
				lastGoodIndex = currIndex;
				//write to swap pool
				for (ind = 0; ind < swapLength; ind++) {
					if (swapEnergy[ind] == 9999.) swapInd = ind;
					if (calculateRMSD(swapChains[ind], chain) < rmsdCutoff) {
						inCache = 1;
						//fprintf(stderr, "inCache");
						if (currTargetEnergy < swapEnergy[ind]) {
							//fprintf(stderr, "inCache %d \n",ind);
							targetBest = currTargetEnergy;
							fprintf(stderr, "swap between best curr %g swap %g best %g\n", currTargetEnergy, swapEnergy[ind], targetBest);
							copybetween(swapChains[ind], chain);
							swapEnergy[ind] = currTargetEnergy;
							//sprintf(swapname, "swap%d.pdb", ind);
							//swapFile = fopen(swapname, "w+");
							//pdbprint(chain->aa, chain->NAA, &(sim_params->protein_model), swapFile, &currTargetEnergy);
						}
						break;
					}
				}
				if (inCache != 1) {
					if (swapEnergy[swapInd] != 9999.) {
						swapInd = rand() % swapLength;
						swapInd2 = rand() % swapLength;
						swapInd = swapEnergy[swapInd] > swapEnergy[swapInd2] ? swapInd : swapInd2;
					}
					
					//sprintf(swapname, "swap%d.pdb", swapInd);
					//swapFile = fopen(swapname, "w+");
					//pdbprint(chain->aa, chain->NAA, &(sim_params->protein_model), swapFile, &currTargetEnergy);
					//fclose(swapname);
					//copybetween(chain, swapChains[swapInd]);	
					if (swapEnergy[swapInd] < 0) {
						fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
						tests(swapChains[swapInd], biasmap, sim_params->tmask, sim_params, 0x11, NULL);
					}

					fprintf(stderr, "swap in best curr %g swap %g best %g\n", currTargetEnergy, swapEnergy[swapInd], targetBest);
					copybetween(swapChains[swapInd], swapChains[swapLength]);
					swapEnergy[swapInd] = swapEnergy[swapLength];
					//fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
					//tests(chain, biasmap, sim_params->tmask, sim_params, 0x11, NULL);
					//copybetween(swapChains[swapInd2], swapChains[swapLength]);
					//swapEnergy[swapInd2] = swapEnergy[swapLength];
				}
				targetBest = currTargetEnergy;
				//write out best solutions to output pdb
				if (currTargetEnergy < 0) {
					fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
					tests(chain, biasmap, sim_params->tmask, sim_params, 0x11, NULL);
				}
				//fprintf(stderr, "best found %g \n", swapEnergy[swapLength]);
				//write to last element of swap pool;
				copybetween(swapChains[swapLength], chain);
				swapEnergy[swapLength] = currTargetEnergy;
				//sprintf(swapname, "swap%d.pdb", swapLength);
				//swapFile = fopen(swapname, "w+");
				//pdbprint(chain->aa, chain->NAA, &(sim_params->protein_model), swapFile, &currTargetEnergy);
				//fclose(swapname);
				//fprintf(stderr, "best found2 %g \n", swapEnergy[swapLength]);
				continue;
			}
			else if ((currIndex - bestIndex) > noImprovStopSteps) {
				fprintf(stderr, "No improvement after %d runs last best %d, stops here.\n", i, bestIndex);
				stopSignal = 1;
				break;
			}
			else if ((currIndex - resetIndex) > noImprovHeatSteps && currIndex - mutateIndex > swapMutateSteps && sim_params->protein_model.external_k[0] == external_k) {
				sim_params->protein_model.external_k[0] = external_k * heatFactor;
				fprintf(stderr, "No improvement after %d, Heat up system\n", noImprovHeatSteps);
				swapInd = rand() % (swapLength + 1);
				while (swapEnergy[swapInd] > targetBest + goodEnergyDiff) swapInd = rand() % (swapLength + 1);
				//fprintf(stderr, "swap in best %d \n", swapInd);
				//swapInd = swapLength;

				copybetween(chain, swapChains[swapInd]);



				if (swapTransMutate == 1 && rand()%10000 > swapBadProb){
					fprintf(stderr, "heat and transmutate curr %g best %g curriter %d \n", currTargetEnergy, targetBest, i);
					transmutate(chain, chaint, biasmap, 0, 0, &temp, sim_params);
					transopt(chain, chaint, biasmap, 0, 0, &temp, sim_params, 1);
					mutateIndex = currIndex;
				} else if (rand()%100 < 0 && targetBest < 0 && swapFlipChain) {
					fprintf(stderr, "before flip %g \n",extenergy(chain));
					rotate_cyclic(chain, chaint, biasmap, 0, 0, &temp, sim_params);
					transopt(chain, chaint, biasmap, 0, 0, &temp, sim_params, 1);
					fprintf(stderr, "after flip %g \n",extenergy(chain));
					mutateIndex = currIndex;
				}
				lastGoodIndex = currIndex;
				resetIndex = currIndex;
				lastIndex = currIndex;
				fprintf(stderr, "swap out no improv curr %g swap %g best %g \n", currTargetEnergy, swapEnergy[swapInd], targetBest);
			}
			// good energy found
			else if (currTargetEnergy - targetBest <= goodEnergyDiff) {
				// check RMSD with the swapping pool
				for (ind = 0; ind < swapLength; ind++) {
					if (swapEnergy[ind] == 9999.) swapInd = ind;
					// it is within the clusters
					if (calculateRMSD(swapChains[ind], chain) < rmsdCutoff) {
						inCache = 1;
						// update the energy and swap in if curr has better energy
						if (currTargetEnergy < swapEnergy[ind]) {
							copybetween(swapChains[ind], chain);
							fprintf(stderr, "swap between good curr %g swap %g best %g\n", currTargetEnergy, swapEnergy[ind], targetBest);
							swapEnergy[ind] = currTargetEnergy;
							//fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
							//tests(chain, biasmap, sim_params->tmask, sim_params, 0x11, NULL);
							//sprintf(swapname, "swap%d.pdb", ind);
							//swapFile = fopen(swapname, "w+");
							//pdbprint(chain->aa, chain->NAA, &(sim_params->protein_model), swapFile, &currTargetEnergy);
							lastGoodIndex = currIndex;
						}
						break;
					}
				}
				// it is a new cluster, swap in
				if (inCache != 1) {
					if (swapEnergy[swapInd] != 9999.) {
						swapInd = rand() % swapLength;
						swapInd2 = rand() % swapLength;
						swapInd = swapEnergy[swapInd] > swapEnergy[swapInd2] ? swapInd : swapInd2;
					}
					if (swapEnergy[swapInd]>currTargetEnergy){
						//copybetween(chain, swapChains[swapInd]);
						fprintf(stderr, "swap in good curr %g swap %g best %g\n", currTargetEnergy, swapEnergy[swapInd], targetBest);
						//sprintf(swapname, "swap%d.pdb", swapInd);
						//swapFile = fopen(swapname, "w+");
						//pdbprint(chain->aa, chain->NAA, &(sim_params->protein_model), swapFile, &currTargetEnergy);
						//fclose(swapname);
						if (swapEnergy[swapInd] < 0) {
							fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
							tests(swapChains[swapInd], biasmap, sim_params->tmask, sim_params, 0x11, NULL);
						}
						copybetween(swapChains[swapInd], chain);
						swapEnergy[swapInd] = currTargetEnergy;
						//fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
						//tests(chain, biasmap, sim_params->tmask, sim_params, 0x11, NULL);
						lastGoodIndex = currIndex;
					}
				}

				if (currIndex - lastGoodIndex > swapGoodSteps && currIndex - mutateIndex > swapMutateSteps && rand()%10000<swapGoodProb) {
					if (swapAneal || external_k == sim_params->protein_model.external_k[0]) {

						swapInd = rand() % (swapLength + 1);
						while (swapEnergy[swapInd] > currTargetEnergy) swapInd = rand() % (swapLength + 1);
						//swapInd = swapLength;
						fprintf(stderr, "swap out good curr %g swap %g best %g\n", currTargetEnergy, swapEnergy[swapInd], targetBest);
						copybetween(chain, swapChains[swapInd]);
						if (rand()%100 < 10 && targetBest < 0 && swapFlipChain) {
							fprintf(stderr, "before flip %g \n",extenergy(chain));
							rotate_cyclic(chain, chaint, biasmap, 0, 0, &temp, sim_params);
							transopt(chain, chaint, biasmap, 0, 0, &temp, sim_params, 1);
							fprintf(stderr, "after flip %g \n",extenergy(chain));
							mutateIndex = currIndex;
						}
						lastGoodIndex = currIndex;
						//sim_params->protein_model.external_k[0] = external_k;
					}
				}
				lastIndex = currIndex;
			}
			else if (currIndex - lastIndex > swapBadSteps && currIndex - mutateIndex > swapMutateSteps) {
				if (rand() % 10000 < swapBadProb) {
					if (swapAneal || external_k == sim_params->protein_model.external_k[0]) {
						swapInd = rand() % (swapLength + 1);
						//energy_matrix_print(chain, biasmap, &(sim_params->protein_model));
						while (swapEnergy[swapInd] - targetBest > goodEnergyDiff) swapInd = rand() % (swapLength + 1);
						//swapInd = swapLength;
						fprintf(stderr, "swap out bad curr %g swap %g best %g curriter %d \n", currTargetEnergy, swapEnergy[swapInd], targetBest, i);
						copybetween(chain, swapChains[swapInd]);
						if (rand()%100 < 0 && targetBest < 0 && swapFlipChain) {
							fprintf(stderr, "before flip %g \n",extenergy(chain));
							rotate_cyclic(chain, chaint, biasmap, 0, 0, &temp, sim_params);
							transopt(chain, chaint, biasmap, 0, 0, &temp, sim_params, 1);
							fprintf(stderr, "after flip %g \n",extenergy(chain));
							mutateIndex = currIndex;
						}
							
						lastIndex = currIndex;
						//sim_params->protein_model.external_k[0] = external_k;
					}
				}
				else if (swapTransMutate == 1) {
					swapInd = rand() % (swapLength + 1);
					//energy_matrix_print(chain, biasmap, &(sim_params->protein_model));
					while (swapEnergy[swapInd] > currTargetEnergy) swapInd = rand() % (swapLength + 1);
					//swapInd = swapLength;
					fprintf(stderr, "transmutate bad curr %g best %g curriter %d \n", currTargetEnergy, targetBest, i);
					copybetween(chain, swapChains[swapInd]);
					transmutate(chain, chaint, biasmap, 0, 0, &temp, sim_params);
					transopt(chain, chaint, biasmap, 0, 0, &temp, sim_params, 1);
					//flipChain(chain, chaint, biasmap, 0, 0, &temp, sim_params);

					lastIndex = currIndex;
					mutateIndex = currIndex;
				}
			}

		}
		//clean up
		for (int i = 0; i < swapLength + 1; i++) {
			fprintf(sim_params->outfile, "-+- %5d CLUSTERS BLOCK %5d -+-\n", swapLength+1, i);
			tests(swapChains[i], biasmap, sim_params->tmask, sim_params, 0x11, NULL);
			freemem_chain(swapChains[i]); free(swapChains[i]);
		}

	}
    /* regular MC with bestE recorded */
	else if (sim_params->protein_model.opt == 3) {
		//double targetBestPrev = targetBest;
		double targetBestTemp = targetBest;
		targetBest = 9999.;
		double currTargetEnergy = 99999.;
		//double lastTargetEnergy = 99999.;
		for (i = 1; i < sim_params->stretch; i++) {
			targetBestTemp = targetBest;
			if (!sim_params->keep_amplitude_fixed) { // potentially alter amplitude
				if ((i % 100 == 1 && i < 1000) || (i % 1000 == 1)) {
					/* This bit ensures the amplitude of the moves
					is independent of the chain's history*/
					copybetween(chain2, chain);
					move(chain2, chaint, biasmap, 0.0, &temp, -1, sim_params);
					for (j = 1; (j < sim_params->pace || j < 1024); j++) {
						move(chain2, chaint, biasmap, 0.0, &temp, 1, sim_params);
					}
				}
			}
			targetBest = targetBestTemp;
			for (j = 1; (j < sim_params->pace || j < 1024); j++) {
				//targetBestPrev = targetBest;
				move(chain, chaint, biasmap, 0, &temp, 0, sim_params);
				currTargetEnergy = sim_params->protein_model.opt_totE_weight*totenergy(chain) 
					+ sim_params->protein_model.opt_extE_weight*extenergy(chain)
					+ sim_params->protein_model.opt_firstlastE_weight*locenergy(chain);
				//currTargetEnergy = targetenergy(chain);
				if (currTargetEnergy - targetBest < 0.0) {
					fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
					tests(chain, biasmap, sim_params->tmask, sim_params, 0x11, NULL);
					//lastTargetEnergy = currTargetEnergy;
					targetBest = currTargetEnergy;
					//pdbprint(chain->aa, chain->NAA, &(sim_params->protein_model), "Gary_Hack.pdb", totenergy(chain));
				}
			}
			fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
			tests(chain, biasmap, sim_params->tmask, sim_params, 0x11, NULL);
			//lastTargetEnergy = currTargetEnergy;			
			/* annealing or tempering swap */
			if (--k == 0) {
				thermoswap(chain, sim_params);
				k = sim_params->intrvl;
			}
		}
	}
	/*original MC implementation*/
	else {
		for (i = 1; i <= sim_params->stretch; i++) {
			if (!sim_params->keep_amplitude_fixed) { // potentially alter amplitude
				if ((i % 100 == 1 && i < 1000) || (i % 1000 == 1)) {
					/* This bit ensures the amplitude of the moves
					is independent of the chain's history*/
					copybetween(chain2, chain);
					move(chain2, chaint, biasmap, 0.0, &temp, -1, sim_params);
					for (j = 1; (j < sim_params->pace || j < 1024); j++) {
						move(chain2, chaint, biasmap, 0.0, &temp, 1, sim_params);
					}
				}
			}
			//		for (j = 1; (j < sim_params->pace || j < 1024); j++) {
			for (j = 0; j < sim_params->pace; j++) {
				//sim_params->NS = 1;temp = totenergy(chain);
				move(chain, chaint, biasmap, 0, &temp, 0, sim_params);
				/* annealing or tempering swap */
				if (--k == 0) {
					for (int t = 0; t < sim_params->nswap_per_try; t++) {
						thermoswap(chain, sim_params);
					}
					k = sim_params->intrvl;
				}
			}
#ifdef PARALLEL
			//uncomment this part if you want the energies of each replica
			static FILE* fptr;
			if (fptr == NULL) {
				char filename[1024];
				sprintf(filename, "%d.temp", rank);
				fptr = fopen(filename, "w");
				fprintf(fptr, "#T: %f\n", 1000.0 / (1.9858775*sim_params->thermobeta) - 273.15);
			}
			fprintf(fptr, "%f %f\n", sim_params->thermobeta, totenergy(chain));

			static FILE* fptr_pdb;
			if (fptr_pdb == NULL) {
				char filename_pdb[1024];
				sprintf(filename_pdb, "%d.pdb", rank);
				fptr_pdb = fopen(filename_pdb, "w");
				fprintf(fptr_pdb, "#T: %f\n", 1000.0 / (1.9858775*sim_params->thermobeta) - 273.15);
				double tttt = totenergy(chain);
				pdbprint(chain->aa, chain->NAA, &(sim_params->protein_model), fptr_pdb, &tttt);
			}
			fprintf(fptr, "%f %f\n", sim_params->thermobeta, totenergy(chain));

			if (sim_params->thermobeta != sim_params->beta1)
				continue;
#endif
			fprintf(sim_params->outfile, "-+- TEST BLOCK %5d -+-\n", i);
			tests(chain, biasmap, sim_params->tmask, sim_params, 0x11, NULL);
		}
		
	}
	freemem_chain(chain2); free(chain2);
}

char *getDataFolder(char *argv0) {
  char *path, *p;
#ifdef _WIN32
  int sep='\\';
#else
  int sep='/';
#endif  
  //int sep = argv0[0];
  path = malloc(sizeof(char)*(strlen(argv0)+8));
  p = strrchr(argv0, sep);
  strncpy(path, argv0, p+1-argv0);
  strcat(path, "data");
  int len = strlen(path);
  path[len] = sep;
  path[len+1]= '\0';
  //printf("FOFO %s %s\n", argv0, path);
  return path;
}

int countToken(char *str, char *sep)
{
  char *p;
  int count=0;
  p = strtok(str, sep);
  while (p != NULL) {
    count +=1 ;
    p =strtok(NULL, sep);
  }
  return count;
}

char **getTokens(char *str, char *sep, int nbToken)
{
  char **tokens = malloc(sizeof(char *)*nbToken);
  tokens[0] = strtok(str, sep);
  for (int i=1; i < nbToken; i++) {
    tokens[i] = strtok(NULL, sep);
  }
  return tokens;
}

char *read_options(int argc, char *argv[], simulation_params *sim_params)
{
	unsigned int pace = 0, stretch = 16, tmask = 0x0, seed = 0;
	double thermobeta;
	int lowtemp;
	double beta1 = 1.0, beta2 = 0.0; //, bstp = 1.0;
	unsigned int intrvl = 16384;
	unsigned int nswap_per_try = 16384;
	int NS = 0; //1 = yes
	int iter_max = 1000;
	int num_NS_per_checkpoint = 0;
	char checkpoint_filename[256];
	int checkpoint_counter = 0;
	int checkpoint = 0;
	int restart_from_checkpoint = 0;

	int i, opt;
	char *retval = NULL;
	double acceptance_rate = sim_params->acceptance_rate;
	double amplitude = sim_params->amplitude;
	int keep_amplitude_fixed = sim_params->keep_amplitude_fixed;
	char error_string[DEFAULT_LONG_STRING_LENGTH]="";

	// MS used for rotamer library filenames
	char *copy;
	int osSep;
#ifdef _WIN32
	osSep = '\\';
#else
	osSep = '/';
#endif

	sim_params->data_folder = getDataFolder(argv[0]);
	printf("dataFolder: %s\n", sim_params->data_folder);
	
	for (i = 1; i < argc; i++) {
	  if (argv[i][0] != '-') {
	    //if (freopen(argv[i], "r", stdin) == NULL)
	    //try to open file
	    if ((sim_params->infile = fopen(argv[i],"r")) == NULL)
	      //failed: it must be a sequence
	      retval = argv[i];
	    else
	      if (sim_params->infile_name) free(sim_params->infile_name);
	    copy_string(&(sim_params->infile_name),argv[i]);
	    continue;
	  }

	  opt = argv[i][1];
	  if (++i >= argc && opt != 'n')
	    opt = 0;

	  switch (opt) {
	    case 'a':
	      sscanf(argv[i], "%lf", &acceptance_rate);
	      if(acceptance_rate > 1.0 || acceptance_rate <= 0.0){
		fprintf(stderr,"Acceptance rate must be between 0.0 and 1.0, setting it to 0.5\n");
		acceptance_rate = 0.5; 	
	      }
	      sim_params->acceptance_rate = acceptance_rate;			
			
	      break;
	  case 'A': //amplitude
	    sscanf(argv[i], "%lf,%d", &amplitude,&keep_amplitude_fixed);
	    if(amplitude > 0.0 || amplitude <= -M_PI){
	      fprintf(stderr,"Amplitude must be between -PI and 0.0, setting it to -0.25\n");
	      amplitude = -0.25;
	    }
	    sim_params->amplitude = amplitude;
	    sim_params->keep_amplitude_fixed = keep_amplitude_fixed;
	    break;
	  case 'b':
	    sscanf(argv[i], "%lf-%lf:%u,%u", &beta1, &beta2, &intrvl, &nswap_per_try);
	    sim_params->beta1 = beta1;
	    sim_params->beta2 = beta2;
	    sim_params->intrvl = intrvl;
	    sim_params->nswap_per_try = nswap_per_try;
	    break;
	  case 'c':
	    sscanf(argv[i],"%lf",&thermobeta);
	    thermobeta = 298.0/(thermobeta+273.0);
	    lowtemp = 1;
	    sim_params->thermobeta = thermobeta;
	    sim_params->lowtemp = lowtemp;
	    break;
	  case 'C':
	    sscanf(argv[i], "%d,%256s",&num_NS_per_checkpoint,checkpoint_filename);
	    checkpoint = 1;
	    sim_params->num_NS_per_checkpoint = num_NS_per_checkpoint;
	    sim_params->checkpoint_filename = realloc(sim_params->checkpoint_filename,DEFAULT_SHORT_STRING_LENGTH);
	    strcpy(sim_params->checkpoint_filename,checkpoint_filename);
	    sim_params->checkpoint = checkpoint;
	    break;
	  case 'f':
	    //if (freopen(argv[i], "r", stdin) == NULL)
	    //try to open file
	    if ((sim_params->infile = fopen(argv[i],"r")) == NULL)
	      //failed: it must be a sequence
	      retval = argv[i];
	    else
	      if (sim_params->infile_name) free(sim_params->infile_name);
	    copy_string(&(sim_params->infile_name),argv[i]);
	    break;
	  case 'm':
	    sscanf(argv[i], "%u", &iter_max);
	    sim_params->iter_max = iter_max;
	    break;
	  case 'M':
	    sscanf(argv[i],"%u",&(sim_params->number_initial_MC));
	    break;
	  case 'n':
	    NS = 1;
	    sim_params->NS = NS;
	    i--;
	    break;
	  case 'o':
	    //freopen(argv[i], "w", stdout);
	    if ((sim_params->outfile = fopen(argv[i],"w")) == NULL)
	      stop("Could not open output file.\n");
	    else
	      if (sim_params->outfile_name) free(sim_params->outfile_name);
	    copy_string(&(sim_params->outfile_name),argv[i]);
	    break;
	  case 'p':
	    if (sim_params->prm) free(sim_params->prm);
	    copy_string(&sim_params->prm,argv[i]);
	    break;
	  case 'd':
	    if (strcmp(argv[i],"hard_cutoff")==0) {
	      sim_params->protein_model.vdw_potential=HARD_CUTOFF_VDW_POTENTIAL;
	      set_hard_cutoff_default_params(&(sim_params->protein_model));
	    } else if (strcmp(argv[i],"lj")==0) {
	      sim_params->protein_model.vdw_potential=LJ_VDW_POTENTIAL;
	      set_lj_default_params(&(sim_params->protein_model));
	    } else if (strcmp(argv[i],"lj_hard_cutoff")==0) {
	      sim_params->protein_model.vdw_potential=LJ_VDW_POTENTIAL;
	      set_lj_default_params(&(sim_params->protein_model));
	      sim_params->protein_model.vdw_lj_neighbour_hard=1;
	      sim_params->protein_model.vdw_lj_hbonded_hard=1;
	    } else {
	      sprintf(error_string,"Unknown value for vdW potential model (%s).  It must be one of hard_cutoff and lj.",argv[i]);
	      stop(error_string);
	    }
	    break;
	  case 'r':
	    sscanf(argv[i], "%ux%u", &pace, &stretch);
	    sim_params->pace = pace;
	    sim_params->stretch = stretch;
	    break;
	  case 'L':
	    copy = malloc(sizeof(char)*strlen(argv[i])+1);
	    strcpy(copy, argv[i]);
	    sim_params->nbRotLibs = countToken(copy, ":");
	    sim_params->rotamer_libs = getTokens(argv[i], ":", sim_params->nbRotLibs);
	    for (int ii=0; ii< sim_params->nbRotLibs; ii++)
	      printf("sysRotLib %d/%d \"%s\"\n", ii, sim_params->nbRotLibs, sim_params->rotamer_libs[ii]);
	    break;
	  case 'l':
	    copy = malloc(sizeof(char)*strlen(argv[i])+1);
	    strcpy(copy, argv[i]);
	    sim_params->nbUserRotLibs = countToken(copy, ":");
	    sim_params->userRotamer_libs = getTokens(copy, ":", sim_params->nbUserRotLibs);
	    break;
	  case 'R':
	    sscanf(argv[i], "%d",&checkpoint_counter);
	    if(checkpoint_counter != -1){
	      restart_from_checkpoint = 1;
	      sim_params->checkpoint_counter = checkpoint_counter;
	      sim_params->restart_from_checkpoint = restart_from_checkpoint;
	    }
	    break;
	  case 's':
	    sscanf(argv[i], "%u", &seed);
	    sim_params->seed = seed;
	    break;
	  case 't':
	    sscanf(argv[i], "%x", &tmask);
	    sim_params->tmask = tmask;
	    break;
	  case 'T':
	    if (sim_params->prm) free(sim_params->target_folder);
	    char folderName[254];
	    strcpy(folderName, argv[i]);
	    if (argv[i][strlen(argv[i])-1]!=osSep)
	      sprintf(folderName, "%s%c", argv[i], osSep);
	    else
	      strcpy(folderName, argv[i]);
	    copy_string(&sim_params->target_folder, folderName);
	    break;
	  default:
	    fprintf(stderr, VER USE PARAM_USE, argv[0]);
	    helps();
	    exit(EXIT_FAILURE);
	  }
	}

	if (sim_params->seq) free(sim_params->seq);
	copy_string(&(sim_params->seq),retval);

	// check that the data folder contains the required files
	char msg[254], buffer[254];
	FILE *f;
	DIR* dir;

	dir = opendir(sim_params->data_folder);
	if (dir) {	  /* Directory exists. */
	  closedir(dir);
	} else if (ENOENT == errno) {
	  sprintf(msg, "data folder %s does not exist", sim_params->data_folder);
	  stop(msg);
	} else {
	  sprintf(msg, "couldn't open data folder %s", sim_params->data_folder);
	  stop(msg);
	  /* opendir() failed for some other reason. */
	}
	strcpy(buffer, sim_params->data_folder);
	strcat(buffer, "ramaprob.data");
	f = fopen(buffer, "r");
	if (f == NULL) {
	  sprintf(msg, "ramaprob.data not found in %s", sim_params->data_folder);
	  stop(msg);
	} else fclose(f);

	// count rotamer table entries
	nbCanAA = 20; // from stdaa.lib

	sim_params->rlfullnames = malloc(sizeof(char *)*(sim_params->nbRotLibs+sim_params->nbUserRotLibs));
	// build full names with paths for system rotamer libraries
	int rli;
	for (rli=0; rli < sim_params->nbRotLibs; rli++) {
	  // length of folder and library + 14 to account for "rotamer/" and ".lib" for security take more
	  sim_params->rlfullnames[rli] = malloc(sizeof(char)*(strlen(sim_params->data_folder)+20+strlen(sim_params->rotamer_libs[rli])));
	  sprintf(sim_params->rlfullnames[rli], "%srotamers%c%s.lib", sim_params->data_folder, osSep, sim_params->rotamer_libs[rli]);
	  f = fopen(sim_params->rlfullnames[rli], "r");
	  if (f == NULL) {
	    sprintf(msg, "rotamer library %s not found\n", sim_params->rlfullnames[rli]);
	    stop(msg);
	  } else {
	    fclose(f);
	    int nbrot = countRotamers(sim_params->rlfullnames[rli]);
	    printf("found %d rotamers in %s\n", nbrot, sim_params->rlfullnames[rli]);
	    nbCanAA += nbrot;
	  }
	}
	for (int jj=0; jj < sim_params->nbUserRotLibs; jj++) {
	  sim_params->rlfullnames[rli] = malloc(sizeof(char)*(strlen(sim_params->userRotamer_libs[jj])+1));
	  strcpy(sim_params->rlfullnames[rli], sim_params->userRotamer_libs[jj]);
	  int nbrot = countRotamers(sim_params->rlfullnames[rli]);
	  printf("found %d rotamers in %s\n", nbrot, sim_params->rlfullnames[rli]);
	  nbCanAA += nbrot;
	  rli++;
	}

	if (sim_params->target_folder == NULL) {
	  sprintf(msg, "No target folder specify please use -T option to provide target data");
	  stop(msg);
	}
	dir = opendir(sim_params->target_folder);
	if (dir) {	  /* Directory exists. */
	  closedir(dir);
	} else if (ENOENT == errno) {
	  sprintf(msg, "data folder %s does not exist", sim_params->target_folder);
	  stop(msg);
	} else {
	  sprintf(msg, "couldn't open data folder %s", sim_params->target_folder);
	  stop(msg);
	  /* opendir() failed for some other reason. */
	}

	strcpy(buffer, sim_params->target_folder);
	strcat(buffer, "constrains");
	f = fopen(buffer, "r");
	if (f == NULL) {
	  sprintf(msg, "constrains file not found in %s", sim_params->target_folder);
	  stop(msg);
	} else fclose(f);

	return sim_params->seq;
}

void graceful_exit(int sig)
{
	exit(EXIT_FAILURE);	/* flushes streams too */
}

void single_point_test(Chain *chain, Chaint *chaint, Biasmap *biasmap, simulation_params *sim_params,unsigned int i)
{

	//fprintf(stderr,"SINGLE POINT\n");
	//fprintf(stderr,"XAA\n");
	//for (int i = 1; i < chain->NAA; i++) {
	//	for (int j=0; j<3; j++) {
	//		for (int k=0; k<3; k++) {
	//			fprintf(stderr,"%g ",chain->xaa[i][j][k]);
	//		}
	//	}
	//	fprintf(stderr,"\n");
	//}
	//fprintf(stderr,"END XAA\n");

	/* correct peptide if needed */
	if(sim_params->protein_model.fixit) fixpeptide(chain->aa, chain->NAA, &(sim_params->protein_model)); // PEPTIDE MODIFICATION!!
	chkpeptide(chain->aa, chain->NAA, &(sim_params->protein_model));
	update_sim_params_from_chain(chain,sim_params); // updating NAA and seq

	biasmap_initialise(chain,biasmap,&(sim_params->protein_model));
	energy_matrix_calculate(chain,biasmap,&(sim_params->protein_model));

	/* tests before projecting the peptide onto the CRANKITE model */
	fprintf(sim_params->outfile,"-+- PLAY BLOCK %5d -+-\n", i);
//	fprintf(sim_params->outfile,"-+- before init %5d -+-\n", i);
	tests(chain,biasmap,sim_params->tmask, sim_params, 0x1, NULL );
//	fprintf(sim_params->outfile,"-+- end before init %5d -+-\n", i);


	initialize(chain,chaint,sim_params); // PEPTIDE MODIFICATION!!
	/* tests before projecting the peptide onto the CRANKITE model */
	energy_matrix_calculate(chain,biasmap,&(sim_params->protein_model));
//	fprintf(sim_params->outfile,"-+- after init %5d -+-\n", i);
	tests(chain,biasmap,sim_params->tmask, sim_params, 0x10, NULL );
//	fprintf(sim_params->outfile,"-+- end after init %5d -+-\n", i);

}

void set_random_seed(simulation_params *sim_params) {

	if (sim_params->seed == 0)
		sim_params->seed = (unsigned int) time(NULL);

	/*random seed*/
#ifdef PARALLEL
	srand(sim_params->seed + 257 * rank);	/* asynchronous RNG */
	/* if MC, not NS, broadcast synchronous RNG */
	if(!sim_params->NS){
	MPI_Bcast(&(sim_params->seed), 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	  srand16(sim_params->seed);		/* synchronous RNG */
	}
#else
	srand(sim_params->seed);
#endif

}

void AD_init(Chain *chain, simulation_params *sim_params) {
	if (sim_params->protein_model.external_potential_type == 5) {
		/* int hasCYS = 0; */
		/* int hasAroC = 0; */
		/* int hasNA = 0; */
		/* for(int i = 1; i < chain->NAA; i++){ */
          	/* 	if(chain->aa[i].id == 'C') */
                /* 		hasCYS = 1; */
		/* 	if(chain->aa[i].id == 'F' || chain->aa[i].id == 'Y' || chain->aa[i].id == 'H') */
                /* 		hasAroC = 1; */
		/* 	if(chain->aa[i].id == 'H') */
                /* 		hasNA = 1; */
    		/* } */
		/* transpts_initialise(); */
		/* gridbox_initialise(); */
		/* /\* elements are 0:C, 1:N, 2:O, 3:HD, 4:SA, 5:CA, 6:NA ,7:elec 8:desolv      *\/ */
		/* for (int i = 0; i< 32; i++) { */
		/*   gridmap_initialise("rigidReceptor.%s.map", i); */
		/* } */
		/* gridmap_initialise("rigidReceptor.C.map", 0); */
		/* gridmap_initialise("rigidReceptor.N.map", 1); */
		/* gridmap_initialise("rigidReceptor.OA.map", 2); */
		/* gridmap_initialise("rigidReceptor.HD.map", 3); */
		/* if (hasCYS) */
		/* 	gridmap_initialise("rigidReceptor.SA.map", 4); */
		/* else */
		/* 	gridmap_initialise("rigidReceptor.C.map", 4); */
		/* if (hasAroC) */
		/* 	gridmap_initialise("rigidReceptor.A.map", 5); */
		/* else */
		/* 	gridmap_initialise("rigidReceptor.C.map", 5); */
		/* if (hasNA) */
		/* 	gridmap_initialise("rigidReceptor.NA.map", 6); */
		/* else */
		/* 	gridmap_initialise("rigidReceptor.C.map", 6); */
		/* gridmap_initialise("rigidReceptor.F.map", 7); */
		/* gridmap_initialise("rigidReceptor.e.map", 8); */
		/* gridmap_initialise("rigidReceptor.d.map", 9); */
		/* //printf("transpoints box initialise success %i %g %g %g \n", transPtsCount, Xpts[0], Ypts[transPtsCount - 1], Zpts[transPtsCount - 1]); */
		/* fprintf(stderr, "AD Grid maps initialisation finished \n"); */

		// MS loop over peptide to find which atom types a present
		// so that we can load the proper maps. This replace the code above
		// that was working for a predefined set or amino acids
		// 31 AutoDock atom types

		// force hasType to 1 for C N OA and HD for backbone atoms
		for (int i = 0; i< 4; i++) hasType[i] = 1;
		// hasType for all other atom types is initially 0
		for (int i = 4; i< 32; i++) hasType[i] = 0;
		// loop over sequence to set hasType for side chain atom types
		for(int i = 1; i < chain->NAA; i++){
		  if (chain->aa[i].sideChainTemplateIndex!=-1) {
		    for (int j=0; j<_AASCRotTable[chain->aa[i].sideChainTemplateIndex].nbAtoms;j++){
		      hasType[_AASCRotTable[chain->aa[i].sideChainTemplateIndex].atypes[j]]= 1;
		    }
		    
		    // since HIS also check HIE and HID that has NA (acceptor, type 6 we also have to load this map
		    if (chain->aa[i].id==8) {
		      hasType[6] = 1;
		    }
		  }
		}
		/* printf("ATOMTYPES: "); */
		/* for (int i = 0; i< 32; i++) printf("%d ",hasType[i]); */
		/* printf("\n"); */

		transpts_initialise(sim_params);
		gridbox_initialise(sim_params);

		char mapname[255];
		sprintf(mapname, "%srigidReceptor.e.map", sim_params->target_folder);
		gridmap_initialise(mapname, -1);
		sprintf(mapname, "%srigidReceptor.d.map", sim_params->target_folder);
		gridmap_initialise(mapname, -2);
		/* elements are 0:C, 1:N, 2:O, 3:HD, 4:SA, 5:CA, 6:NA ,7:elec 8:desolv      */
		for (int i = 0; i< 32; i++) {
		  if (hasType[i]==1) {
		    sprintf(mapname, "%srigidReceptor.%s.map", sim_params->target_folder, atypes[i]);
		    printf("loading map %s\n", mapname);
		    gridmap_initialise(mapname, i);
		  }
		}
		fprintf(stderr, "AD Grid maps initialisation finished \n");
	}
}

int main(int argc, char *argv[])
{
	//set a timer
	time_t startTime = time(NULL);

	//char *seq;
	simulation_params sim_params;

	signal(SIGTERM, graceful_exit);

	int osSep;
#ifdef _WIN32
	osSep = '\\';
#else
	osSep = '/';
#endif

#ifdef PARALLEL
	/* initialise MPI */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	
	/* SET SIMULATION PARAMS */
	param_initialise(&sim_params); //set default
	set_lj_default_params(&(sim_params.protein_model)); // set the default parameters for the LJ model
	read_options(argc, argv, &sim_params);

	//param_print(sim_params,sim_params.outfile); //default + read-in

	set_random_seed(&sim_params);

	/* set different default simulation params for NS and MC */
	if(sim_params.NS){
	  if(sim_params.tmask == 0x0) sim_params.tmask = 0x18803;
	  if(sim_params.pace == 0) {sim_params.pace = 1; sim_params.stretch = 1000;}
	} else {
	  if(sim_params.tmask == 0x0) sim_params.tmask = 0x3;
	  sim_params.tmask |= 0x8800;
	  if(sim_params.pace==0) {sim_params.pace = 4096; sim_params.stretch = 16;} 
	}

	/* set temperature for MC */
	if(!sim_params.NS){
	  thermoset(&sim_params);
	}

	/* INITIALISE PROTEIN MODEL */

	//model_param_initialise(&(sim_params.protein_model));
	model_param_read(sim_params.prm,&(sim_params.protein_model),&(sim_params.flex_params));

	ramaprob_initialise(sim_params.data_folder);

	//intialize_AASCRotTable();
	char buffer[254];

	// allocate table of rotamers structures
	_AASCRotTable = malloc(nbCanAA*sizeof(struct _AASCRot));
	if (_AASCRotTable == NULL) {
	  sprintf(buffer, "failed to allocate rotamer table for %d rotamers\n", nbCanAA);
	  stop(buffer);
	} else {
	  printf("allocated space for rotamers for %d amino acids\n", nbCanAA);
	}
	// load rotamers
	// load stdaa always
	sprintf(buffer, "%s%crotamers%cstdaa.lib", sim_params.data_folder, osSep, osSep);
	int lastInd =0;
	lastInd = initialize_AASCRotTable_from_file(buffer, lastInd);
	    
	// load additional rotamer libraries
	for (int rli=0; rli < sim_params.nbRotLibs+sim_params.nbUserRotLibs; rli++)
	  lastInd = initialize_AASCRotTable_from_file(sim_params.rlfullnames[rli], lastInd);

	if (lastInd-1 > nbCanAA) {
	  sprintf(buffer, "allocate memory for %d entries in rotamer table but read %d\n", nbCanAA, lastInd-1);
	  stop(buffer);
	}
	
	initialize_sidechain_properties(&(sim_params.protein_model));
	vdw_cutoff_distances_calculate(&sim_params, stderr, 0);
	peptide_init();
	param_print(sim_params,sim_params.outfile); //read-in    

	/* HERE STARTS THE ACTUAL SIMULATION */

	if(!sim_params.NS){
	    /* allocate memory for the peptide */
	    Chain *chain = (Chain *)malloc(sizeof(Chain)); chain->NAA = 0;
            Chaint* chaint = (Chaint *)malloc(sizeof(Chaint));
      	    chain->NAA =0;
      	    chain->aa = NULL; chain->xaa = NULL; chain->erg = NULL; chain->xaa_prev = NULL;
      	    chaint->aat = NULL; chaint->xaat = NULL; chaint->ergt = NULL; chaint->xaat_prev = NULL;
	    /* allocate memory for the biasmap */
            Biasmap *biasmap = (Biasmap *)malloc(sizeof(Biasmap));
      	    biasmap->distb = NULL;
	    
	    /* read in / generate the peptide */
	    if (sim_params.seq != NULL) {

		//if (sim_params.protein_model.external_constrained_aalist_file) {
		//	stop("constraints unimplemented!");
		//}
		/* build peptide from scratch, do not do tests */
		build_peptide_from_sequence(chain,chaint,sim_params.seq, &sim_params);
		//FILE* fptr1;
		//fptr1 = fopen("afterBuit.pdb", "w");
		//pdbprint(chain->aa, chain->NAA, &(sim_params.protein_model), fptr1, NULL);
		//fclose(fptr1);
		AD_init(chain,&sim_params);
		mark_fixed_aa_from_file(chain,&sim_params);
		mark_constrained_aa_from_file(chain,&sim_params);
		chkpeptide(chain->aa, chain->NAA, &(sim_params.protein_model));
		update_sim_params_from_chain(chain,&sim_params); // updating NAA and seq
		biasmap_initialise(chain,biasmap,&(sim_params.protein_model));		
		energy_matrix_calculate(chain,biasmap,&(sim_params.protein_model));
	    
	    
	   } else { /* read in peptide */

		unsigned int i;
		int readin = 0;

		if (sim_params.stretch == 0) { /* single point testing */

#ifdef PARALLEL
		   if (rank==0) {  // only do the tests and printing on the master node
		   fprintf(stderr,"WARNING! Single point tests are done on the master node only!  A parallel call with stretch=0 might waste resources.\n");
#endif
			fprintf(stderr,"INFO: 0 MC step was asked for.  Only doing tests on the PDB entries.\n");
			/* while reading, do test one by one, to save memory */

			for (i = 1; pdbin(chain,&sim_params,sim_params.infile) != EOF; i++) {
			    mark_fixed_aa_from_file(chain,&sim_params);
			    mark_constrained_aa_from_file(chain,&sim_params);
			    single_point_test(chain,chaint,biasmap,&sim_params,i);
			}
#ifdef PARALLEL
		    }
#endif

		} else { /* reading in for an MC simulation */
#ifdef PARALLEL
			fprintf(stderr, "INFO: Attempting a parallel tempering simulation.\n");
			/* parallel tempering: upto max(rank) initial configs */
			int i;
			int retv;
			for (i = 0; (retv = pdbin(chain, &sim_params, sim_params.infile)) != EOF && i < rank; i++)  readin = 1;
			if (!(retv > 0) && readin == 0)  fprintf(stderr, "WARNING! Rank %d did not read in a PDB\n", rank);
			else readin = 1; //fix rank 0 (readin==0) warning

#else
				 /* serial MC:  only last entry */
			fprintf(stderr, "INFO: Attempting a serial MC simulation.\n");
			unsigned int i;

			for (i = 1; pdbin(chain, &sim_params, sim_params.infile) != EOF; i++);
			if (i>2) fprintf(stderr, "WARNING! First %d entries in the input PDB will be ignored, only using last one for the MC simulation.\n", i - 1);
			if (i == 1) {
				stop("ERROR! EOF while reading in from input PDB file.");
			}
			readin = 1;
#endif

			/* project the peptide onto the CRANKITE model and do the initial tests for serial MC */
			if (readin) {
				if (sim_params.protein_model.fixit) fixpeptide(chain->aa, chain->NAA, &(sim_params.protein_model));
				chkpeptide(chain->aa, chain->NAA, &(sim_params.protein_model));
				mark_fixed_aa_from_file(chain, &sim_params);
				mark_constrained_aa_from_file(chain, &sim_params);
			}
			update_sim_params_from_chain(chain, &sim_params); // updating NAA and seq
															  //fprintf(stderr,"Updated sim_params->seq: %s\n",sim_params.seq); //this also has a starting A which is not part of the polypeptide

#ifndef PARALLEL
			fprintf(sim_params.outfile, "-+- PLAY BLOCK %5d -+-\n", --i);
			/* initialise the biasmap and energy matrix before testing in case some pre-initialisation tests need energies */
			biasmap_initialise(chain, biasmap, &(sim_params.protein_model));
			energy_matrix_calculate(chain, biasmap, &(sim_params.protein_model));
			tests(chain, biasmap, sim_params.tmask, &sim_params, 0x1, NULL);
#endif
			initialize(chain, chaint, &sim_params); // peptide modification
													/* allocate energy matrix and read in biasmap */
			biasmap_initialise(chain, biasmap, &(sim_params.protein_model));
			AD_init(chain,&sim_params);
			energy_matrix_calculate(chain, biasmap, &(sim_params.protein_model));
#ifndef PARALLEL
			/* initial test */
			tests(chain, biasmap, sim_params.tmask, &sim_params, 0x10, NULL);
#endif
		   }
	    }

	/* MC */
	simulate(chain,chaint,biasmap,&sim_params);
	/* print last snapshots for restart */
#ifdef PARALLEL
	FILE* fptr1;
	char filename1[1024];
	sprintf(filename1, "%d_last.pdb", rank);
	fptr1 = fopen(filename1, "w");
	pdbprint(chain->aa, chain->NAA, &(sim_params.protein_model), fptr1, NULL);
	fclose(fptr1);
#endif

	finalize(chain,chaint,biasmap); //free memory allocated in initialize
     
	} else { /* Nested Sampling, parallel or serial */
	  nestedsampling(sim_params.pace,sim_params.stretch,&sim_params);
	}
	
	param_finalise(&sim_params);

#ifdef PARALLEL
    /* finalise MPI */
    MPI_Finalize();
	if (rank == 0);
#endif

	fprintf(stderr, "best target energy %g\n", targetBest);
	// free memory in AutoPK
	if (sim_params.protein_model.external_potential_type == 5) {
		free(Xpts);
		free(Ypts);
		free(Zpts);
		//for (int atype = 0; atype < sizeof(gridmapvalues) / sizeof(gridmapvalues)[0]; atype++)
		// 	free(gridmapvalues[atype]);
		for (int atypeInd=0; atypeInd<MAX_ATOM_TYPES; atypeInd++)
		    if (hasType[atypeInd]==1) free(gridmapvalues[atypeInd]);
		for (int rli=0; rli < sim_params.nbRotLibs+sim_params.nbUserRotLibs; rli++)
		  free(sim_params.rlfullnames[rli]);
		free(sim_params.rlfullnames);
	}
	free(ramaprob);
	free(alaprob);
	free(glyprob);

	//print out the timing
	fprintf(stderr,"The program has successfully finished in %ld seconds. :)  Bye-bye!\n", time(NULL)- startTime);
	return EXIT_SUCCESS;
}
