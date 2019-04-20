/*=============================================================================
topol : molecular topology
Copyright (C) 2002-2018 Franca Fraternali
Copyright (C) 2008-2018 Jens Kleinjung
Copyright (C) 2002 Luigi Cavallo
Copyright (C) 2002 Kuang Lin and Valerie Hindie
Read the COPYING file for license information.
=============================================================================*/

#include "config.h"
#include "topol.h"

#ifdef MPI
#include <mpi.h>
#endif
extern int nodes;
extern int my_rank;

/*___________________________________________________________________________*/
/** print bonded and non-bonded atom pair */
__inline__ static void print_pair(Str *pdb, int t1, int t2)
{
	fprintf(stderr, "pair %d:%d-%d:%d\n",
		t1, pdb->atom[t1].atomNumber,
		t2, pdb->atom[t2].atomNumber);
}

/*___________________________________________________________________________*/
/** print angle */
__inline__ static void print_angle(Str *pdb, int t1, int t2, int t3)
{
	fprintf(stderr, "angle %d:%d-%d:%d-%d:%d\n",
		t1, pdb->atom[t1].atomNumber,
		t2, pdb->atom[t2].atomNumber,
		t3, pdb->atom[t3].atomNumber);
}

/*___________________________________________________________________________*/
/** print torsion */
__inline__ static void print_torsion(Str *pdb, int t1, int t2, int t3, int t4)
{
	fprintf(stderr, "torsion %d:%d-%d:%d-%d:%d-%d:%d\n",
		t1, pdb->atom[t1].atomNumber,
		t2, pdb->atom[t2].atomNumber,
		t3, pdb->atom[t3].atomNumber,
		t4, pdb->atom[t4].atomNumber);
}

/*___________________________________________________________________________*/
/** init topology */
void init_topology(Str *pdb, Topol *topol)
{
	unsigned int i;
#ifdef MPI
    int dimx = (int)ceil(pdb->nAtom / nodes); /* dimension of MPI loop */
    int nodeBondState[dimx]; /* vectors holding results of this node */
    int nodeNeighbourState[dimx];
    float nodeNeighbourPar[dimx];

	/* assuming an upper limit of 63 bonded interactions per atom */
	topol->bondState = alloc_mat2D_int(topol->bondState, (dimx * nodes), 64);
	/* assuming an upper limit of 255 non-bonded interactions per atom */
	topol->neighbourState = alloc_mat2D_int(topol->neighbourState, (dimx * nodes), 1024);
	/* assuming an upper limit of 255 neighbours per atom */
	topol->neighbourPar = alloc_mat2D_float(topol->neighbourPar, (dimx * nodes), 1024);
	topol->neighbourPar = alloc_mat2D_float(topol->neighbourPar, (dimx * nodes), 1024);

    for (i = 0; i < dimx; ++ i) {
        /* skip excess loops */
        /*if (((my_rank * dimx) + i) < pdb->nAtom) {*/
		nodeBondState[i] = 0;
		nodeNeighbourState[i] = 0;
		nodeNeighbourPar[i] = 0.;
	}

	/* communicate data between all nodes */
    MPI_Allgather(nodeBondState, dimx, MPI_FLOAT, topol->bondState, dimx, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Allgather(nodeNeighbourState, dimx, MPI_FLOAT, topol->neighbourState, dimx, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Allgather(nodeNeighbourPar, dimx, MPI_FLOAT, topol->neighbourPar, dimx, MPI_FLOAT, MPI_COMM_WORLD);
#else
	/* assuming an upper limit of 63 bonded interactions per atom */
	topol->bondState = alloc_mat2D_int(topol->bondState, pdb->nAtom, 64);
	/* assuming an upper limit of 255 non-bonded interactions per atom */
	topol->neighbourState = alloc_mat2D_int(topol->neighbourState, pdb->nAtom, 1024);
	/* assuming an upper limit of 255 neighbours per atom */
	topol->neighbourPar = alloc_mat2D_float(topol->neighbourPar, pdb->nAtom, 1024);

	/* maximally one interface nearest neighbour per atom */
	topol->interfaceNn = safe_malloc(pdb->nAtom * sizeof(int));
	topol->interfaceNnDist = safe_malloc(pdb->nAtom * sizeof(float));

	for (i = 0; i < pdb->nAtom; ++ i) {
		topol->bondState[i][0] = 0; /* no bonded pairs recorded */
		topol->neighbourState[i][0] = 0; /* no non-bonded pairs recorded */
		topol->neighbourPar[i][0] = 0.; /* no neighbours recorded */
		topol->interfaceNn[i] = -1; /* no nearest neighbour recorded */
		topol->interfaceNnDist[i] = FLT_MAX; /* distance to nearest neighbour */
	}
#endif
}

/*___________________________________________________________________________*/
/** free topology */
void free_topology(Str *pdb, Topol *topol)
{
#ifdef MPI
	int dimx = (int)ceil(pdb->nAtom / nodes);
	int dimMat2D = dimx * nodes;
#else
	int dimMat2D = pdb->nAtom;
#endif
	/* topol */
	free(topol->ib); /* bond */
	free(topol->jb);
	free(topol->it); /* angle */
	free(topol->jt);
	free(topol->kt);
	free(topol->ip); /* torsion */
	free(topol->jp);
	free(topol->kp);
	free(topol->lp);
	free(topol->in); /* non-bonded */
	free(topol->jn);

	free_mat2D_int(topol->bondState, dimMat2D); /* bond status of atom pairs */
	free_mat2D_int(topol->neighbourState, dimMat2D); /* neighbour status of atom pairs */
	free_mat2D_float(topol->neighbourPar,dimMat2D); /* POPS parameters of neighbours */
	free(topol->interfaceNn);
	free(topol->interfaceNnDist);
}

/*___________________________________________________________________________*/
/** print bond state */
void print_bondState(Str *pdb, Topol *topol)
{
	unsigned int i, j;

	for (i = 0; i < pdb->nAtom; ++ i) {
		fprintf(stderr, "atom %6d, bonds %2d, ", i, topol->bondState[i][0]);

		for (j = 1; j <= topol->bondState[i][0]; ++ j)
			fprintf(stderr, " %d", topol->bondState[i][j]);

		fprintf(stderr, "\n");
	}
}

/*____________________________________________________________________________*/
/** pairwise atom distance */
float atom_distance(Str *pdb, int i, int j)
{
	return v_rmsd(&(pdb->atom[i].pos), &(pdb->atom[j].pos));
}

/*____________________________________________________________________________*/
/** cutoff radius for non-bonded interaction */
float cutoff_radius(Type *type, ConstantSasa *constant_sasa, int i, int j, float rSolvent)
{
	return (constant_sasa[0].atomDataSasa[type->residueType[i]][type->atomType[i]].radius + \
		    constant_sasa[0].atomDataSasa[type->residueType[j]][type->atomType[j]].radius + \
			(2. * rSolvent));
}

/*____________________________________________________________________________*/
/** get bond state of two atoms */
__inline__ static int get_bondState(Topol *topol, int i, int j)
{
	unsigned int k, l;
	int state = 0;

	for (k = 1; k <= topol->bondState[i][0]; ++ k)
		if (topol->bondState[i][k] == j)
			return ++ state;

	for (l = 1; l <= topol->bondState[j][0]; ++ l)
		if (topol->bondState[j][l] == i)
			return ++ state;

	return state;
}

/*____________________________________________________________________________*/
/** calculate bonds */
int get_bonds(Str *pdb, Type *type, Topol *topol, ConstantSasa *constant_sasa, Argpdb *argpdb)
{
	unsigned int allocated = 64;
	float atomDistance; /*atom distance */
	float cutoffRadius; /* cutoff radius for bonded state */
	float cutoffFactor; /* pre-factor for cutoff radius calculation */

	/*___________________________________________________________________________*/
	/* MPI */
#ifdef MPI
	int xy;
	int *node_ib = safe_malloc(allocated * sizeof(int));
	int *node_jb = safe_malloc(allocated * sizeof(int));
	int node_nBond = 0;
	int nodes_nBond[nodes];
	/* total loop space is 1-dimensional square matrix of dimension (pdb->nAtom - 1) */
	int dimxy = (int)ceil((pdb->nAtom - 1) * (pdb->nAtom - 1) / nodes);

	/* for all pairwise atom combinations (as one-dimensional loop) */
	for (xy = 0; xy < dimxy; ++ xy) {
		/* to read in atom data, we need to refer to the original indices */
		/* matrix row index i */
		int i = floor(((my_rank * (pdb->nAtom - 1)) + xy) / (pdb->nAtom - 1));
		/* matrix column index j */
		int j = xy % (pdb->nAtom - 1);

		/* coarse grained 'P' needs a more generous cutoff */
		if (argpdb->coarse && (strncmp(pdb->atom[i].atomName, " P  ", 4) == 0))
			cutoffFactor = 0.7;
		else
			cutoffFactor = 0.5;

		/* if atoms i,j in the same or proximate residue */
		if (((pdb->atom[j].residueNumber == pdb->atom[i].residueNumber) || \
			(pdb->atom[j].residueNumber == pdb->atom[i].residueNumber + 1)) && \
			 strcmp(pdb->atom[j].chainIdentifier, pdb->atom[i].chainIdentifier) == 0) {

			/* add bond if atom distance shorter than cutoff */
			/* atoms bonded if dist =< 0.5 * (atomRadius_i + atomRadius_j) */
			atomDistance = atom_distance(pdb, i, j);
			cutoffRadius = (cutoffFactor * \
							(constant_sasa[0].atomDataSasa[type->residueType[i]][type->atomType[i]].radius + \
							 constant_sasa[0].atomDataSasa[type->residueType[j]][type->atomType[j]].radius));

			if (atomDistance < cutoffRadius) {
				/* assign arrays of bonded atoms ib-jb */
				node_ib[node_nBond] = i;
				node_jb[node_nBond] = j;

				++ node_nBond; /* increment bond index */

				/* record 1,2-bond in bondState matrix */
				/* first index is atom number, second index is mixed:
					array element '0' records the total number of bonded atoms,
					then each bonded atom's ID number is added to that array position;
					the result is a list of all bonded atom IDs and their total number at '0' */
/*
				++ topol->bondState[i][0];
				topol->bondState[i][topol->bondState[i][0]] = j;
				++ topol->bondState[j][0];
				topol->bondState[j][topol->bondState[j][0]] = i;
*/
				/*print_pair(pdb, i, j);*/

				/* add memory if needed */ 
				if (node_nBond == allocated) {
					allocated += 64;
					node_ib = safe_realloc(node_ib, allocated * sizeof(int));
					node_jb = safe_realloc(node_jb, allocated * sizeof(int));
				}

				/* warn if atoms too close */
				if (atomDistance < 0.5)
					fprintf(stderr, "Warning: Atoms %d %d too close\n", i, j);
			}
		}
	}

	MPI_Allgather(&node_nBond, 1, MPI_INT, nodes_nBond, 1, MPI_INT, MPI_COMM_WORLD);

	/* assert that allocation is equal across nodes,
		otherwise make it equal */
	MPI_Allgather(node_ib, 1, MPI_INT, nodes_ib, 1, MPI_INT, MPI_COMM_WORLD);
	MPI_Allgather(node_jb, 1, MPI_INT, nodes_jb, 1, MPI_INT, MPI_COMM_WORLD);

	/* construct 'topol' arrays from results on nodes */

	/*___________________________________________________________________________*/
	/* non-MPI */
#else
	unsigned int i, j;
	/* allocate memory */
	topol->ib = safe_malloc(allocated * sizeof(int));
	topol->jb = safe_malloc(allocated * sizeof(int));
	topol->nBond = 0;

	/* for all pairwise atom combinations */
	for (i = 0; i < pdb->nAtom - 1; ++ i) {
		/* coarse grained 'P' needs more generous cutoff */
		if (argpdb->coarse && (strncmp(pdb->atom[i].atomName, " P  ", 4) == 0))
			cutoffFactor = 0.7;
		else
			cutoffFactor = 0.5;
		for (j = i + 1; j < pdb->nAtom; ++ j) {

			/* if atoms i,j in the same or proximate residue */
			if (((pdb->atom[j].residueNumber == pdb->atom[i].residueNumber) || \
				(pdb->atom[j].residueNumber == pdb->atom[i].residueNumber + 1)) && \
				 strcmp(pdb->atom[j].chainIdentifier, pdb->atom[i].chainIdentifier) == 0) {

				/* add bond if atom distance shorter than cutoff */
				/* atoms bonded if dist =< 0.5 * (atomRadius_i + atomRadius_j) */
				atomDistance = atom_distance(pdb, i, j);
				cutoffRadius = (cutoffFactor * \
								(constant_sasa[0].atomDataSasa[type->residueType[i]][type->atomType[i]].radius + \
					             constant_sasa[0].atomDataSasa[type->residueType[j]][type->atomType[j]].radius));

				if (atomDistance < cutoffRadius) {
					/* assign arrays of bonded atoms ib-jb */
					topol->ib[topol->nBond] = i;
					topol->jb[topol->nBond] = j;

					++ topol->nBond; /* increment bond index */

					/* record 1,2-bond in bondState matrix */
					/* first index is atom number, second index is mixed:
						array element '0' records the total number of bonded atoms,
						then each bonded atom's ID number is added to that array position;
						the result is a list of all bonded atom IDs and their total number at '0' */
					++ topol->bondState[i][0];
					topol->bondState[i][topol->bondState[i][0]] = j;
					++ topol->bondState[j][0];
					topol->bondState[j][topol->bondState[j][0]] = i;

					/*print_pair(pdb, i, j);*/

					/* add memory if needed */ 
					if (topol->nBond == allocated) {
						allocated += 64;
						topol->ib = safe_realloc(topol->ib, allocated * sizeof(int));
						topol->jb = safe_realloc(topol->jb, allocated * sizeof(int));
					}

					/* warn if atoms too close */
					if (atomDistance < 0.5)
						fprintf(stderr, "Warning: Atoms %d %d too close\n", i, j);
				}
			}
        }
    }
#endif

	return(0);
}

/*___________________________________________________________________________*/
/** calculate angles (from bonds) */
/** this routine tests for pairs of bonds with an identical atom;
 * if that is the case, one of the bond ends with the identical atom 
 * is skipped and the other three atoms are defined as forming the bond angle */
int get_angles(Str *pdb, Topol *topol)
{
	unsigned int i, j;
	unsigned int allocated = 64;

	/* allocate memory */
	topol->it = safe_malloc(allocated * sizeof(int));
	topol->jt = safe_malloc(allocated * sizeof(int));
	topol->kt = safe_malloc(allocated * sizeof(int));

    topol->nAngle = 0;

	/** for all pairwise bond combinations */
    for (i = 0; i < topol->nBond - 1; ++ i) {
        for (j = i + 1; j < topol->nBond; ++ j) {
			/* add memory if needed */ 
			if (topol->nAngle == allocated) {
				allocated += 64;
				topol->it = safe_realloc(topol->it, allocated * sizeof(int));
				topol->jt = safe_realloc(topol->jt, allocated * sizeof(int));
				topol->kt = safe_realloc(topol->kt, allocated * sizeof(int));
			}

			/*____________________________________________________________________________*/
			/* 1 */
			/** shared atoms between bonds: first atom of bond i, first atom of bond j */
            if (topol->ib[i] == topol->ib[j]) {
				/* assign arrays of bonded atoms it-jt-kt */
				topol->it[topol->nAngle] = topol->jb[i];
                topol->jt[topol->nAngle] = topol->ib[i];
                topol->kt[topol->nAngle] = topol->jb[j];
				/*print_angle(pdb, topol->jb[i], topol->ib[i], topol->jb[j]);*/

				/* record 1,3-bond in bondState matrix */
				++ topol->bondState[topol->jb[i]][0];
				topol->bondState[topol->jb[i]][topol->bondState[topol->jb[i]][0]] = topol->jb[j];
				++ topol->bondState[topol->jb[j]][0];
				topol->bondState[topol->jb[j]][topol->bondState[topol->jb[j]][0]] = topol->jb[i];

                ++ topol->nAngle; /* increment angle index */
				continue;
            }

			/*____________________________________________________________________________*/
			/* 2 */
			/** shared atoms between bonds: first atom of bond i, second atom of bond j */
            if (topol->ib[i] == topol->jb[j]) {
				/* assign arrays of bonded atoms it-jt-kt */
				topol->it[topol->nAngle] = topol->jb[i];
                topol->jt[topol->nAngle] = topol->ib[i];
                topol->kt[topol->nAngle] = topol->ib[j];
				/*print_angle(pdb, topol->jb[i], topol->ib[i], topol->ib[j]);*/

				/* record 1,3-bond in bondState matrix */
				++ topol->bondState[topol->jb[i]][0];
				topol->bondState[topol->jb[i]][topol->bondState[topol->jb[i]][0]] = topol->ib[j];
				++ topol->bondState[topol->ib[j]][0];
				topol->bondState[topol->ib[j]][topol->bondState[topol->ib[j]][0]] = topol->jb[i];

                ++ topol->nAngle; /* increment angle index */
				continue;
            }

			/*____________________________________________________________________________*/
			/* 3 */
			/** shared atoms between bonds: second atom of bond i, first atom of bond j */
            if (topol->jb[i] == topol->ib[j]) {
				/* assign arrays of bonded atoms it-jt-kt */
                topol->it[topol->nAngle] = topol->ib[i];
                topol->jt[topol->nAngle] = topol->jb[i];
                topol->kt[topol->nAngle] = topol->jb[j];
				/*print_angle(pdb, topol->ib[i], topol->jb[i], topol->jb[j]);*/

				/* record 1,3-bond in bondState matrix */
				++ topol->bondState[topol->ib[i]][0];
				topol->bondState[topol->ib[i]][topol->bondState[topol->ib[i]][0]] = topol->jb[j];
				++ topol->bondState[topol->jb[j]][0];
				topol->bondState[topol->jb[j]][topol->bondState[topol->jb[j]][0]] = topol->ib[i];

                ++ topol->nAngle; /* increment angle index */
				continue;
            }

			/*____________________________________________________________________________*/
			/* 4 */
			/** shared atoms between bonds: second atom of bond i, second atom of bond j */
            if (topol->jb[i] == topol->jb[j]) {
				/* assign arrays of bonded atoms it-jt-kt */
                topol->it[topol->nAngle] = topol->ib[i];
                topol->jt[topol->nAngle] = topol->jb[i];
                topol->kt[topol->nAngle] = topol->ib[j];
				/*print_angle(pdb, topol->ib[i], topol->jb[i], topol->ib[j]);*/

				/* record 1,3-bond in bondState matrix */
				++ topol->bondState[topol->ib[i]][0];
				topol->bondState[topol->ib[i]][topol->bondState[topol->ib[i]][0]] = topol->ib[j];
				++ topol->bondState[topol->ib[j]][0];
				topol->bondState[topol->ib[j]][topol->bondState[topol->ib[j]][0]] = topol->ib[i];

                ++ topol->nAngle; /* increment angle index */
				continue;
            }
        }
    }

	return(0);
}

/*___________________________________________________________________________*/
/** calculate torsions (from angles) */
/** this routine tests for pairs of angles with two identical atoms;
 * if that is the case, two of the identical (redundant) atoms 
 * are skipped and the other four atoms are defined as forming the torsion angle */
/** Checking for ring topologies by triple-angle redundancy checks has been
 * replaced by ring annotation in 'sasa_data.h'. The redundancy check has
 * cubic complexity and was by far the slowest step of this program.
 * Other improvements are:
 * 1. Some torsions in the TRP ring were not excluded,
 * because the endpoints could not be linked by an angle or torsion.
 * 2. Some torsions in any six-membered ring were not excluded,
 * because the linking torsion may not yet be defined at that point in
 * the algorithm.
 * Therefore, this program gives slightly different (but correct) results
 * if compared to older POPS versions.
 * The redundancy check should be performed for all HET residues. */
int get_torsions(Str *pdb, Type *type, Topol *topol, ConstantSasa *constant_sasa)
{
	unsigned int i, j/*, k*/;
	unsigned int allocated = 64;
	/*int angle_redundancy, torsion_redundancy;*/
	int ring_i, ring_j;
	
	/* allocate memory */
	topol->ip = safe_malloc(allocated * sizeof(int));
	topol->jp = safe_malloc(allocated * sizeof(int));
	topol->kp = safe_malloc(allocated * sizeof(int));
	topol->lp = safe_malloc(allocated * sizeof(int));

    topol->nTorsion = 0;
	
	/* for all angle pair combinations */
	for (i = 0; i < topol->nAngle - 1; ++ i) {
		/* ring status of angle 'i' */
		ring_i = constant_sasa[0].atomDataSasa[type->residueType[topol->it[i]]][type->atomType[topol->it[i]]].ring + \
				 constant_sasa[0].atomDataSasa[type->residueType[topol->jt[i]]][type->atomType[topol->jt[i]]].ring + \
				 constant_sasa[0].atomDataSasa[type->residueType[topol->kt[i]]][type->atomType[topol->kt[i]]].ring;

		for (j = i + 1; j < topol->nAngle; ++ j) {
			/* ring status of angle 'j' */
			ring_j = constant_sasa[0].atomDataSasa[type->residueType[topol->it[j]]][type->atomType[topol->it[j]]].ring + \
					 constant_sasa[0].atomDataSasa[type->residueType[topol->jt[j]]][type->atomType[topol->jt[j]]].ring + \
					 constant_sasa[0].atomDataSasa[type->residueType[topol->kt[j]]][type->atomType[topol->kt[j]]].ring;


			/* if angles 'i' and 'j' entirely in ring, skip torsion assignment */
			if ((ring_i == 3) && (ring_j == 3))
				continue;

			/* add memory if needed */ 
			if (topol->nTorsion == allocated) {
				allocated += 64;
				topol->ip = safe_realloc(topol->ip, allocated * sizeof(int));
				topol->jp = safe_realloc(topol->jp, allocated * sizeof(int));
				topol->kp = safe_realloc(topol->kp, allocated * sizeof(int));
				topol->lp = safe_realloc(topol->lp, allocated * sizeof(int));
			}

			/*____________________________________________________________________________*/
			/* 1 */
			/** shared atoms between angles:
				first atom of angle i, second atom of angle j,
				second atom of angle i, first atom of angle j */
            if ((topol->it[i] == topol->jt[j]) && (topol->jt[i] == topol->it[j])) {

				/*___________________________________________________________________________*/
				/* assign arrays of bonded atoms kt[j]-it[i]-jt[i]-kt[i] */
				topol->ip[topol->nTorsion] = topol->kt[j];
				topol->jp[topol->nTorsion] = topol->it[i];
				topol->kp[topol->nTorsion] = topol->jt[i];
				topol->lp[topol->nTorsion] = topol->kt[i];

				/*print_torsion(pdb, topol->kt[j], topol->it[i], topol->jt[i], topol->kt[i]);*/

				/* record 1,4-bond in bondState matrix */
				++ topol->bondState[topol->kt[i]][0];
				topol->bondState[topol->kt[i]][topol->bondState[topol->kt[i]][0]] = topol->kt[j];
				++ topol->bondState[topol->kt[j]][0];
				topol->bondState[topol->kt[j]][topol->bondState[topol->kt[j]][0]] = topol->kt[i];

				++ topol->nTorsion; /* increment torsion index */
				continue;
			}

			/*____________________________________________________________________________*/
            /* 2 */
			/** shared atoms between angles:
				first atom of angle i, second atom of angle j,
				second atom of angle i, third atom of angle j */
            if (((topol->it[i] == topol->jt[j])) && ((topol->jt[i] == topol->kt[j]))) {

				/*___________________________________________________________________________*/
				/* assign arrays of bonded atoms it[j]-it[i]-jt[i]-kt[i] */
				topol->ip[topol->nTorsion] = topol->it[j];
				topol->jp[topol->nTorsion] = topol->it[i];
				topol->kp[topol->nTorsion] = topol->jt[i];
				topol->lp[topol->nTorsion] = topol->kt[i];

				/*print_torsion(pdb, topol->it[j], topol->it[i], topol->jt[i], topol->kt[i]);*/

				/* record 1,4-bond in bondState matrix */
				++ topol->bondState[topol->kt[i]][0];
				topol->bondState[topol->kt[i]][topol->bondState[topol->kt[i]][0]] = topol->it[j];
				++ topol->bondState[topol->ib[j]][0];
				topol->bondState[topol->it[j]][topol->bondState[topol->it[j]][0]] = topol->kt[i];

				++ topol->nTorsion; /* increment torsion index */
				continue;
            }

			/*____________________________________________________________________________*/
            /* 3 */
			/** shared atoms between angles:
				third atom of angle i, second atom of angle j,
				second atom of angle i, first atom of angle j */
            if ((topol->kt[i] == topol->jt[j]) && (topol->jt[i] == topol->it[j])) {

				/*___________________________________________________________________________*/
				/* assign arrays of bonded atoms it[j]-it[i]-jt[i]-kt[j] */
				topol->ip[topol->nTorsion] = topol->it[i];
				topol->jp[topol->nTorsion] = topol->jt[i];
				topol->kp[topol->nTorsion] = topol->kt[i];
				topol->lp[topol->nTorsion] = topol->kt[j];

				/*print_torsion(pdb, topol->it[i], topol->jt[i], topol->kt[i], topol->kt[j]);*/

				/* record 1,4-bond in bondState matrix */
				++ topol->bondState[topol->it[i]][0];
				topol->bondState[topol->it[i]][topol->bondState[topol->it[i]][0]] = topol->kt[j];
				++ topol->bondState[topol->kt[j]][0];
				topol->bondState[topol->kt[j]][topol->bondState[topol->kt[j]][0]] = topol->it[i];

				++ topol->nTorsion; /* increment torsion index */
				continue;
			}

			/*____________________________________________________________________________*/
            /* 4 */
			/** shared atoms between angles:
				third atom of angle i, second atom of angle j,
				second atom of angle i, third atom of angle j */
            if ((topol->kt[i] == topol->jt[j]) && (topol->jt[i] == topol->kt[j])) {

				/*___________________________________________________________________________*/
				/* assign arrays of bonded atoms it[i]-jt[i]-kt[i]-it[j] */
				topol->ip[topol->nTorsion] = topol->it[i];
				topol->jp[topol->nTorsion] = topol->jt[i];
				topol->kp[topol->nTorsion] = topol->kt[i];
				topol->lp[topol->nTorsion] = topol->it[j];

				/*print_torsion(pdb, topol->it[i], topol->jt[i], topol->kt[i], topol->it[j]);*/

				/* record 1,4-bond in bondState matrix */
				++ topol->bondState[topol->it[i]][0];
				topol->bondState[topol->it[i]][topol->bondState[topol->it[i]][0]] = topol->it[j];
				++ topol->bondState[topol->it[j]][0];
				topol->bondState[topol->it[j]][topol->bondState[topol->it[j]][0]] = topol->it[i];

				++ topol->nTorsion; /* increment torsion index */
				continue;
			}
        }
    }

	return(0);
}

/*___________________________________________________________________________*/
/** calculate non-bonded overlapping atoms */ 
/** atoms overlapping if dist < RADATM(i) + RADATM(j) + 2*RSOLV */
int nonbonded_overlaps(Str *pdb, Type *type, Topol *topol, ConstantSasa *constant_sasa, Arg *arg)
{
	unsigned int i, j;
	unsigned int allocated = 64;
	float atomDistance;
	float cutoffRadius;
	int bondState;

	/* allocate memory */
	topol->in = safe_malloc(allocated * sizeof(int));
	topol->jn = safe_malloc(allocated * sizeof(int));

    topol->nNonBonded = 0;

	/* for all atom pair combinations */
    for (i = 0; i < pdb->nAtom - 1; ++ i) {
		for (j = i + 1; j < pdb->nAtom; ++ j) {
			atomDistance = atom_distance(pdb, i, j);
			cutoffRadius = cutoff_radius(type, constant_sasa, i, j, arg->rProbe);
			bondState = get_bondState(topol, i, j);

			/* if atoms closer than non-bonded cutoff and not bonded */
			if ((atomDistance < cutoffRadius) && (bondState == 0)) {
				/*___________________________________________________________________________*/
				/* assign arrays of non-bonded atoms */
				topol->in[topol->nNonBonded] = i;
				topol->jn[topol->nNonBonded] = j;

				/*print_pair(pdb, i, j);*/

				/* record non-bonded pair in neighbourState matrix */
				++ topol->neighbourState[i][0];
				assert(topol->neighbourState[i][0] < 1024);
				topol->neighbourState[i][topol->neighbourState[i][0]] = j;
				++ topol->neighbourState[j][0];
				assert(topol->neighbourState[j][0] < 1024);
				topol->neighbourState[j][topol->neighbourState[j][0]] = i;

				/* record nearest neighbour */
				if (atomDistance < topol->interfaceNnDist[i]) {
					topol->interfaceNnDist[i] = atomDistance;
					topol->interfaceNn[i] = j;
				}
				if (atomDistance < topol->interfaceNnDist[j]) {
					topol->interfaceNnDist[j] = atomDistance;
					topol->interfaceNn[j] = i;
				}

				++ topol->nNonBonded; /* increment non-bonded atom index */

				/*___________________________________________________________________________*/
				/* add memory if needed */ 
				if (topol->nNonBonded == allocated) {
					allocated += 64;
					topol->in = safe_realloc(topol->in, allocated * sizeof(int));
					topol->jn = safe_realloc(topol->jn, allocated * sizeof(int));
				}
			}
		}
	}

	return(0);
}

/*____________________________________________________________________________*/
/** derive molecular topology */
int get_topology(Str *pdb, Type *type, Topol *topol, ConstantSasa *constant_sasa, Argpdb *argpdb, Arg *arg)
{
	if (pdb->nAtom < 2)
		Error("Need at least 2 atoms! Try to run POPS* using the '--coarse' switch.");

	/* the bondState matrix records, which atom pairs are bonded 
	   (bond distances 1,2 1,3 1,4); the rest has non-bonded interactions 
		if distance < cutoffRadius */

#if DEBUG>1
	fprintf(stderr, "%s:%d: nAtom = %d\n", __FILE__, __LINE__, pdb->nAtom);
#endif

	get_bonds(pdb, type, topol, constant_sasa, argpdb); /* calculate bonds (from atoms) */
#if DEBUG>1
	fprintf(stderr, "%s:%d: nBond = %d\n", __FILE__, __LINE__, topol->nBond);
#endif
	if (topol->nBond < 2)
		Error("Need at least 2 bonds! Try to run POPS* using the '--coarse' switch.");

	get_angles(pdb, topol); /* calculate angles (from bonds) */
#if DEBUG>1
	fprintf(stderr, "%s:%d: nAngle = %d\n", __FILE__, __LINE__, topol->nAngle);
#endif
	if (topol->nAngle < 2)
		Error("Need at least 2 angles! Try to run POPS* using the '--coarse' switch.\n");

	get_torsions(pdb, type, topol, constant_sasa); /* calculate torsions (from angles) */
#if DEBUG>1
	fprintf(stderr, "%s:%d: nTorsion = %d\n", __FILE__, __LINE__, topol->nTorsion);
#endif
	if (topol->nTorsion < 2)
		Error("Need at least 2 torsions! Try to run POPS* using the '--coarse' switch.\n");

	nonbonded_overlaps(pdb, type, topol, constant_sasa, arg); /* calculate overlapping atoms */
#if DEBUG>1
	fprintf(stderr, "%s:%d: nNonBonded = %d\n", __FILE__, __LINE__, topol->nNonBonded);
#endif

	/*print_bondState(pdb, topol);*/

	return 0;
}

