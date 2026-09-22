/*=============================================================================
topol : molecular topology
Copyright (C) 2002 Franca Fraternali
Copyright (C) 2008 Jens Kleinjung
Copyright (C) 2002 Luigi Cavallo
Copyright (C) 2002 Kuang Lin and Valerie Hindie
Read the COPYING file for license information.
=============================================================================*/

#include "config.h"
#include "topol.h"

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
/** record a bonded (1-2, 1-3 or 1-4) atom pair in the bondState matrix */
__inline__ static void add_bondState(Topol *topol, int a, int b)
{
	if (topol->bondState[a][0] >= 63 || topol->bondState[b][0] >= 63)
		Error("More than 63 bonded interactions per atom: check input structure for overlapping atoms");

	++ topol->bondState[a][0];
	topol->bondState[a][topol->bondState[a][0]] = b;
	++ topol->bondState[b][0];
	topol->bondState[b][topol->bondState[b][0]] = a;
}

/*___________________________________________________________________________*/
/** init topology */
void init_topology(Arg *arg, Str *pdb, Topol *topol)
{
	unsigned int i, j;
	/* assuming an upper limit of 63 bonded interactions per atom */
	topol->bondState = alloc_mat2D_int(topol->bondState, pdb->nAtom, 64);
	/* assuming an upper limit of 255 non-bonded interactions per atom */
	topol->neighbourState = alloc_mat2D_int(topol->neighbourState, pdb->nAtom, 1024);
	/* assuming an upper limit of 255 neighbours per atom */
	topol->neighbourPar = alloc_mat2D_float(topol->neighbourPar, pdb->nAtom, 1024);

	/* maximally one interface nearest neighbour per atom */
	topol->interfaceNn = safe_malloc(pdb->nAtom * sizeof(int));
	topol->interfaceNnDist = safe_malloc(pdb->nAtom * sizeof(float));

	/* between-chain or between domain Calpha distance matrix */
	/* determine number of Calpha atoms per chain or domain */
	/* rows: CA atoms of the first chain, columns: CA atoms of the second chain */
	for (i = 0, topol->nCA1 = 0, topol->nCA2 = 0; i < pdb->nAtom; ++ i) {
		if (strcmp(pdb->atom[i].atomName, "CA") == 0) {
			if (pdb->atom[i].chainIndex == 0)
				++ topol->nCA1;
			else if (pdb->atom[i].chainIndex == 1)
				++ topol->nCA2;
		}
	}
	if (! arg->silent && pdb->nChain == 2)
		fprintf(stdout, "CA distance matrix has dimensions %d x %d \n", topol->nCA1, topol->nCA2);
	/* allocate Calpha distance matrix */
	topol->distMatCA = alloc_mat2D_float(topol->distMatCA, topol->nCA1, topol->nCA2);
	init_mat2D_float(topol->distMatCA, topol->nCA1, topol->nCA2, 0.);
	/* array of residue numbers for each chain */
	topol->resCA1 = safe_malloc(topol->nCA1 * sizeof(int));
	for (i = 0; i < topol->nCA1; ++ i) {
		topol->resCA1[i] = -999;
	}
	topol->resCA2 = safe_malloc(topol->nCA2 * sizeof(int));
	for (j = 0; j < topol->nCA2; ++ j) {
		topol->resCA2[j] = -999;
	}

	for (i = 0; i < pdb->nAtom; ++ i) {
		topol->bondState[i][0] = 0; /* no bonded pairs recorded */
		topol->neighbourState[i][0] = 0; /* no non-bonded pairs recorded */
		topol->neighbourPar[i][0] = 0.; /* no neighbours recorded */
		topol->interfaceNn[i] = -1; /* no nearest neighbour recorded */
		topol->interfaceNnDist[i] = FLT_MAX; /* distance to nearest neighbour */
	}
}

/*___________________________________________________________________________*/
/** free topology */
void free_topology(Str *pdb, Topol *topol)
{
	int dimMat2D = pdb->nAtom;

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
	free_mat2D_float(topol->neighbourPar, dimMat2D); /* POPS parameters of neighbours */
	free(topol->interfaceNn);
	free(topol->interfaceNnDist);
	free_mat2D_float(topol->distMatCA, topol->nCA1);
	free(topol->resCA1);
	free(topol->resCA2);
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
	unsigned int i, j;
	/* allocate memory */
	topol->ib = safe_malloc(allocated * sizeof(int));
	topol->jb = safe_malloc(allocated * sizeof(int));
	topol->nBond = 0;

	/* for all pairwise atom combinations */
	for (i = 0; i < pdb->nAtom - 1; ++ i) {
		for (j = i + 1; j < pdb->nAtom; ++ j) {
			/* coarse grained 'P' needs more generous cutoff */
			if (argpdb->coarse && ((strcmp(pdb->atom[i].atomName, "P") == 0) ||
									(strcmp(pdb->atom[j].atomName, "P") == 0)))
				cutoffFactor = 0.7;
			else
				cutoffFactor = 0.5;

			/* if atoms i,j in the same or the next residue of the same chain;
				'next' is the position in the input (residueIndex), not the residue number,
				so that numbering gaps, insertion codes and descending numbering
				do not break the chain topology; a real chain break fails the distance test */
			if ((pdb->atom[j].chainIndex == pdb->atom[i].chainIndex) &&
				((pdb->atom[j].residueIndex == pdb->atom[i].residueIndex) ||
				 (pdb->atom[j].residueIndex == pdb->atom[i].residueIndex + 1))) {

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
					add_bondState(topol, i, j);

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
				add_bondState(topol, topol->jb[i], topol->jb[j]);

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
				add_bondState(topol, topol->jb[i], topol->ib[j]);

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
				add_bondState(topol, topol->ib[i], topol->jb[j]);

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
				add_bondState(topol, topol->ib[i], topol->ib[j]);

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
				add_bondState(topol, topol->kt[i], topol->kt[j]);

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
				add_bondState(topol, topol->kt[i], topol->it[j]);

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
				add_bondState(topol, topol->it[i], topol->kt[j]);

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
				add_bondState(topol, topol->it[i], topol->it[j]);

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
				if (topol->neighbourState[i][0] >= 1023 || topol->neighbourState[j][0] >= 1023)
					Error("More than 1023 neighbours per atom: check input structure for overlapping atoms");
				++ topol->neighbourState[i][0];
				topol->neighbourState[i][topol->neighbourState[i][0]] = j;
				++ topol->neighbourState[j][0];
				topol->neighbourState[j][topol->neighbourState[j][0]] = i;

				/* record nearest neighbour on a different chain */
				if (pdb->atom[i].chainIndex != pdb->atom[j].chainIndex &&
					atomDistance < topol->interfaceNnDist[i]) {
					topol->interfaceNnDist[i] = atomDistance;
					topol->interfaceNn[i] = j;
				}
				if (pdb->atom[i].chainIndex != pdb->atom[j].chainIndex &&
					atomDistance < topol->interfaceNnDist[j]) {
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
		exit_structure_error(arg, pdb, "< 2 atoms: try '--coarse' switch");

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
		exit_structure_error(arg, pdb, "< 2 bonds: try '--coarse' switch");

	get_angles(pdb, topol); /* calculate angles (from bonds) */
#if DEBUG>1
	fprintf(stderr, "%s:%d: nAngle = %d\n", __FILE__, __LINE__, topol->nAngle);
#endif
	if (topol->nAngle < 2)
		exit_structure_error(arg, pdb, "< 2 angles: try '--coarse' switch");

	get_torsions(pdb, type, topol, constant_sasa); /* calculate torsions (from angles) */
#if DEBUG>1
	fprintf(stderr, "%s:%d: nTorsion = %d\n", __FILE__, __LINE__, topol->nTorsion);
#endif
	if (topol->nTorsion < 2)
		exit_structure_error(arg, pdb, "< 2 torsions: try '--coarse' switch");

	nonbonded_overlaps(pdb, type, topol, constant_sasa, arg); /* calculate overlapping atoms */
#if DEBUG>1
	fprintf(stderr, "%s:%d: nNonBonded = %d\n", __FILE__, __LINE__, topol->nNonBonded);
#endif

	/*print_bondState(pdb, topol);*/

	return 0;
}

/*____________________________________________________________________________*/
/** Calpha distances between different chains */
/** The raw distances are converted in residue-size-weighted soft distances */
int calpha_distances(Arg *arg, Str *pdb, Topol *topol, ConstantSasa *res_sasa) {
	unsigned int i, j, k;
	int n_cai, n_caj;
	float distance = 0.;
	float soft_distance = 0.;
	int found_CA1 = 0;
	int found_CA2 = 0;
	float r_CA1 = 0.;
	float r_CA2 = 0.;
	float soft_distance_norm = 0.;

	n_cai = 0;
	n_caj = 0;

	/* rows: CA atoms of the first chain */
	for (i = 0; i < pdb->nAtom; ++i) {
    	if (strcmp(pdb->atom[i].atomName, "CA") != 0 || pdb->atom[i].chainIndex != 0)
        	continue;

    	n_caj = 0;

		/* columns: CA atoms of the second chain */
    	for (j = 0; j < pdb->nAtom; ++j) {
			if (strcmp(pdb->atom[j].atomName, "CA") != 0 || pdb->atom[j].chainIndex != 1)
            	continue;

			distance = atom_distance(pdb, i, j);
			soft_distance = 1 / (1 + exp((distance - 8) / 0.5));
			
			/* get residue radius (coarse grained) */
			found_CA1 = found_CA2 = 0;
			r_CA1 = r_CA2 = 0.;
			for (k = 0; k < res_sasa->nResidueType; ++ k) {
 				if (strcmp(pdb->atom[i].residueName, res_sasa->atomDataSasa[k][0].residueName) == 0) {
        			/* found residue of atom i at index k */
					r_CA1 = res_sasa->atomDataSasa[k][0].radius;
					++ found_CA1;
				}

 				if (strcmp(pdb->atom[j].residueName, res_sasa->atomDataSasa[k][0].residueName) == 0) {
        			/* found residue of atom j at index k */
					r_CA2 = res_sasa->atomDataSasa[k][0].radius;
					++ found_CA2;
        		}
			
				if (found_CA1 && found_CA2) {
					break;
    			}
			}

			/* the mean radius across 20 amino acids is 4.35 */
			soft_distance_norm = soft_distance *  (sqrt(r_CA1 * r_CA2) / 4.35);
        	topol->distMatCA[n_cai][n_caj] = soft_distance_norm;
			
			topol->resCA1[n_cai] = pdb->atom[i].residueNumber;
			topol->resCA2[n_caj] = pdb->atom[j].residueNumber;		

        	++n_caj;
    	}

    	++n_cai;
	}

	if (n_cai != topol->nCA1 || (n_cai > 0 && n_caj != topol->nCA2))
		ErrorLoc("Inconsistent CA atom count", __FILE__, __LINE__);

	return(0);
}

/*____________________________________________________________________________*/
/*____________________________________________________________________________*/

