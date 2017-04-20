/*
  A simple 2D hydro code
  (C) Romain Teyssier : CEA/IRFU           -- original F90 code
  (C) Pierre-Francois Lavallee : IDRIS      -- original F90 code
  (C) Guillaume Colin de Verdiere : CEA/DAM -- for the C version
*/

#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>

#include "parametres.h"
#include "hydro_godunov.h"
#include "hydro_funcs.h"
#include "utils.h"
#include "make_boundary.h"

#include "cmpflx.h"
#include "conservar.h"
#include "equation_of_state.h"
#include "qleftright.h"
#include "constoprim.h"
#include "riemann.h"
#include "trace.h"
#include "slope.h"

static void
Dmemset(double *dq, double motif, size_t nbr)
{
    long i;
    for (i = 0; i < nbr; i++) {
        dq[i] = motif;
    }
}

void
hydro_godunov(long idim, double dt, const hydroparam_t H, hydrovar_t * Hv,
              hydrowork_t * Hw, hydrovarwork_t * Hvw)
{

    // Local variables
    long i, j;
    double dtdx;

    double *dq;
    double *e;
    double *flux;
    double *q;
    double *qleft, *qright;
    double *qxm, *qxp;
    double *u;
    double *c;
    double *uold;
    double *rl, *ul, *pl, *cl, *wl, *rr, *ur, *pr, *cr, *wr, *ro, *uo, *po, *co, *wo,
           *rstar, *ustar, *pstar, *cstar, *spin, *spout, *ushock, *frac, *scr, *delp,
           *pold;
    long   *sgnm, *ind, *ind2;
    double *qgdnv;
    double *qID;
    double *qIP;

    WHERE("hydro_godunov");

    // constant
    dtdx = dt / H.dx;

#pragma omp single
{
    // Update boundary conditions
    if (H.prt) {
      fprintf(stdout, "godunov %ld\n", idim);
      PRINTUOLD(H, Hv);
      }
    make_boundary(idim, H, Hv);
    //PRINTUOLD(H, Hv);
} // End omp single
    
    // Allocate work space for 1D sweeps
    allocate_work_space(H, Hw, Hvw);
    
    uold = Hv->uold;
    qgdnv = Hvw->qgdnv;
    flux = Hvw->flux;
    c = Hw->c;
    q = Hvw->q;
    e = Hw->e;
    u = Hvw->u;
    qxm = Hvw->qxm;
    qxp = Hvw->qxp;
    qleft = Hvw->qleft;
    qright = Hvw->qright;
    dq = Hvw->dq;
    rl = Hw->rl;
    ul = Hw->ul;
    pl = Hw->pl;
    cl = Hw->cl;
    wl = Hw->wl;
    rr = Hw->rr;
    ur = Hw->ur;
    pr = Hw->pr;
    cr = Hw->cr;
    wr = Hw->wr;
    ro = Hw->ro;
    uo = Hw->uo;
    po = Hw->po;
    co = Hw->co;
    wo = Hw->wo;
    rstar = Hw->rstar;
    ustar = Hw->ustar;
    pstar = Hw->pstar;
    cstar = Hw->cstar;
    sgnm = Hw->sgnm;
    spin = Hw->spin;
    spout = Hw->spout;
    ushock = Hw->ushock;
    frac = Hw->frac;
    scr = Hw->scr;
    delp = Hw->delp;
    pold = Hw->pold;
    ind = Hw->ind;
    ind2 = Hw->ind2;

    if (idim == 1) {
#pragma omp for schedule(dynamic,2)
      for (j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++) {
	qID = &q[IHvw(0, ID)];
	qIP = &q[IHvw(0, IP)];
	
	gatherConservativeVars(idim, j, uold, u, H.imin, H.imax, H.jmin,
			       H.jmax, H.nvar, H.nxt, H.nyt, H.nxyt);
	
	// Convert to primitive variables
	constoprim(u, q, e, H.nxt, H.nxyt, H.nvar, H.smallr);
	
	equation_of_state(qID, e, qIP, c, 0, H.nxt, H.smallc, H.gamma);
	
	Dmemset(dq, 0, (H.nxyt + 2) * H.nvar);
	
	// Characteristic tracing
	if (H.iorder != 1) {
	  slope(q, dq, H.nxt, H.nvar, H.nxyt, H.slope_type);
	}
	trace(q, dq, c, qxm, qxp, dtdx, H.nxt, H.scheme, H.nvar, H.nxyt);
	
	qleftright(idim, H.nx, H.ny, H.nxyt, H.nvar, qxm, qxp, qleft, qright);
	
	// Solve Riemann problem at interfaces
	riemann(qleft, qright, qgdnv,
		rl, ul, pl, cl, wl, rr, ur, pr, cr, wr, ro, uo, po, co, wo,
		rstar, ustar, pstar, cstar,
		sgnm, spin, spout, ushock, frac,
		scr, delp, pold, ind, ind2,
		H.nx + 1, H.smallr, H.smallc, H.gamma, H.niter_riemann,
		H.nvar, H.nxyt);

	// Compute fluxes
	cmpflx(qgdnv, flux, H.nxt, H.nxyt, H.nvar, H.gamma);
	updateConservativeVars(idim, j, dtdx, uold, u, flux, H.imin,
			       H.imax, H.jmin, H.jmax, H.nvar, H.nxt, H.nyt, H.nxyt);
      }                       // end for j

#pragma omp single
{
      if (H.prt) {
	printf("After pass %ld\n", idim);
	PRINTUOLD(H, Hv);
      }
} // End omp single
    } else {
#pragma omp for schedule(dynamic,2)
      for (i = H.imin + ExtraLayer; i < H.imax - ExtraLayer; i++) {
	qID = &Hvw->q[IHvw(0, ID)];
	qIP = &Hvw->q[IHvw(0, IP)];
	gatherConservativeVars(idim, i, uold, u, H.imin, H.imax, H.jmin,
			       H.jmax, H.nvar, H.nxt, H.nyt, H.nxyt);
	PRINTARRAYV(Hvw->u, H.nyt, "uY", H);
	
	// Convert to primitive variables
	constoprim(u, q, e, H.nyt, H.nxyt, H.nvar, H.smallr);
	
	equation_of_state(qID, e, qIP, c, 0, H.nyt, H.smallc, H.gamma);
	PRINTARRAY(Hw->c, H.nyt, "cY", H);
	
	// Characteristic tracing
	// compute slopes
	Dmemset(dq, 0, H.nyt * H.nvar);
	if (H.iorder != 1) {
	  slope(q, dq, H.nyt, H.nvar, H.nxyt, H.slope_type);
	}
	PRINTARRAYV(Hvw->dq, H.nyt, "dqY", H);
	trace(q, dq, c, qxm, qxp, dtdx, H.nyt, H.scheme, H.nvar, H.nxyt);
	qleftright(idim, H.nx, H.ny, H.nxyt, H.nvar, qxm, qxp, qleft, qright);
	PRINTARRAYV(Hvw->qleft, H.ny + 1, "qleftY", H);
	PRINTARRAYV(Hvw->qright, H.ny + 1, "qrightY", H);
	
	// Solve Riemann problem at interfaces
	riemann(qleft, qright, qgdnv, rl, ul,
		pl, cl, wl, rr, ur, pr,
		cr, wr, ro, uo, po, co,
		wo, rstar, ustar, pstar, cstar,
		sgnm, spin, spout, ushock, frac,
		scr, delp, pold, ind, ind2,
		H.ny + 1, H.smallr, H.smallc, H.gamma, H.niter_riemann,
		H.nvar, H.nxyt);

	// Compute fluxes
	cmpflx(qgdnv, flux, H.nyt, H.nxyt, H.nvar, H.gamma);
	PRINTARRAYV(Hvw->flux, H.ny + 1, "fluxY", H);
	// updateConservativeVars(idim, i, dtdx, H, Hv, Hvw);
	updateConservativeVars(idim, i, dtdx, uold, u, flux, H.imin,
			       H.imax, H.jmin, H.jmax, H.nvar, H.nxt, H.nyt, H.nxyt);
      }                       // end for i
#pragma omp single
{
      if (H.prt) {
	  printf("After pass %ld\n", idim);
	  PRINTUOLD(H, Hv);
      }
} // End omp single
    }
    // Deallocate work space
    deallocate_work_space(H, Hw, Hvw);
}                               // hydro_godunov


// EOF
