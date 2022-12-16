/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * gkmetis.c
 *
 * This is the entry point of parallel geometry based partitioning
 * routines
 *
 * Started 10/19/96
 * George
 *
 * $Id: gkmetis.c,v 1.8 2003/07/31 16:23:30 karypis Exp $
 *
 */

#include <parmetislib.h>




/***********************************************************************************
* This function is the entry point of the parallel kmetis algorithm that uses
* coordinates to compute an initial graph distribution.
************************************************************************************/
void ParMETIS_V3_PartGeomKway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy,
              idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *ndims, 
	      float *xyz, int *ncon, int *nparts, float *tpwgts, float *ubvec, 
	      int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  int h, i, j;
  int nvtxs = -1, npes, mype;
  int uwgtflag, cut, gcut, maxnvtxs;
  int ltvwgts[MAXNCON];
  int moptions[10];
  CtrlType ctrl;
  idxtype *uvwgt;
  WorkSpaceType wspace;
  GraphType *graph, *mgraph;
  float avg, maximb, balance, *mytpwgts;
  int seed, dbglvl = 0;
  int iwgtflag, inumflag, incon, inparts, ioptions[10];
  float *itpwgts, iubvec[MAXNCON];

  MPI_Comm_size(*comm, &npes);
  MPI_Comm_rank(*comm, &mype);

  /* If too many processors switch to non-geometric partitioning.
     This is to take care the current poor implementation of sorting
     that has an npes*npes memory complexity. The following fix assumes
     that the machine can allocate about 128MB of memory per node just 
     for sorting alone. 
     Also, if each processor does not have npes vertices, switch to the 
     non-geometric version of the code.
  */
  ctrl.comm = *comm;
  if (npes > 4096 || GlobalSEMin(&ctrl, vtxdist[mype+1]-vtxdist[mype]) < npes) {
    return ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, 
               numflag, ncon, nparts, tpwgts, ubvec, options, edgecut, part, comm);
  }


  /********************************/
  /* Try and take care bad inputs */
  /********************************/
  if (options != NULL && options[0] == 1)
    dbglvl = options[PMV3_OPTION_DBGLVL];

  CheckInputs(STATIC_PARTITION, npes, dbglvl, wgtflag, &iwgtflag, numflag, &inumflag,
              ncon, &incon, nparts, &inparts, tpwgts, &itpwgts, ubvec, iubvec, 
	      NULL, NULL, options, ioptions, part, comm);


  /*********************************/
  /* Take care the nparts = 1 case */
  /*********************************/
  if (inparts <= 1) {
    idxset(vtxdist[mype+1]-vtxdist[mype], 0, part);
    *edgecut = 0;
    return;
  }

  /******************************/
  /* Take care of npes = 1 case */
  /******************************/
  if (npes == 1 && inparts > 1) {
    moptions[0] = 0;
    nvtxs = vtxdist[1];

    if (incon == 1) {
      METIS_WPartGraphKway(&nvtxs, xadj, adjncy, vwgt, adjwgt, &iwgtflag, &inumflag, 
            &inparts, itpwgts, moptions, edgecut, part);
    }
    else {
      /* ADD: this is because METIS does not support tpwgts for all constraints */
      mytpwgts = fmalloc(inparts, "mytpwgts");
      for (i=0; i<inparts; i++)
        mytpwgts[i] = itpwgts[i*incon];

      moptions[7] = -1;
      METIS_mCPartGraphRecursive2(&nvtxs, &incon, xadj, adjncy, vwgt, adjwgt, &iwgtflag, 
            &inumflag, &inparts, mytpwgts, moptions, edgecut, part);

      GKfree((void **)&mytpwgts, LTERM);
    }

    return;
  }


  if (inumflag == 1)
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  /*****************************/
  /* Set up control structures */
  /*****************************/
  if (ioptions[0] == 1) {
    dbglvl = ioptions[PMV3_OPTION_DBGLVL];
    seed = ioptions[PMV3_OPTION_SEED];
  }
  else {
    dbglvl = GLOBAL_DBGLVL;
    seed = GLOBAL_SEED;
  }
  SetUpCtrl(&ctrl, npes, dbglvl, *comm);
  ctrl.CoarsenTo   = amin(vtxdist[npes]+1, 25*incon*amax(npes, inparts));
  ctrl.seed        = (seed == 0) ? mype : seed*mype;
  ctrl.sync        = GlobalSEMax(&ctrl, seed);
  ctrl.partType    = STATIC_PARTITION;
  ctrl.ps_relation = -1;
  ctrl.tpwgts      = itpwgts;
  scopy(incon, iubvec, ctrl.ubvec);

  uwgtflag = iwgtflag|2;
  uvwgt = idxsmalloc(vtxdist[mype+1]-vtxdist[mype], 1, "uvwgt");
  graph = Mc_SetUpGraph(&ctrl, 1, vtxdist, xadj, uvwgt, adjncy, adjwgt, &uwgtflag);
  GKfree((void **)&graph->nvwgt, &uvwgt, LTERM); 

  AllocateWSpace(&ctrl, graph, &wspace);

  /*=================================================================
   * Compute the initial npes-way partitioning geometric partitioning
   =================================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  Coordinate_Partition(&ctrl, graph, *ndims, xyz, 1, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));

  /*=================================================================
   * Move the graph according to the partitioning
   =================================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.MoveTmr));

  graph->vwgt = ((iwgtflag&2) != 0) ? vwgt : idxsmalloc(graph->nvtxs*incon, 1, "vwgt");
  graph->ncon = incon;
  j = ctrl.nparts;
  ctrl.nparts = ctrl.npes;
  mgraph = Mc_MoveGraph(&ctrl, graph, &wspace);
  ctrl.nparts = j;

  /**********************************************************/
  /* Do the same functionality as Mc_SetUpGraph for mgraph */
  /**********************************************************/
  /* compute tvwgts */
  for (j=0; j<incon; j++)
    ltvwgts[j] = 0;

  for (i=0; i<graph->nvtxs; i++)
    for (j=0; j<incon; j++)
      ltvwgts[j] += mgraph->vwgt[i*incon+j];

  for (j=0; j<incon; j++)
    ctrl.tvwgts[j] = GlobalSESum(&ctrl, ltvwgts[j]);

  /* check for zero wgt constraints */
  for (i=0; i<incon; i++) {
    /* ADD: take care of the case in which tvwgts is zero */
    if (ctrl.tvwgts[i] == 0) {
      if (ctrl.mype == 0) printf("ERROR: sum weight for constraint %d is zero\n", i);
      MPI_Finalize();
      exit(-1);
    }
  }

  /* compute nvwgt */
  mgraph->nvwgt = fmalloc(mgraph->nvtxs*incon, "mgraph->nvwgt");
  for (i=0; i<mgraph->nvtxs; i++)
    for (j=0; j<incon; j++)
      mgraph->nvwgt[i*incon+j] = (float)(mgraph->vwgt[i*incon+j]) / (float)(ctrl.tvwgts[j]);


  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.MoveTmr));

  if (ctrl.dbglvl&DBG_INFO) {
    cut = 0;
    for (i=0; i<graph->nvtxs; i++)
      for (j=graph->xadj[i]; j<graph->xadj[i+1]; j++)
        if (graph->where[i] != graph->where[graph->adjncy[j]])
          cut += graph->adjwgt[j];
    gcut = GlobalSESum(&ctrl, cut)/2;
    maxnvtxs = GlobalSEMax(&ctrl, mgraph->nvtxs);
    balance = (float)(maxnvtxs)/((float)(graph->gnvtxs)/(float)(npes));
    rprintf(&ctrl, "XYZ Cut: %6d \tBalance: %6.3f [%d %d %d]\n",
      gcut, balance, maxnvtxs, graph->gnvtxs, npes);
  }

  /*=================================================================
   * Set up the newly moved graph
   =================================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  ctrl.nparts = inparts;
  AdjustWSpace(&ctrl, mgraph, &wspace);

  /*=======================================================
   * Now compute the partition of the moved graph
   =======================================================*/
  if (vtxdist[npes] < SMALLGRAPH || 
      vtxdist[npes] < npes*20 || 
      GlobalSESum(&ctrl, mgraph->nedges) == 0) {
    IFSET(ctrl.dbglvl, DBG_INFO, rprintf(&ctrl, "Partitioning a graph of size %d serially\n", vtxdist[npes]));
    PartitionSmallGraph(&ctrl, mgraph, &wspace);
  }
  else {
    Mc_Global_Partition(&ctrl, mgraph, &wspace);
  }
  ParallelReMapGraph(&ctrl, mgraph, &wspace);

  /* Invert the ordering back to the original graph */
  ctrl.nparts = npes;
  ProjectInfoBack(&ctrl, graph, part, mgraph->where, &wspace);

  *edgecut = mgraph->mincut;

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));

  /*******************/
  /* Print out stats */
  /*******************/
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));

  if (ctrl.dbglvl&DBG_INFO) {
    rprintf(&ctrl, "Final %d-way CUT: %6d \tBalance: ", inparts, mgraph->mincut);
    avg = 0.0;
    for (h=0; h<incon; h++) {
      maximb = 0.0;
      for (i=0; i<inparts; i++)
        maximb = amax(maximb, mgraph->gnpwgts[i*incon+h]/itpwgts[i*incon+h]);
      avg += maximb;
      rprintf(&ctrl, "%.3f ", maximb);
    }
    rprintf(&ctrl, "  avg: %.3f\n", avg/(float)incon);
  }

  GKfree((void **)&itpwgts, LTERM);
  FreeGraph(mgraph);
  FreeInitialGraphAndRemap(graph, iwgtflag, 1);
  FreeWSpace(&wspace);
  FreeCtrl(&ctrl);

  if (inumflag == 1)
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);

}



/***********************************************************************************
* This function is the entry point of the parallel ordering algorithm.
* This function assumes that the graph is already nice partitioned among the
* processors and then proceeds to perform recursive bisection.
************************************************************************************/
void ParMETIS_V3_PartGeom(idxtype *vtxdist, int *ndims, float *xyz, idxtype *part, MPI_Comm *comm)
{
  int i, npes, mype, nvtxs, firstvtx, dbglvl;
  idxtype *xadj, *adjncy;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph;
  int zeroflg = 0;

  MPI_Comm_size(*comm, &npes);
  MPI_Comm_rank(*comm, &mype);

  if (npes == 1) {
    idxset(vtxdist[mype+1]-vtxdist[mype], 0, part);
    return;
  }

  /* Return without computing a partitioning under the following cases: 
     - The number of processors is greater than 4096 (due to npes^2
       memory complexity of the sorting algorithm implemented).
     - When each processor does not have at least npes elements.
     These retrictions will be fixed in 4.0.
  */
  ctrl.comm = *comm;
  if (npes > 4096 || GlobalSEMin(&ctrl, vtxdist[mype+1]-vtxdist[mype]) < npes) {
    if (mype == 1)
      printf("ParMETIS_V3_PartGeom can only be used for less than 4096 processors "
             "and when each processor has at least npes elements.\n");
    return;
  }


  /* Setup a fake graph to allow the rest of the code to work unchanged */
  dbglvl = 0;

  nvtxs = vtxdist[mype+1]-vtxdist[mype];
  firstvtx = vtxdist[mype];
  xadj = idxmalloc(nvtxs+1, "ParMETIS_PartGeom: xadj");
  adjncy = idxmalloc(nvtxs, "ParMETIS_PartGeom: adjncy");
  for (i=0; i<nvtxs; i++) {
    xadj[i] = i;
    adjncy[i] = firstvtx + (i+1)%nvtxs;
  }
  xadj[nvtxs] = nvtxs;

  /* Proceed with the rest of the code */
  SetUpCtrl(&ctrl, npes, dbglvl, *comm);
  ctrl.seed      = mype;
  ctrl.CoarsenTo = amin(vtxdist[npes]+1, 25*npes);

  graph = Mc_SetUpGraph(&ctrl, 1, vtxdist, xadj, NULL, adjncy, NULL, &zeroflg);

  AllocateWSpace(&ctrl, graph, &wspace);

  /*=======================================================
   * Compute the initial geometric partitioning
   =======================================================*/
  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  Coordinate_Partition(&ctrl, graph, *ndims, xyz, 0, &wspace);

  idxcopy(graph->nvtxs, graph->where, part);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));

  FreeInitialGraphAndRemap(graph, 0, 1);
  FreeWSpace(&wspace);
  FreeCtrl(&ctrl);

  GKfree((void **)&xadj, (void **)&adjncy, LTERM);
}




