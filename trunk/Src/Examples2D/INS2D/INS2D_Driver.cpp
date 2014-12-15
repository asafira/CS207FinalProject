// INS2D_Driver.cpp
// Driver for computing the solution of the
// incompressible Navier-Stokes equations
// by Kevin Chen
// 2014/11/22
//---------------------------------------------------------
#include "NDGLib_headers.h"
#include "INS2D.h"


//---------------------------------------------------------
void INS2D::Driver()
//---------------------------------------------------------
{
  umLOG(1, "INS2D::Driver()\n");

  //-------------------------------------
  // Choose simulation type, domain, 
  // initial solution and BC functions
  //-------------------------------------
  sim_type = eVolkerCylinder;

  switch (sim_type) {

  case eVolkerCylinder:
  //FileName        = "Grid/CFD/cylinderA00075b.neu";
  //FileName        = "Grid/CFD/cylinderCA0015.neu";
  //FileName        = "Grid/CFD/Volker_306.neu";
  // FileName        = "Grid/CFD/Volker_374.neu";
    FileName        = "Grid/CFD/cylinderFlowFineForSolver.neu";
    ExactSolution   = &INS2D::INScylinderIC2D;
    ExactSolutionBC = &INS2D::INScylinderBC2D;
    FinalTime = 40.0;  nu = 1e-0;
    break;

  }

  //---------------------------------------------
  // Order of polynomial approximation (N)
  //---------------------------------------------
//N =  5;
  N =  1;    // Note: Volker_374.neu - very coarse at outflow!
//N = 12;

//FinalTime = 1.0;
   FinalTime = 40.0;
//FinalTime = 0.005;

  // Read in Mesh: [vertices, elements, materials, BC's]
  if (!MeshReaderGambit2D(FileName)) {
    umWARNING("INS2D::Driver", "Error loading mesh (file: %s)\nExiting.\n", FileName.c_str());
    return;
  }

  try {
    Run();  // Solve Problem
  } catch (...) {
    umWARNING("INS2D::Driver", "Caught exception from Run()");
    return;
  }
}
