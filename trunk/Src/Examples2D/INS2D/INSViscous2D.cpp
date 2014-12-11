// INSViscous2D.m
// compute right hand side for viscous step solves
// by Kevin Chen
// 2014/11/22
//---------------------------------------------------------
#include "NDGLib_headers.h"
#include "INS2D.h"


//---------------------------------------------------------
void INS2D::INSViscous2D()
//---------------------------------------------------------
{
  double t1 = timer.read(), t2,t3;

  ////////////////////commented out by Kevin ////////////////
  /*
  // compute right hand side for viscous step solves
  DMat mmUxTT = (m_cub.VT)*(m_cub.W.dm(m_cub.V*UxTT));
  Uxrhs  = mmUxTT*(g0/(nu*dt)) + rhsbcUx;
  
  DMat mmUyTT = (m_cub.VT)*(m_cub.W.dm(m_cub.V*UyTT));
  Uyrhs  = mmUyTT*(g0/(nu*dt)) + rhsbcUy;
  */
  ///////////////////////////////////////////////////////////

  /////////////////added by Kevin ///////////////////////////
  // compute right hand side for viscous step solves
  DMat mmUxTT = MassMatrix*(J.dm(UxTT));
  Uxrhs = mmUxTT * (g0/(nu*dt)) + rhsbcUx;

  DMat mmUyTT = MassMatrix*(J.dm(UyTT));
  Uyrhs = mmUyTT * (g0/(nu*dt)) + rhsbcUy;
  ///////////////////////////////////////////////////////////
  
  // save Ux,Uy
  Uxold = Ux; Uyold = Uy;

  // viscous solves (Cholesky, CG, LU, GMRES solvers)
  t2 = timer.read();
  Ux = VELsystemC->solve(Uxrhs);
  Uy = VELsystemC->solve(Uyrhs);
  t3 = timer.read();

  //---------------------------
  time_viscous_sol += t3 - t2;
  time_viscous     += t3 - t1;
  //---------------------------
}
