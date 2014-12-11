// PoissonIPDG2D.m
// function [OP,MM] = PoissonIPDG2D(gauss, cub)
// 2014/11/22
// by Kevin Chen
//---------------------------------------------------------
#include "NDGLib_headers.h"
#include "NDG2D.h"


//---------------------------------------------------------
void NDG2D::PoissonIPDG2D
(
  CSd&      spOP,  // [out] sparse
  CSd&      spMM   // [out] sparse
)
//---------------------------------------------------------
{
  // function [OP,MM] = PoissonIPDG2D(gauss, cub)
  //
  // Purpose: Set up the discrete Poisson matrix directly
  //          using LDG. The operator is set up in the weak form
  DVec faceR("faceR"), faceS("faceS");
  DMat V1D;  IVec Fm("Fm");  IVec i1_Nfp = Range(1,Nfp);
  double opti1=0.0, opti2=0.0; int i=0; 

  umLOG(1, "\n ==> {OP,MM} assembly: ");
  opti1 = timer.read(); // time assembly

  // build local face matrices
  DMat massEdge[4]; // = zeros(Np,Np,Nfaces);
  for (i=1; i<=Nfaces; ++i) {
    massEdge[i].resize(Np,Np);
  }

  // face mass matrix 1
  Fm = Fmask(All,1); faceR = r(Fm); 
  V1D = Vandermonde1D(N, faceR);
  massEdge[1](Fm,Fm) = inv(V1D*trans(V1D));

  // face mass matrix 2
  Fm = Fmask(All,2); faceR = r(Fm); 
  V1D = Vandermonde1D(N, faceR);
  massEdge[2](Fm,Fm) = inv(V1D*trans(V1D));
  
  // face mass matrix 3
  Fm = Fmask(All,3); faceS = s(Fm); 
  V1D = Vandermonde1D(N, faceS); 
  massEdge[3](Fm,Fm) = inv(V1D*trans(V1D));

  // build local volume mass matrix
  MassMatrix = trans(invV)*invV;

  DMat Dx("Dx"),Dy("Dy"),Dx2("Dx2"),Dy2("Dy2");
  DMat Dn1("Dn1"),Dn2("Dn2"), mmE("mmE"), OP11("OP11"), OP12("OP12");
  DMat mmE_All_Fm1, mmE_Fm1_Fm1, Dn2_Fm2_All;
  IMat rows1,cols1,rows2,cols2;  int k1=0,f1=0,k2=0,f2=0,id=0;
  Index1D entries, entriesMM, idsM;  IVec fidM,vidM,Fm1,vidP,Fm2;
  double lnx=0.0,lny=0.0,lsJ=0.0,hinv=0.0,gtau=0.0;
  double N1N1 = double((N+1)*(N+1)); int NpNp = Np*Np;

  // build DG derivative matrices
  int max_OP = (K*Np*Np*(1+Nfaces));
  int max_MM = (K*Np*Np);

  // "OP" triplets (i,j,x), extracted to {Ai,Aj,Ax}
  IVec OPi(max_OP), OPj(max_OP), Ai,Aj; DVec OPx(max_OP), Ax;
  // "MM" triplets (i,j,x)
  IVec MMi(max_MM), MMj(max_MM); DVec MMx(max_MM);
  IVec OnesNp = Ones(Np);

  // global node numbering
  entries.reset(1,NpNp); entriesMM.reset(1,NpNp);

  OP12.resize(Np,Np);
  
  for (k1=1; k1<=K; ++k1)
  {
    if (! (k1%250)) { umLOG(1, "%d, ",k1); }

    rows1 = outer( Range((k1-1)*Np+1,k1*Np), OnesNp );
    cols1 = trans(rows1);

    // Build local operators  
    Dx = rx(1,k1)*Dr + sx(1,k1)*Ds;  
    Dy = ry(1,k1)*Dr + sy(1,k1)*Ds;

    OP11 = J(1,k1)*(trans(Dx)*MassMatrix*Dx + 
                    trans(Dy)*MassMatrix*Dy );

    // Build element-to-element parts of operator
    for (f1=1; f1<=Nfaces; ++f1) {
      k2 = EToE(k1,f1); f2 = EToF(k1,f1); 

      rows2 = outer( Range((k2-1)*Np+1, k2*Np), OnesNp );
      cols2 = trans(rows2);

      fidM  = (k1-1)*Nfp*Nfaces + (f1-1)*Nfp + i1_Nfp;
      vidM = vmapM(fidM); Fm1 = mod(vidM-1,Np)+1;
      vidP = vmapP(fidM); Fm2 = mod(vidP-1,Np)+1;

      id = 1+(f1-1)*Nfp + (k1-1)*Nfp*Nfaces;
      lnx = nx(id);  lny = ny(id); lsJ = sJ(id); 
      hinv = std::max(Fscale(id), Fscale(1+(f2-1)*Nfp, k2));    

      Dx2 = rx(1,k2)*Dr + sx(1,k2)*Ds;   
      Dy2 = ry(1,k2)*Dr + sy(1,k2)*Ds;
      
      Dn1 = lnx*Dx  + lny*Dy;
      Dn2 = lnx*Dx2 + lny*Dy2;

      mmE = lsJ*massEdge[f1];

      gtau = 10.0 * N1N1 * hinv; // set penalty scaling

      switch (BCType(k1,f1)){
      case BC_Dirichlet:
        OP11 += ( gtau*mmE - mmE*Dn1 - trans(Dn1)*mmE ); // ok
        break;
      case BC_Neuman:
        //do nothing
        break;
      default:
         // interior face variational terms
        OP11 += 0.5*( gtau*mmE - mmE*Dn1 - trans(Dn1)*mmE );

        // extract mapped regions:
        mmE_All_Fm1 = mmE(All,Fm1);
        mmE_Fm1_Fm1 = mmE(Fm1,Fm1);
        Dn2_Fm2_All = Dn2(Fm2,All);

        OP12 = 0.0;   // reset to zero
        OP12(All,Fm2)  = -0.5*(       gtau*mmE_All_Fm1 );
        OP12(Fm1,All) -=  0.5*(            mmE_Fm1_Fm1*Dn2_Fm2_All );
        //OP12(All,Fm2) -=  0.5*(-trans(Dn1)*mmE_All_Fm1 );
        OP12(All,Fm2) +=  0.5*( trans(Dn1)*mmE_All_Fm1 );

        // load this set of triplets
        OPi(entries)=rows1; OPj(entries)=cols2, OPx(entries)=OP12;
        entries += (NpNp);
      }
    }

      
    OPi(entries  )=rows1; OPj(entries  )=cols1, OPx(entries  )=OP11;
    MMi(entriesMM)=rows1; MMj(entriesMM)=cols1; MMx(entriesMM)=J(1,k1)*MassMatrix;
    entries += (NpNp); entriesMM += (NpNp);
  }
  umLOG(1, "\n ==> {OP,MM} to sparse\n");

  entries.reset(1, entries.hi()-Np*Np);

  // Extract triplets from the large buffers. Note: this 
  // requires copying each array, and since these arrays 
  // can be HUGE(!), we force immediate deallocation:

  Ai=OPi(entries);  OPi.Free();
  Aj=OPj(entries);  OPj.Free();
  Ax=OPx(entries);  OPx.Free();
  umLOG(1, " ==> triplets ready (OP) nnz = %10d\n", entries.hi());

  // adjust triplet indices for 0-based sparse operators
  Ai -= 1; Aj -= 1; MMi -= 1; MMj -= 1;  int npk=Np*K;

//#if defined(NDG_USE_CHOLMOD) || defined(NDG_New_CHOLINC)
  // load only the lower triangle tril(OP)        free args?
  //spOP.load(npk,npk, Ai,Aj,Ax, sp_LT, false,1e-15, true);  // {LT, false} -> TriL
//#else
  // select {upper,lower,both} triangles
  spOP.load(npk,npk, Ai,Aj,Ax, sp_LT, true,1e-15,true);   // LT -> enforce symmetry
//spOP.load(npk,npk, Ai,Aj,Ax, sp_All,true,1e-15,true);   // All-> includes "noise"
//spOP.load(npk,npk, Ai,Aj,Ax, sp_UT, false,1e-15,true);  // UT -> triu(OP) only
//#endif
  
  Ai.Free();  Aj.Free();  Ax.Free();

  umLOG(1, " ==> triplets ready (MM) nnz = %10d\n", entriesMM.hi());

  //-------------------------------------------------------
  // The mass matrix operator will NOT be factorised, 
  // Load ALL elements (both upper and lower triangles):
  //-------------------------------------------------------
  spMM.load(npk,npk, MMi,MMj,MMx, sp_All,false,1.00e-15,true);
  MMi.Free(); MMj.Free(); MMx.Free();

  opti2 = timer.read(); // time assembly
  umLOG(1, " ==> {OP,MM} converted to csc.  (%g secs)\n", opti2-opti1);
}


//---------------------------------------------------------
void NDG2D::PoissonIPDG2D
(
  Gauss2D&  gauss, // [in]
  Cub2D&    cub,   // [in]
  CSd&      spOP,  // [out] sparse
  CSd&      spMM   // [out] sparse
  ){
}
/*
//---------------------------------------------------------
{
  // function [OP,MM] = PoissonIPDG2D(gauss, cub)
  //
  // Purpose: Set up the discrete Poisson matrix directly
  //          using LDG. The operator is set up in the weak form

  NGauss = gauss.NGauss;

  // build DG derivative matrices
  int max_OP = (K*Np*Np*(1+Nfaces));
  int max_MM = (K*Np*Np);

  // "OP" triplets (i,j,x), extracted to {Ai,Aj,Ax}
  IVec OPi(max_OP), OPj(max_OP), Ai,Aj; DVec OPx(max_OP), Ax;
  // "MM" triplets (i,j,x)
  IVec MMi(max_MM), MMj(max_MM); DVec MMx(max_MM);

  DMat gDxM, gDyM, gDxP, gDyP, gDnM, gDnP, OP11, OP12;
  IVec idsPR;  IMat rows1, cols1, cols2;  Index1D entries, entriesMM, idsM;
  DVec xk1, yk1, xk2, yk2, locmm;
  DMat cDx, cDy, gVM, gVMT, gVP, gVPR;
  DMat_Diag cw, gnx,gny,gw; double opti1=0.0, opti2=0.0;
  int k1=0,f1=0, k2=0,f2=0; double hinv=0.0, gtau=0.0;

  // global node numbering
  entries.reset(1,Np*Np); entriesMM.reset(1,Np*Np);

  umMSG(1, "\n ==> {OP,MM} assembly: ");
  opti1 = timer.read(); // time assembly

  for (k1=1; k1<=K; ++k1)
  {
    if (! (k1%100)) { umMSG(1, "%4d, ",k1); }

    rows1 = outer( Range((k1-1)*Np+1,k1*Np), Ones(Np) );
    cols1 = trans(rows1);

    // Build local operators  
    locmm = cub.mm(All,k1);
    cw  = cub.W(All,k1);

    xk1 = x(All,k1); yk1 = y(All,k1);
    PhysDmatrices2D(xk1, yk1, cub.V, cDx, cDy);
    OP11 = trans(cDx)*cw*cDx + trans(cDy)*cw*cDy;
    
    // Build element-to-element parts of operator
    for (f1=1; f1<=Nfaces; ++f1)
    {
      k2 = EToE(k1,f1); f2 = EToF(k1,f1);

      idsM.reset((f1-1)*NGauss+1, f1*NGauss);
      idsPR.range(NGauss,1);

      gVM = gauss.finterp[f1];
      gVP = gauss.finterp[f2];  gVP.reverse_rows();
      gVMT= trans(gVM);   // store transpose

      xk1 = x(All,k1); yk1 = y(All,k1);
      xk2 = x(All,k2); yk2 = y(All,k2);
      PhysDmatrices2D(xk1,yk1,gVM,  gDxM,gDyM);
      PhysDmatrices2D(xk2,yk2,gVP,  gDxP,gDyP);
      gnx = gauss.nx(idsM, k1);
      gny = gauss.ny(idsM, k1);
      gw  = gauss.W (idsM, k1);

      gDnM = gnx*gDxM + gny*gDyM;
      gDnP = gnx*gDxP + gny*gDyP;

      cols2 = outer( Ones(Np), Range((k2-1)*Np+1,k2*Np));
      hinv = std::max(Fscale(1+(f1-1)*Nfp, k1), Fscale(1+(f2-1)*Nfp, k2));
    //gtau = 100*2*(N+1)*(N+1)*hinv; // set penalty scaling
      gtau = ( 5*2*(N+1)*(N+1))*hinv; // set penalty scaling

      switch (BCType(k1,f1)) {
      case BC_Dirichlet:
        OP11 += ( gVMT*(gw*gtau)*gVM - gVMT*gw*gDnM - trans(gDnM)*gw*gVM );
        break;
      case BC_Neuman:
        // nada 
        break;
      default:
        // interior face variational terms
        OP11 +=  0.5*( gVMT*(gw*gtau)*gVM - gVMT*gw*gDnM - trans(gDnM)*gw*gVM );
        OP12  = -0.5*( gVMT*(gw*gtau)*gVP + gVMT*gw*gDnP - trans(gDnM)*gw*gVP );
        OPi(entries)=rows1; OPj(entries)=cols2; OPx(entries)=OP12;
        entries += (Np*Np);
        break;
      }
    }

    OPi(entries  )=rows1; OPj(entries  )=cols1; OPx(entries  )=OP11;
    MMi(entriesMM)=rows1; MMj(entriesMM)=cols1; MMx(entriesMM)=locmm;
    entries += (Np*Np);  entriesMM += (Np*Np);
  }
  umMSG(1, "\n ==> {OP,MM} to sparse\n");

  entries.reset(1, entries.hi()-Np*Np);

  // Extract triplets from the large buffers. Note: this 
  // requires copying each array, and since these arrays 
  // can be HUGE(!), we force immediate deallocation:

  Ai=OPi(entries);  OPi.Free();
  Aj=OPj(entries);  OPj.Free();
  Ax=OPx(entries);  OPx.Free();

  // return 0-based sparse results
  // Note: create strictly symmetric matrices by loading 
  // the upper triangle, then adding (transpose-diag)
  Ai -= 1; Aj -= 1; MMi -= 1; MMj -= 1;  int npk=Np*K;

#ifdef NDG_USE_CHOLMOD
  // pass only the lower triangle
  spOP.load(npk,npk, Ai,Aj,Ax, sp_LT,false,1e-15,true);  // {LT, false} -> TriL
#else 
  // pass both upper and lower triangles
  spOP.load(npk,npk, Ai,Aj,Ax, sp_LT, true,1e-15,true);  // LT -> enforce symmetry
//spOP.load(npk,npk, Ai,Aj,Ax, sp_All,true,1e-15,true);  // All-> includes "noise"
#endif

  Ai.Free();  Aj.Free();  Ax.Free();

  //-------------------------------------------------------
  // The mass matrix operator will NOT be factorised, 
  // only used to build velocity system 
  //
  //  VELsystem += (*mm) * (g0/(dt*nu));
  //
  // Load ALL elements (both upper and lower triangles):
  //-------------------------------------------------------
  spMM.load(npk,npk, MMi,MMj,MMx, sp_All,false,1e-15,true);
  MMi.Free(); MMj.Free(); MMx.Free();


#if (0)
  //#######################################################
  // NBN: checking symmetry  (A-A')./abs(A)
  //#######################################################
  // CSd A2 = (spOP-trans(spOP)); A2.droptol();
  // A2.scale(abs(spOP));
  // double d_n = spOP.nonsymmetry();
  // umLOG(1, "\n *** max non-symmetry (scaled): %g ***\n\n", d_n);
  // FILE* fp=fopen("errA.dat", "w"); errA.write_ML(fp); fclose(fp);
  // umERROR("Nigel", "checking symmetry");
  //#######################################################
#endif

  opti2 = timer.read(); // time assembly
  umMSG(1, " ==> {OP,MM} ready.  (%g secs)\n", opti2-opti1);

#if (1)
  // check on original estimates for nnx
  umMSG(1, " ==> max_OP: %12d\n", max_OP);
  umMSG(1, " ==> nnz_OP: %12d\n", entries.hi());
  umMSG(1, " ==> max_MM: %12d\n", max_MM);
#endif
}
*/
