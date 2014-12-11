// PoissonIPDGbc2D.m
// function [OP] = PoissonIPDGbc2D()
// written by Kevin Chen
// 2014/11/22
//---------------------------------------------------------
#include "NDGLib_headers.h"
#include "NDG2D.h"

//---------------------------------------------------------
DVec& NDG2D::PoissonIPDGbc2D
(DVec& ubc, //[in]
 DVec& qbc  //[in]
)
//---------------------------------------------------------
{
  // function [OP] = PoissonIPDGbc2D()
  // Purpose: Set up the discrete Poisson matrix directly
  //          using LDG. The operator is set up in the weak form

  // build DG derivative matrices
  int max_OP = (K*Np*Np*(1+Nfaces));

  // initialize parameters
  DVec faceR("faceR"), faceS("faceS");
  DMat V1D("V1D"), Dx("Dx"),Dy("Dy"), Dn1("Dn1"), mmE_Fm1("mmE(:,Fm1)");
  IVec Fm("Fm"), Fm1("Fm1"), fidM("fidM");
  double lnx=0.0,lny=0.0,lsJ=0.0,hinv=0.0,gtau=0.0;
  int i=0,k1=0,f1=0,id=0;
  IVec i1_Nfp = Range(1,Nfp);
  double N1N1 = double((N+1)*(N+1));

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
 
  // build DG right hand side
  DVec* pBC = new DVec(Np*K, "bc", OBJ_temp); 
  DVec& bc = (*pBC);  // reference, for syntax
  ////////////////////////////////////////////////////////////////

  umMSG(1, "\n ==> {OP} assembly [bc]: ");
  for (k1=1; k1<=K; ++k1)
  {
    if (! (k1%100)) { umMSG(1, "%d, ",k1); }

    // rows1 = outer(Range((k1-1)*Np+1,k1*Np), Ones(NGauss));

    // Build element-to-element parts of operator
    for (f1=1; f1<=Nfaces; ++f1)
    {
      if (BCType(k1,f1))
      {
        ////////////////////////added by Kevin ///////////////////////////////
        Fm1 = Fmask(All,f1); 
        fidM  = (k1-1)*Nfp*Nfaces + (f1-1)*Nfp + i1_Nfp;
        id = 1+(f1-1)*Nfp + (k1-1)*Nfp*Nfaces;

        lnx = nx(id); lny = ny(id); 
        lsJ = sJ(id); hinv = Fscale(id);

        Dx = rx(1,k1)*Dr + sx(1,k1)*Ds;  
        Dy = ry(1,k1)*Dr + sy(1,k1)*Ds;
        Dn1 = lnx*Dx + lny*Dy;

      //mmE = lsJ*massEdge(:,:,f1);
      //bc(All,k1) += (gtau*mmE(All,Fm1) - Dn1'*mmE(All,Fm1))*ubc(fidM);

        mmE_Fm1 = massEdge[f1](All,Fm1);  mmE_Fm1 *= lsJ;

        gtau = 10*N1N1*hinv; // set penalty scaling
        //bc(All,k1) += (gtau*mmE_Fm1 - trans(Dn1)*mmE_Fm1) * ubc(fidM);

        switch(BCType(k1,f1)){
	  case BC_Dirichlet: 
            bc(Np*k1+Range(1,Np)) += (gtau*mmE_Fm1 - trans(Dn1)*mmE_Fm1)*ubc(fidM);  
            break;
          case BC_Neuman:
            bc(Np*k1+Range(1,Np)) += mmE_Fm1*qbc(fidM);
            break;
	default:
	  std::cout<<"warning: boundary condition is incorrect"<<std::endl;
	}
      }
    }
  }
  return bc;
}


void NDG2D::PoissonIPDGbc2D(
  CSd& spOP //[out] sparse operator 
  )
{
  // function [OP] = PoissonIPDGbc2D()
  // Purpose: Set up the discrete Poisson matrix directly
  //          using LDG. The operator is set up in the weak form

  // build DG derivative matrices
  int max_OP = (K*Np*Np*(1+Nfaces));

  //initialize parameters
  DVec faceR("faceR"), faceS("faceS");
  IVec Fm("Fm"), Fm1("Fm1"), fidM("fidM");
  DMat V1D("V1D"); int i=0;

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

  //continue initialize parameters
  DMat Dx("Dx"),Dy("Dy"), Dn1("Dn1"), mmE_Fm1("mmE(:,Fm1)");
  double lnx=0.0,lny=0.0,lsJ=0.0,hinv=0.0,gtau=0.0;
  int k1=0,f1=0,id=0;
  IVec i1_Nfp = Range(1,Nfp);
  double N1N1 = double((N+1)*(N+1));
  
  // "OP" triplets (i,j,x), extracted to {Ai,Aj,Ax}
  IVec OPi(max_OP),OPj(max_OP), Ai,Aj; DVec OPx(max_OP), Ax;
  IMat rows1, cols1;  Index1D entries; DMat OP11(Np,Nfp, 0.0);

  // global node numbering
  entries.reset(1,Np*Nfp); 
  cols1 = outer(Ones(Np), Range(1,Nfp));

  umMSG(1, "\n ==> {OP} assembly [bc]: ");
  for (k1=1; k1<=K; ++k1)
  {
    if (! (k1%100)) { umMSG(1, "%d, ",k1); }
    rows1 = outer(Range((k1-1)*Np+1,k1*Np), Ones(Nfp));

    // Build element-to-element parts of operator
    for (f1=1; f1<=Nfaces; ++f1)
    {
      if (BCType(k1,f1))
      {   
        ////////////////////////added by Kevin ///////////////////////////////
        Fm1 = Fmask(All,f1); 
        fidM  = (k1-1)*Nfp*Nfaces + (f1-1)*Nfp + i1_Nfp;
        id = 1+(f1-1)*Nfp + (k1-1)*Nfp*Nfaces;

        lnx = nx(id); lny = ny(id); 
        lsJ = sJ(id); hinv = Fscale(id);

        Dx = rx(1,k1)*Dr + sx(1,k1)*Ds;  
        Dy = ry(1,k1)*Dr + sy(1,k1)*Ds;
        Dn1 = lnx*Dx + lny*Dy;

      //mmE = lsJ*massEdge(:,:,f1);
      //bc(All,k1) += (gtau*mmE(All,Fm1) - Dn1'*mmE(All,Fm1))*ubc(fidM);

        mmE_Fm1 = massEdge[f1](All,Fm1);  mmE_Fm1 *= lsJ;

        gtau = 10*N1N1*hinv; // set penalty scaling
        //bc(All,k1) += (gtau*mmE_Fm1 - trans(Dn1)*mmE_Fm1) * ubc(fidM);

        switch(BCType(k1,f1)){
	  case BC_Dirichlet: 
            OP11 = gtau*mmE_Fm1 - trans(Dn1)*mmE_Fm1;  
            break;
          case BC_Neuman:
            OP11 = mmE_Fm1;
            break;
	default:
	  std::cout<<"warning: boundary condition is incorrect"<<std::endl;
	}

        OPi(entries)=rows1; OPj(entries)=cols1; OPx(entries)=OP11; 
        entries += (Np*Nfp);
      }
      cols1 += Nfp;
    }
  }

  umMSG(1, "\n ==> {OPbc} to sparse\n");
  entries.reset(1, entries.hi()-(Np*Nfp));

  // extract triplets from large buffers
  Ai=OPi(entries); Aj=OPj(entries); Ax=OPx(entries);

  // These arrays can be HUGE, so force deallocation
  OPi.Free(); OPj.Free(); OPx.Free();

  // return 0-based sparse result
  Ai -= 1; Aj -= 1;

  //-------------------------------------------------------
  // This operator is not symmetric, and will NOT be 
  // factorised, only used to create reference RHS's:
  //
  //    refrhsbcPR = spOP1 * bcPR;
  //    refrhsbcUx = spOP2 * bcUx;
  //    refrhsbcUy = spOP2 * bcUy;
  //
  // Load ALL elements (both upper and lower triangles):
  //-------------------------------------------------------
  spOP.load(Np*K, Nfp*Nfaces*K, Ai,Aj,Ax, sp_All,false, 1e-15,true);

  Ai.Free();  Aj.Free();  Ax.Free();
  umMSG(1, " ==> {OPbc} ready.\n");

#if (1)
  // check on original estimates for nnx
  umMSG(1, " ==> max_OP: %12d\n", max_OP);
  umMSG(1, " ==> nnz_OP: %12d\n", entries.hi());
#endif
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
////////////////////// TO BE DELETED /////////////////////////////////////////////
void NDG2D::PoissonIPDGbc2D
(
 Gauss2D&  gauss,  // [in]
 CSd&      spOP    // [out] sparse
  ){
}


/*
void NDG2D::PoissonIPDGbc2D
(
 Gauss2D&  gauss,  // [in]
 CSd&      spOP    // [out] sparse
)
//---------------------------------------------------------
{
  // function [OP] = PoissonIPDGbc2D()
  // Purpose: Set up the discrete Poisson matrix directly
  //          using LDG. The operator is set up in the weak form
  NGauss = gauss.NGauss;

  // build DG derivative matrices
  int max_OP = (K*Np*Np*(1+Nfaces));
  
#if (1)
  //#####################################
  // NBN: how many boundary nodes?
  //#####################################
  max_OP /= 10;
#endif


  // "OP" triplets (i,j,x), extracted to {Ai,Aj,Ax}
  IVec OPi(max_OP),OPj(max_OP), Ai,Aj; DVec OPx(max_OP), Ax;
  IMat rows1, cols1;  Index1D entries, idsM;
  DMat VM, dVdxM, dVdyM, DnM; DVec xk1, yk1;
  DMat_Diag gnx, gny, gw;  DMat OP11(Np,Nfp, 0.0);
  int k1=0,f1=0; double hinv=0.0, gtau=0.0;

  // global node numbering
  entries.reset(1,Np*NGauss); 
  cols1 = outer(Ones(Np), Range(1,NGauss));


  umMSG(1, "\n ==> {OP} assembly [bc]: ");
  for (k1=1; k1<=K; ++k1)
  {
    if (! (k1%100)) { umMSG(1, "%d, ",k1); }

    rows1 = outer(Range((k1-1)*Np+1,k1*Np), Ones(NGauss));

    // Build element-to-element parts of operator
    for (f1=1; f1<=Nfaces; ++f1)
    {
      if (BCType(k1,f1))
      {
       
        idsM.reset((f1-1)*NGauss+1,f1*NGauss);
        VM = gauss.finterp[f1];

        xk1 = x(All,k1); yk1 = y(All,k1);
        PhysDmatrices2D(xk1,yk1,VM, dVdxM,dVdyM);
        gnx = gauss.nx(idsM, k1);
        gny = gauss.ny(idsM, k1);
        gw  = gauss.W (idsM, k1);
        
        DnM = gnx*dVdxM + gny*dVdyM;
        hinv = Fscale(1+(f1-1)*Nfp, k1);
      //gtau = (100*2*(N+1)*(N+1))*hinv; // set penalty scaling
        gtau = (  5*2*(N+1)*(N+1))*hinv; // set penalty scaling

        switch (BCType(k1,f1)) {
        case BC_Dirichlet: OP11 = ( trans(VM)*gw*gtau - trans(DnM)*gw ); break;
        case BC_Neuman:    OP11 = ( trans(VM)*gw ); break; // recent sign change
        default:           OP11 = 0.0;  // OP11=zeros(Np,Nfp);

        }

        OPi(entries)=rows1; OPj(entries)=cols1; OPx(entries)=OP11; 
        entries += (Np*NGauss);
      }
      cols1 += NGauss;
    }
  }


  umMSG(1, "\n ==> {OPbc} to sparse\n");
  entries.reset(1, entries.hi()-(Np*NGauss));

  // extract triplets from large buffers
  Ai=OPi(entries); Aj=OPj(entries); Ax=OPx(entries);

  // These arrays can be HUGE, so force deallocation
  OPi.Free(); OPj.Free(); OPx.Free();

  // return 0-based sparse result
  Ai -= 1; Aj -= 1;

  //-------------------------------------------------------
  // This operator is not symmetric, and will NOT be 
  // factorised, only used to create reference RHS's:
  //
  //    refrhsbcPR = spOP1 * bcPR;
  //    refrhsbcUx = spOP2 * bcUx;
  //    refrhsbcUy = spOP2 * bcUy;
  //
  // Load ALL elements (both upper and lower triangles):
  //-------------------------------------------------------
  spOP.load(Np*K, NGauss*Nfaces*K, Ai,Aj,Ax, sp_All,false, 1e-15,true);

  Ai.Free();  Aj.Free();  Ax.Free();
  umMSG(1, " ==> {OPbc} ready.\n");

#if (1)
  // check on original estimates for nnx
  umMSG(1, " ==> max_OP: %12d\n", max_OP);
  umMSG(1, " ==> nnz_OP: %12d\n", entries.hi());
#endif

}
*/
