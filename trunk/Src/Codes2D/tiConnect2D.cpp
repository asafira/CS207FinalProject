// tiConnect2D.m
// function [EToE,EToF]= tiConnect2D(EToV)
// 2007/08/02
//---------------------------------------------------------
#include "NDGLib_headers.h"

//---------------------------------------------------------
void tiConnect2D
(
  IMat& EToV,   // [in]
  IMat& EToE,   // [out]
  IMat& EToF    // [out]
)
//---------------------------------------------------------
{
  // function [EToE,EToF]= tiConnect2D(EToV)
  // Purpose: triangle face connect algorithm due to Toby Isaac

  int Nfaces=3, K = EToV.num_rows(), Nnodes = EToV.max_val();
  int Nr = 3*K;  IMat spNodeToNode, sorted;
  
  IMat fnodes(Nr,2,"fnodes");
  Index1D I1(1,K), I2(K+1,2*K), I3(2*K+1,3*K);
  IVec J1(gVecData, 2, "1 2"), J2(gVecData, 2, "2 3"), J3(gVecData, 2, "3 1");

  // create list of all faces 1, then 2, & 3
  fnodes(I1, All) = EToV(All, J1);
  fnodes(I2, All) = EToV(All, J2);
  fnodes(I3, All) = EToV(All, J3);

  // sort each row in ascending order, then decrement each entry
  fnodes.sort(2);  fnodes -= 1;

  // set up default element to element and Element to faces connectivity
  EToE = outer(Range(1,K), Ones(Nfaces));
  EToF = outer(Ones(K), Range(1,Nfaces));

  // uniquely number each set of three faces by their node numbers 
  IVec id = fnodes.get_col(1)*Nnodes + fnodes.get_col(2)+1;

//spNodeToNode=[id, (1:Nfaces*K)', EToE(:), EToF(:)];
  spNodeToNode.resize(Nr, 4);
  spNodeToNode.set_col(1, id);
  spNodeToNode.set_col(2, Range(1,Nfaces*K));
  spNodeToNode.set_col(3, EToE);
  spNodeToNode.set_col(4, EToF);

  // Now we sort by global face number.
  sorted = sortrows(spNodeToNode,1);

  // find matches in the sorted face list
  // [indices,dummy]=find( sorted(1:(end-1),1)==sorted(2:end,1) );

  Index1D II1(1, Nr-1), II2(2, Nr);
  IVec ids1=sorted(II1,1), ids2=sorted(II2,1);
  IVec indices = find(ids1, '=', ids2);
  
  int Ni=indices.size();

  if (Ni>0) {

    II1.reset(1,Ni), II2.reset(Ni+1,2*Ni);
    IMat matchL(2*Ni, 4), matchR(2*Ni, 4);

    // make links reflexive 
    matchL(II1, All) = sorted(indices  ,All);
    matchL(II2, All) = sorted(indices+1,All);

    matchR(II1, All) = sorted(indices+1,All);
    matchR(II2, All) = sorted(indices  ,All);

    // insert matches
    EToE(matchL(All,2)) = matchR(All,3);
    EToF(matchL(All,2)) = matchR(All,4);
  }
}
