#include "ComMod.h"

#include <iostream>

//--------
// ComMod 
//--------
//
ComMod::ComMod() 
{
  mvMsh        = false;

  stFileFlag   = false;
  stFileRepl   = false;

  bin2VTK      = false;
  saveAve      = false;
  sepOutput    = false;
  saveATS      = 1;
  saveIncr     = 10;
  nITs         = 0;
  startTS      = 0;
  stFileName   = "stFile";
  iniFilePath  = "";
  stopTrigName = "STOP_SIM";
  ichckIEN     = true; 
  zeroAve      = false;
  cmmInit      = false;
  cmmVarWall   = false;
  shlEq        = false;
  pstEq        = false;
  sstEq        = false;
  ibFlag       = false;

}

//---------
// ~ComMod 
//---------
//
ComMod::~ComMod() 
{
}

///////////////////////
//   a d j T y p e   //
///////////////////////

void adjType::destroy()
{
  nnz = 0;

  pcol.clear();

  prow.clear();
}

///////////////////////
// c p l B C T y p e //
///////////////////////

cplBCType::cplBCType()
{
  schm = consts::CplBCType::cplBC_NA;
}

///////////////////
// d m n T y p e //
///////////////////

dmnType::dmnType()
{
}

dmnType::~dmnType()
{
}

/////////////////
// e q T y p e //
/////////////////

eqType::eqType()
{
  roInf = 0.2;
}

eqType::~eqType()
{
}

//////////////////////
// f a c e  T y p e //
//////////////////////

faceType::faceType()
{
}

faceType::~faceType()
{
}

//---------
// destroy
//---------
// Free memory and reset some data members.
//
// This replicates the Fortran 'SUBROUTINE DESTROYFACE(lFa)' 
// implemented in ALLFUN.f. 
//
void faceType::destroy()
{
  gE.clear();     
  gN.clear();    
  lN.clear();   
  IEN.clear();  
  gebc.clear();
  w.clear();  
  x.clear(); 
  xi.clear(); 
  N.clear();      
  nV.clear();    
  Nx.clear();   
  Nxx.clear(); 

  nAdj.destroy();
  eAdj.destroy();

  for (int i = 0; i < nFs; i++) {
    fs[i].destroy();
  } 

  eType = consts::ElementType::NA;
  nEl = 0;
  nNo = 0;
  gnEl= 0;
}

///////////////////
// f s T y p e   //
///////////////////

fsType::fsType()
{
}

//---------
// destroy
//---------
//
// SUBROUTINE DESTROYFS(fs)
//
void fsType::destroy()
{
  eType = consts::ElementType::NA;

  w.clear();
  xi.clear();
  xib.clear();
  N.clear();
  Nb.clear();
  Nx.clear();
  Nxx.clear();
}

///////////////////
// m s h T y p e //
///////////////////

mshType::mshType()
{
  eType = consts::ElementType::NA;
}

/////////////////////
// r m s h T y p e //
/////////////////////

rmshType::rmshType()
{
  isReqd  = false;
}


