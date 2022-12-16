
///////////////////////////////////////////////////////////
//      S e t   M a t e r i a l   P r o p e r t i e s    //
///////////////////////////////////////////////////////////
//
// The 'set_material_props' map defined here sets material
// properties from values read in from a file.

using SeMaterialPropertiesMapType = std::map<consts::ConstitutiveModelType, 
    std::function<void(DomainParameters*, double, double, double, dmnType&)>>;

SeMaterialPropertiesMapType set_material_props = {

//---------------------------//
//       stIso_nHook         //
//---------------------------//
//
{consts::ConstitutiveModelType::stIso_nHook, [](DomainParameters* domain_params, double mu, double kap, double lam, 
    dmnType& lDmn) -> void
{
  lDmn.stM.isoType = consts::ConstitutiveModelType::stIso_nHook;
  lDmn.stM.C10 = 0.5 * mu;
} },

//---------------------------//
//       stIso_lin           //
//---------------------------//
//
{consts::ConstitutiveModelType::stIso_lin, [](DomainParameters* domain_params, double mu, double kap, double lam,
    dmnType& lDmn) -> void
{
  lDmn.stM.isoType = consts::ConstitutiveModelType::stIso_lin;
  lDmn.stM.C10 = mu;
} },

//---------------------------//
//       stIso_StVK          //
//---------------------------//
//
{consts::ConstitutiveModelType::stIso_StVK, [](DomainParameters* domain_params, double mu, double kap, double lam,
    dmnType& lDmn) -> void
{
  lDmn.stM.isoType = consts::ConstitutiveModelType::stIso_StVK;
  lDmn.stM.C10 = lam;
  lDmn.stM.C01 = mu;
  lDmn.stM.Kpen = kap;
} },

//---------------------------//
//       stIso_mStVK         //
//---------------------------//
//
{consts::ConstitutiveModelType::stIso_mStVK, [](DomainParameters* domain_params, double mu, double kap, double lam,
    dmnType& lDmn) -> void
{
  lDmn.stM.isoType = consts::ConstitutiveModelType::stIso_mStVK;
  lDmn.stM.C10 = kap;
  lDmn.stM.C01 = mu;
  lDmn.stM.Kpen = kap;
} },

//---------------------------//
//       stIso_MR            //
//---------------------------//
//
{consts::ConstitutiveModelType::stIso_MR, [](DomainParameters* domain_params, double mu, double kap, double lam,
    dmnType& lDmn) -> void
{
  lDmn.stM.isoType = consts::ConstitutiveModelType::stIso_MR;
  lDmn.stM.C10 = domain_params->constitutive_model.mooney_rivlin.c1.value();
  lDmn.stM.C01 = lDmn.stM.C10 = domain_params->constitutive_model.mooney_rivlin.c2.value();
} },

//---------------------------//
//       stIso_HGO           //
//---------------------------//
//
{consts::ConstitutiveModelType::stIso_HGO, [](DomainParameters* domain_params, double mu, double kap, double lam,
    dmnType& lDmn) -> void
{
  lDmn.stM.isoType = consts::ConstitutiveModelType::stIso_HGO;
  auto& params = domain_params->constitutive_model.holzapfel_gasser_ogden;

  lDmn.stM.C10 = mu * 0.5;
  lDmn.stM.aff = params.a4.value();
  lDmn.stM.bff = params.b4.value();
  lDmn.stM.ass = params.a6.value();
  lDmn.stM.bss = params.b6.value();
  lDmn.stM.kap = params.kappa.value();
} },

//---------------------------//
//       stIso_Gucci         //
//---------------------------//
//
{consts::ConstitutiveModelType::stIso_Gucci, [](DomainParameters* domain_params, double mu, double kap, double lam,
    dmnType& lDmn) -> void
{
  lDmn.stM.isoType = consts::ConstitutiveModelType::stIso_Gucci;
  auto& params = domain_params->constitutive_model.guccione;

  lDmn.stM.C10 = params.c.value();
  lDmn.stM.bff = params.bf.value();
  lDmn.stM.bss = params.bt.value();
  lDmn.stM.bfs =params.bfs.value(); 
} },

//---------------------------//
//       stIso_HO            //
//---------------------------//
//
{consts::ConstitutiveModelType::stIso_HO, [](DomainParameters* domain_params, double mu, double kap, double lam,
    dmnType& lDmn) -> void
{
  lDmn.stM.isoType = consts::ConstitutiveModelType::stIso_HO;
  auto& params = domain_params->constitutive_model.holzapfel;

  lDmn.stM.a = params.a.value();
  lDmn.stM.b = params.b.value();
  lDmn.stM.aff = params.a4f.value(); 
  lDmn.stM.bff = params.b4f.value(); 
  lDmn.stM.ass = params.a4s.value();
  lDmn.stM.bss = params.b4s.value();
  lDmn.stM.afs = params.afs.value();
  lDmn.stM.bfs = params.bfs.value();
} },

};
