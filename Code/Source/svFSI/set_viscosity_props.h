///////////////////////////////////////////////////////////
//    S e t   V i s c o s i t y   P r o p e r t i e s    //
///////////////////////////////////////////////////////////
//
// The 'set_viscosity_props ' map defined here sets fluid
// viscosity properties from values read in from a file.

using SetViscosityPropertiesMapType = std::map<
    consts::FluidViscosityModelType,
    std::function<void(Simulation*, ViscosityParameters&, dmnType& lDmn)>>;

SetViscosityPropertiesMapType set_viscosity_props = {

    //---------------------------//
    //      viscType_Const       //
    //---------------------------//
    //
    {consts::FluidViscosityModelType::viscType_Const,
     [](Simulation* simulation, ViscosityParameters& params,
        dmnType& lDmn) -> void {
       using namespace consts;
       auto& com_mod = simulation->get_com_mod();

       lDmn.visc.viscType = FluidViscosityModelType::viscType_Const;
       lDmn.visc.mu_i = params.newtonian_model.constant_value.value();
     }},

    //---------------------------//
    //      viscType_CY          //
    //---------------------------//
    //
    {consts::FluidViscosityModelType::viscType_CY,
     [](Simulation* simulation, ViscosityParameters& params,
        dmnType& lDmn) -> void {
       using namespace consts;
       auto& com_mod = simulation->get_com_mod();
       auto& model_params = params.carreau_yasuda_model;

       lDmn.visc.viscType = FluidViscosityModelType::viscType_CY;

       lDmn.visc.mu_i = model_params.limiting_high_shear_rate_viscosity.value();
       lDmn.visc.mu_o = model_params.limiting_low_shear_rate_viscosity.value();
       lDmn.visc.lam = model_params.shear_rate_tensor_multipler.value();
       lDmn.visc.a = model_params.shear_rate_tensor_exponent.value();
       lDmn.visc.n = model_params.power_law_index.value();

       if (lDmn.visc.mu_i > lDmn.visc.mu_o) {
         throw std::runtime_error(
             "Unexpected inputs for Carreau-Yasuda model. "
             "High shear-rate viscosity value (" +
             std::to_string(lDmn.visc.mu_i) +
             " should be larger than low shear-rate value " +
             std::to_string(lDmn.visc.mu_i) + ".");
       }
     }},

    //---------------------------//
    //      viscType_Cass        //
    //---------------------------//
    //
    {consts::FluidViscosityModelType::viscType_Cass,
     [](Simulation* simulation, ViscosityParameters& params,
        dmnType& lDmn) -> void {
       using namespace consts;
       auto& com_mod = simulation->get_com_mod();
       auto& model_params = params.cassons_model;

       lDmn.visc.viscType = FluidViscosityModelType::viscType_Cass;
       lDmn.visc.mu_i = model_params.asymptotic_viscosity();
       lDmn.visc.mu_o = model_params.yield_stress();

       if (model_params.low_shear_rate_threshold.defined()) {
         lDmn.visc.lam = model_params.low_shear_rate_threshold();
       } else {
         lDmn.visc.lam = 0.5;
       }
     }},

};
