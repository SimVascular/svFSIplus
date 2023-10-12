/*
 This code implements the Darcy Equation for 2D and 3D
 problems in perfusion of porus media.
 -------------------------------------------------------------
 Assumptions:
    - Homogeneous Permeability
    - Homogeneous Density
    - Isotropic Permeability
    - Assumptions of Stokes Flow
    - Steady-State
 -------------------------------------------------------------
 Strong form of the Single-Compartment Darcy equation:
    u = -K∇(P)
    ∇⋅u = β0(P_source - P) - β1(P - P_sink)
 where:
    u -> Volume flux vector
    K -> Permeability tensor
    P -> Pressure
 -------------------------------------------------------------
 Weak form of the Single-Compartment Darcy equation:
    -∫(∇q∇P)dΩ - λ∫qPdΩ = ∫qFdΩ - ∫q∇P⋅nvdΓ
 where:
    q -> Test function
    λ -> (β0 + β1)/K
    F -> -(β0(P_source) + β1(P_sink))/K
    n -> Normal vector to the boundary

*/

#include "darcy.h"

#include "all_fun.h"
#include "lhsa.h"
#include "mat_fun.h"
#include "nn.h"
#include "utils.h"

#ifdef WITH_TRILINOS
#include "trilinos_linear_solver.h"
#endif

namespace darcy {
    //---------
    // b_darcy
    //---------
    //
    void b_darcy(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const double h, Array<double>& lR)
    {
        for (int a = 0; a < eNoN; a++) {
            lR(0,a) = lR(0,a) + w * N(a) * h;
        }
    }

    void construct_darcy(ComMod& com_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Yg) {
#define n_debug_construct_darcy
#ifdef debug_construct_darcy
        DebugMsg dmsg(__func__, com_mod.cm.idcm());
        dmsg.banner();
#endif

        using namespace consts;

        const int nsd = com_mod.nsd;
        const int tDof = com_mod.tDof;
        const int dof = com_mod.dof;
        const int cEq = com_mod.cEq;
        const auto &eq = com_mod.eq[cEq];
        auto &cDmn = com_mod.cDmn;


        int eNoN = lM.eNoN;
#ifdef debug_construct_darcy
        dmsg << "cEq: " << cEq;
        dmsg << "cDmn: " << cDmn;
#endif

        Vector<int> ptr(eNoN);
        Vector<double> N(eNoN);
        Array<double> xl(nsd, eNoN), al(tDof, eNoN), yl(tDof, eNoN);
        Array<double> Nx(nsd, eNoN), lR(dof, eNoN);
        Array3<double> lK(dof * dof, eNoN, eNoN);
        Array<double> ksix(nsd, nsd);
        Vector<double> local_source(eNoN), local_sink(eNoN), local_b0(eNoN), local_b1(eNoN);

        for (int e = 0; e < lM.nEl; e++) {
            cDmn = all_fun::domain(com_mod, lM, cEq, e);
            auto cPhys = eq.dmn[cDmn].phys;
            if (cPhys != EquationType::phys_darcy) {
                std::cout << "cPhys: " << cPhys << std::endl;
                continue;
            }

            // Update shape function for NURBS
            if (lM.eType == ElementType::NRB) {
                //CALL NRMNNX(lm, e)
            }

            // Create local copies
            for (int a = 0; a < eNoN; a++) {
                int Ac = lM.IEN(a, e);
                ptr(a) = Ac;

                for (int i = 0; i < nsd; i++) {
                    xl(i, a) = com_mod.x(i, Ac);
                }

                for (int i = 0; i < tDof; i++) {
                    al(i, a) = Ag(i, Ac);
                    yl(i, a) = Yg(i, Ac);
                }

                local_source(a) = com_mod.perfusion_pressure_source(Ac);
                local_sink(a) = com_mod.perfusion_pressure_sink(Ac);
                local_b0(a) = com_mod.perfusion_beta0(Ac);
                local_b1(a) = com_mod.perfusion_beta1(Ac);
            }

            // Gauss integration

            lR = 0.0;
            lK = 0.0;
            double Jac{0.0};

            for (int g = 0; g < lM.nG; g++) {
                if (g == 0 || !lM.lShpF) {
                    auto Nx_g = lM.Nx.slice(g);
                    nn::gnn(eNoN, nsd, nsd, Nx_g, xl, Nx, Jac, ksix);
                    if (utils::is_zero(Jac)) {
                        throw std::runtime_error(
                                "[construct_darcy] Jacobian for element " + std::to_string(e) + " is < 0.");
                    }
                }

                double w = lM.w(g) * Jac;
                N = lM.N.col(g);

                if (nsd == 3) {
                    darcy_3d(com_mod, eNoN, w, N, Nx, al, yl, lR, lK);
                } else if (nsd == 2) {
                    darcy_2d(com_mod, eNoN, w, N, Nx, al, yl, lR, lK, local_source, local_sink, local_b0, local_b1);
                } else {
                    throw std::runtime_error("[construct_darcy] nsd must be 2 or 3.");
                }
            }

            // Assembly
#ifdef WITH_TRILINOS
            if (eq.assmTLS) {
                trilinos_doassem_(const_cast<int&>(eNoN), const_cast<int*>(ptr.data()), lK.data(), lR.data());
            }
#endif
            lhsa_ns::do_assem(com_mod, eNoN, ptr, lK, lR);
        }
    }
    //----------
    // darcy_2d
    //----------
    //
    void darcy_2d(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx,
                  const Array<double>& al, const Array<double>& yl, Array<double>& lR, Array3<double>& lK,
                  Vector<double>& local_source, Vector<double>& local_sink, Vector<double>& b0, Vector<double>& b1) {
#define n_debug_darcy_2d
#ifdef debug_darcy_2d
        DebugMsg dmsg(__func__, com_mod.cm.idcm());
        dmsg.banner();
        dmsg << "w: " << w;
#endif
        using namespace consts;
        using namespace mat_fun;

        const int nsd = com_mod.nsd;
        const int cEq = com_mod.cEq;
        auto &eq = com_mod.eq[cEq];
        const int cDmn = com_mod.cDmn;
        auto &dmn = eq.dmn[cDmn];
        const double dt = com_mod.dt;
        const int i = eq.s;

        Array<double> K = dmn.prop.at(PhysicalProperyType::permeability); // should be a tensor
        double source = dmn.prop.at(PhysicalProperyType::source_term); //is this just a double?
        double rho = dmn.prop.at(PhysicalProperyType::solid_density); // need fluid density instead
        //double phi = dmn.prop.at(PhysicalProperyType::porosity);
        //double p0 = dmn.prop.at(PhysicalProperyType::porosity_pressure);
        //double beta = dmn.prop.at(PhysicalProperyType::compressibility);
        //double p1 = dmn.prop.at(PhysicalProperyType::density_pressure);
        //double gamma = dmn.prop.at(PhysicalProperyType::gamma);
        //double B = dmn.prop.at(PhysicalProperyType::B);
        //double mu = dmn.prop.at(PhysicalProperyType::viscosity);

        // need to figure out how these are edited for the darcy case
        double T1 = eq.af * eq.gam * dt;
        double amd = eq.am * rho/ T1;
        double wl = w * T1;

#ifdef debug_darcy_2d
        dmsg;
        dmsg << "nu: " << nu;
        dmsg << "source: " << source;
        dmsg << "rho: " << rho ;
        dmsg << "T1: " << T1;
        dmsg << "i: " << i;
        dmsg << "wl: " << wl;
#endif
        double Pd = -source;
        Vector<double> Px(nsd);

        // need to formulate this for the compressible fluid/porous media with respect to pressure only
        for (int a = 0; a < eNoN; a++) {
            Pd = Pd + N(a)*al(i,a);
            Px(0) = Px(0) + Nx(0,a)*yl(i,a);
            Px(1) = Px(1) + Nx(1,a)*yl(i,a);
        }


        Pd = Pd * rho;

#ifdef debug_darcy_2d
        dmsg;
        dmsg << "Td: " << Td;
        dmsg << "Tx: " << Tx;
#endif

        for (int a = 0; a < eNoN; a++) {
            lR(0,a) = lR(0,a) + w*(N(a)*Pd +
                    (Nx(0,a)*(K(0,0)*Px(0)+K(0,1)*Px(1)) +
                     Nx(1,a)*(K(1,0)*Px(0)+K(1,1)*Px(1))));
            for (int b = 0; b < eNoN; b++) {
                lK(0,a,b) = lK(0,a,b) + wl*(N(a)*N(b)*amd +
                        (Nx(0,a)*(K(0,0)*Nx(0,b) + K(0,1)*Nx(1,b)) +
                         Nx(1,a)*(K(1,0)*Nx(0,b) + K(1,1)*Nx(1,b))));
            }
        }
    }

    //----------
    // darcy_3d
    //----------
    //
    void darcy_3d(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx,
                  const Array<double>& al, const Array<double>& yl, Array<double>& lR, Array3<double>& lK) {
#define n_debug_darcy_3d
#ifdef debug_darcy_3d
        DebugMsg dmsg(__func__, com_mod.cm.idcm());
        dmsg.banner();
        dmsg << "w: " << w;
#endif
        using namespace consts;
        using namespace mat_fun;

        const int nsd = com_mod.nsd;
        const int cEq = com_mod.cEq;
        auto &eq = com_mod.eq[cEq];
        const int cDmn = com_mod.cDmn;
        auto &dmn = eq.dmn[cDmn];
        const double dt = com_mod.dt;
        const int i = eq.s;

        double k = dmn.prop.at(PhysicalProperyType::permeability); // should be a tensor
        double source = dmn.prop.at(PhysicalProperyType::source_term); //is this just a double?
        double rho = dmn.prop.at(PhysicalProperyType::solid_density);

        double T1 = eq.af * eq.gam * dt;
        double amd = eq.am * rho / T1;
        double wl = w * T1;

#ifdef debug_darcy_3d
        dmsg;
        dmsg << "k: " << k;
        dmsg << "source: " << source;
        dmsg << "rho: " << rho ;
        dmsg << "T1: " << T1;
        dmsg << "i: " << i;
        dmsg << "wl: " << wl;
#endif
        double Td = -source;
        Vector<double> Tx(nsd);

        for (int a = 0; a < eNoN; a++){
            Td = Td + N(a) * al(i,a);
            Tx(0) = Tx(0) + Nx(0,a) * yl(i,a);
            Tx(1) = Tx(1) + Nx(1,a) * yl(i,a);
            Tx(2) = Tx(2) + Nx(2,a) * yl(i,a);
        }

        Td = Td * rho;

        for (int a = 0; a < eNoN; a++){
            lR(0,a) = lR(0, a) + w*(N(a)*Td +
                    k*(Nx(0,a)*Tx(0) + Nx(1,a)*Tx(1) + Nx(2,a)*Tx(2)));
            for (int b = 0; b < eNoN; b++){
                lK(0,a,b) = lK(0,a,b) + wl*(N(a)*N(b)*amd +
                        k*(Nx(0,a)*Nx(0,b) + Nx(1,a)*Nx(1,b) +
                        Nx(2,a)*Nx(2,b)));
            }
        }
    }
}