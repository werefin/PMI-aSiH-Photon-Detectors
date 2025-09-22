#include <math.h>
#include "PMIModels.h"


// Implementation of Poole-Frenkel mobility model using the PMI interface (Canali --> HighField_Mobility)

class Canali_HighFieldMobility : public PMI_HighFieldMobility {

protected:
  const double T0;
  double mu0, E0, k, T, exp_1, F, beta, gamma, exp_2;

public:
  Canali_HighFieldMobility (const PMI_Environment& env,
                            const PMI_HighFieldDrivingForce force,
                            const PMI_AnisotropyType anisotype);

  ~Canali_HighFieldMobility ();
  
 void Compute_internal (double E0, double k, double T, double& exp_1, double F, double beta, double gamma, double& exp_2);

  void Compute_mu
    (const double pot,        // electrostatic potential
     const double n,          // electron density
     const double p,          // hole density
     const double t,          // lattice temperature
     const double ct,         // carrier temperature
     const double mulow,      // low field mobility
     const double F,          // driving force
     double& mu);             // mobility

  void Compute_dmudpot
    (const double pot,        // electrostatic potential
     const double n,          // electron density
     const double p,          // hole density
     const double t,          // lattice temperature
     const double ct,         // carrier temperature
     const double mulow,      // low field mobility
     const double F,          // driving force
     double& dmudpot);        // derivative of mobility
                              // with respect to electrostatic potential

  void Compute_dmudn
    (const double pot,        // electrostatic potential
     const double n,          // electron density
     const double p,          // hole density
     const double t,          // lattice temperature
     const double ct,         // carrier temperature
     const double mulow,      // low field mobility
     const double F,          // driving force
     double& dmudn);          // derivative of mobility
                              // with respect to electron density

  void Compute_dmudp
    (const double pot,        // electrostatic potential
     const double n,          // electron density
     const double p,          // hole density
     const double t,          // lattice temperature
     const double ct,         // carrier temperature
     const double mulow,      // low field mobility
     const double F,          // driving force
     double& dmudp);          // derivative of mobility
                              // with respect to hole density

  void Compute_dmudt
    (const double pot,        // electrostatic potential
     const double n,          // electron density
     const double p,          // hole density
     const double t,          // lattice temperature
     const double ct,         // carrier temperature
     const double mulow,      // low field mobility
     const double F,          // driving force
     double& dmudt);          // derivative of mobility
                              // with respect to lattice temperature

  void Compute_dmudct
    (const double pot,        // electrostatic potential
     const double n,          // electron density
     const double p,          // hole density
     const double t,          // lattice temperature
     const double ct,         // carrier temperature
     const double mulow,      // low field mobility
     const double F,          // driving force
     double& dmudct);         // derivative of mobility
                              // with respect to carrier temperature

  void Compute_dmudmulow
    (const double pot,        // electrostatic potential
     const double n,          // electron density
     const double p,          // hole density
     const double t,          // lattice temperature
     const double ct,         // carrier temperature
     const double mulow,      // low field mobility
     const double F,          // driving force
     double& dmudmulow);      // derivative of mobility
                              // with respect to low field mobility

  void Compute_dmudF
    (const double pot,        // electrostatic potential
     const double n,          // electron density
     const double p,          // hole density
     const double t,          // lattice temperature
     const double ct,         // carrier temperature
     const double mulow,      // low field mobility
     const double F,          // driving force
     double& dmudF);          // derivative of mobility
                              // with respect to driving force

};

Canali_HighFieldMobility::
Canali_HighFieldMobility (const PMI_Environment& env,
                          const PMI_HighFieldDrivingForce force,
                          const PMI_AnisotropyType anisotype) :
  PMI_HighFieldMobility (env, force, anisotype),
  k (8.6167e-05) // [eV/K]
{
}



Canali_HighFieldMobility::
~Canali_HighFieldMobility ()

{
}

void Canali_HighFieldMobility::
Compute_internal (double E0, double k, double T, double& exp_1, double F, double beta, double gamma, double& exp_2) {
exp_1 = E0 / (k*T);
exp_2 = sqrt(F) * (beta/T - gamma);
}


void Canali_HighFieldMobility::
Compute_mu (const double pot, const double n,
            const double p, const double t, const double ct,
            const double mulow, const double F, double& mu)

{
 Compute_internal (E0, k, T, exp_1, F, beta, gamma, exp_2);
 mu = mu0 * exp(-exp_1) * exp(exp_2);
}



void Canali_HighFieldMobility::
Compute_dmudpot (const double pot, const double n,
                 const double p, const double t, const double ct,
                 const double mulow, const double F, double& dmudpot)

{ dmudpot = 0.0;
}



void Canali_HighFieldMobility::
Compute_dmudn (const double pot, const double n,
               const double p, const double t, const double ct,
               const double mulow, const double F, double& dmudn)

{ dmudn = 0.0;
}



void Canali_HighFieldMobility::
Compute_dmudp (const double pot, const double n,
               const double p, const double t, const double ct,
               const double mulow, const double F, double& dmudp)

{ dmudp = 0.0;
}



void Canali_HighFieldMobility::
Compute_dmudt (const double pot, const double n,
               const double p, const double t, const double ct,
               const double mulow, const double F, double& dmudt)

{ dmudt = 0.0;
}



void Canali_HighFieldMobility::
Compute_dmudct (const double pot, const double n,
                const double p, const double t, const double ct,
                const double mulow, const double F, double& dmudct)

{ dmudct = 0.0;
}



void Canali_HighFieldMobility::
Compute_dmudmulow (const double pot, const double n,
                   const double p, const double t, const double ct,
                   const double mulow, const double F, double& dmudmulow)
{ dmudmulow = 0.0;
}



void Canali_HighFieldMobility::
Compute_dmudF (const double pot, const double n,
               const double p, const double t, const double ct,
               const double mulow, const double F, double& dmudF)

{ dmudF = 0.0;
}


class Canali_e_HighFieldMobility : public Canali_HighFieldMobility {

public:
  Canali_e_HighFieldMobility (const PMI_Environment& env,
                              const PMI_HighFieldDrivingForce force,
                              const PMI_AnisotropyType anisotype);

  ~Canali_e_HighFieldMobility () {}

};


Canali_e_HighFieldMobility::
Canali_e_HighFieldMobility (const PMI_Environment& env,
                            const PMI_HighFieldDrivingForce force,
                            const PMI_AnisotropyType anisotype) :
  Canali_HighFieldMobility (env, force, anisotype)

{ // default values
  mu0 = InitParameter ("mu0_e", 100.00);
  T = InitParameter ("T_e", 313.15);
  E0 = InitParameter ("E0_e", 0.5);
  F = InitParameter ("F_e", 2.0e+05);
  beta = InitParameter ("beta_e", 1.1);
  gamma = InitParameter ("gamma_e", 0.1);
}


class Canali_h_HighFieldMobility : public Canali_HighFieldMobility {

public:
  Canali_h_HighFieldMobility (const PMI_Environment& env,
                              const PMI_HighFieldDrivingForce force,
                              const PMI_AnisotropyType anisotype);

  ~Canali_h_HighFieldMobility () {}

};


Canali_h_HighFieldMobility::
Canali_h_HighFieldMobility (const PMI_Environment& env,
                            const PMI_HighFieldDrivingForce force,
                            const PMI_AnisotropyType anisotype) :
  Canali_HighFieldMobility (env, force, anisotype)

{ // default values
  mu0 = InitParameter ("mu0_h", 100.00);
  T = InitParameter ("T_h", 313.15);
  E0 = InitParameter ("E0_h", 0.5);
  F = InitParameter ("F_h", 2.0e+05);
  beta = InitParameter ("beta_h", 1.1);
  gamma = InitParameter ("gamma_h", 0.1);
  
}


extern "C"
PMI_HighFieldMobility* new_PMI_HighField_e_Mobility
  (const PMI_Environment& env, const PMI_HighFieldDrivingForce force,
   const PMI_AnisotropyType anisotype)

{ return new Canali_e_HighFieldMobility (env, force, anisotype);
}


extern "C"
PMI_HighFieldMobility* new_PMI_HighField_h_Mobility
  (const PMI_Environment& env, const PMI_HighFieldDrivingForce force,
   const PMI_AnisotropyType anisotype)

{ return new Canali_h_HighFieldMobility (env, force, anisotype);
}
