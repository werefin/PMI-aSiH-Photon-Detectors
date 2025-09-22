#include <math.h>
#include "PMIModels.h"

// Implementation of the Poole-Frenkel mobility model using the PMI interface

class PF_HighFieldMobility : public PMI_HighFieldMobility {

protected:
  const double k;
  double E0, beta, gamma;

public:
  PF_HighFieldMobility (const PMI_Environment& env,
                        const PMI_HighFieldDrivingForce force,
                        const PMI_AnisotropyType anisotype);

  ~PF_HighFieldMobility ();
 

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


PF_HighFieldMobility::
PF_HighFieldMobility (const PMI_Environment& env,
                      const PMI_HighFieldDrivingForce force,
                      const PMI_AnisotropyType anisotype) :

PMI_HighFieldMobility (env, force, anisotype),

k (8.617331e-05) // [eV/K]
{
}

PF_HighFieldMobility::
~PF_HighFieldMobility ()

{
}

void PF_HighFieldMobility::
Compute_mu (const double pot, const double n,
            const double p, const double t, const double ct,
            const double mulow, const double F, double& mu)

{
 double exp_1 = E0 / (k * t);
 double exp_2 = sqrt(F) * ((beta / t) - gamma);
 mu = mulow * exp(-exp_1) * exp(exp_2);
}


void PF_HighFieldMobility::
Compute_dmudpot (const double pot, const double n,
                 const double p, const double t, const double ct,
                 const double mulow, const double F, double& dmudpot)

{ dmudpot = 0.0;
}


void PF_HighFieldMobility::
Compute_dmudn (const double pot, const double n,
               const double p, const double t, const double ct,
               const double mulow, const double F, double& dmudn)

{ dmudn = 0.0;
}


void PF_HighFieldMobility::
Compute_dmudp (const double pot, const double n,
               const double p, const double t, const double ct,
               const double mulow, const double F, double& dmudp)

{ dmudp = 0.0;
}


void PF_HighFieldMobility::
Compute_dmudt (const double pot, const double n,
               const double p, const double t, const double ct,
               const double mulow, const double F, double& dmudt)

{ dmudt = 0.0;
}


void PF_HighFieldMobility::
Compute_dmudct (const double pot, const double n,
                const double p, const double t, const double ct,
                const double mulow, const double F, double& dmudct)

{ dmudct = 0.0;
}


void PF_HighFieldMobility::
Compute_dmudmulow (const double pot, const double n,
                   const double p, const double t, const double ct,
                   const double mulow, const double F, double& dmudmulow)
{ dmudmulow = 0.0;
}


void PF_HighFieldMobility::
Compute_dmudF (const double pot, const double n,
               const double p, const double t, const double ct,
               const double mulow, const double F, double& dmudF)

{ dmudF = 0.0;
}

class PF_e_HighFieldMobility : public PF_HighFieldMobility {

public:
PF_e_HighFieldMobility (const PMI_Environment& env,
                        const PMI_HighFieldDrivingForce force,
                        const PMI_AnisotropyType anisotype);
~PF_e_HighFieldMobility () {}

};

PF_e_HighFieldMobility::
PF_e_HighFieldMobility (const PMI_Environment& env,
                        const PMI_HighFieldDrivingForce force,
                        const PMI_AnisotropyType anisotype) :

PF_HighFieldMobility (env, force, anisotype)
{ // default values
  E0 = InitParameter ("E0_e", 0.5);
  beta = InitParameter ("beta_e", 1.1);
  gamma = InitParameter ("gamma_e", 0.1);
}

class PF_h_HighFieldMobility : public PF_HighFieldMobility {

public:
PF_h_HighFieldMobility (const PMI_Environment& env,
                        const PMI_HighFieldDrivingForce force,
                        const PMI_AnisotropyType anisotype);
~PF_h_HighFieldMobility () {}

};

PF_h_HighFieldMobility::
PF_h_HighFieldMobility (const PMI_Environment& env,
                        const PMI_HighFieldDrivingForce force,
                        const PMI_AnisotropyType anisotype) :

PF_HighFieldMobility (env, force, anisotype)
{ // default values
  E0 = InitParameter ("E0_h", 0.5);
  beta = InitParameter ("beta_h", 1.1);
  gamma = InitParameter ("gamma_h", 0.1);
}

extern "C"
PMI_HighFieldMobility* new_PMI_HighField_e_Mobility
  (const PMI_Environment& env, const PMI_HighFieldDrivingForce force,
   const PMI_AnisotropyType anisotype)

{ return new PF_e_HighFieldMobility (env, force, anisotype);
}

extern "C"
PMI_HighFieldMobility* new_PMI_HighField_h_Mobility
  (const PMI_Environment& env, const PMI_HighFieldDrivingForce force,
   const PMI_AnisotropyType anisotype)

{ return new PF_h_HighFieldMobility (env, force, anisotype);
}
