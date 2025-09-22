#include <math.h>
#include "PMIModels.h"

// Implementation Morozzi/Polzoni mobility model using the PMI interface

class PMI_Morozzi_Polzoni : public PMI_HighFieldMobility {

protected:
    const double q; // elementary charge
    const double F_min;
    double A_star, b; // A_star and b will be specified as the average of values obtained by simulations

public:
  PMI_Morozzi_Polzoni (const PMI_Environment& env,
                       const PMI_HighFieldDrivingForce force,
                       const PMI_AnisotropyType anisotype);

  ~PMI_Morozzi_Polzoni();
 

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

PMI_Morozzi_Polzoni::
PMI_Morozzi_Polzoni (const PMI_Environment& env,
                     const PMI_HighFieldDrivingForce force,
                     const PMI_AnisotropyType anisotype) :

PMI_HighFieldMobility (env, force, anisotype),

F_min (1.0e-06) // [V/cm]
{
}

PMI_Morozzi_Polzoni::
~PMI_Morozzi_Polzoni ()
{
}

void PMI_Morozzi_Polzoni::
Compute_mu (const double pot, const double n,
            const double p, const double t, const double ct,
            const double mulow, const double F, double& mu)

{
  double temp_mu = A_star * pow(pot, 4) * pow(t, 2);
  double exponent = b * (sqrt(abs(F) + F_min) / t);
  mu = temp_mu * exp(exponent);
}

void PMI_Morozzi_Polzoni::
Compute_dmudpot (const double pot, const double n,
                 const double p, const double t, const double ct,
                 const double mulow, const double F, double& dmudpot)

{
 double temp_dmudpot = 4 * A_star * pow(t, 2) * pow(pot, 3);
 double exp_1 = b * (sqrt(abs(F) + F_min) / t);
 dmudpot = temp_dmudpot * exp(exp_1);
}

void PMI_Morozzi_Polzoni::
Compute_dmudn (const double pot, const double n,
               const double p, const double t, const double ct,
               const double mulow, const double F, double& dmudn)

{
 dmudn = 0.0;
}

void PMI_Morozzi_Polzoni::
Compute_dmudp (const double pot, const double n,
               const double p, const double t, const double ct,
               const double mulow, const double F, double& dmudp)

{ dmudp = 0.0;
}

void PMI_Morozzi_Polzoni::
Compute_dmudt (const double pot, const double n,
               const double p, const double t, const double ct,
               const double mulow, const double F, double& dmudt)

{
 double temp_dmudt_1 = 2 * A_star * t * pow(pot, 4);
 double temp_dmudt_2 = A_star * sqrt(abs(F) + F_min) * pow(pot, 4);
 double exp_2 = b * (sqrt(abs(F) + F_min) / t);
 dmudt = (temp_dmudt_1 * exp(exp_2)) - (temp_dmudt_2 * exp(exp_2));
}

void PMI_Morozzi_Polzoni::
Compute_dmudct (const double pot, const double n,
                const double p, const double t, const double ct,
                const double mulow, const double F, double& dmudct)

{ dmudct = 0.0;
}

void PMI_Morozzi_Polzoni::
Compute_dmudmulow (const double pot, const double n,
                   const double p, const double t, const double ct,
                   const double mulow, const double F, double& dmudmulow)
{ dmudmulow = 0.0;
}

void PMI_Morozzi_Polzoni::
Compute_dmudF (const double pot, const double n,
               const double p, const double t, const double ct,
               const double mulow, const double F, double& dmudF)

{
 double temp_dmudF = (A_star * b * t * pow(pot, 4)) / (2 * sqrt(abs(F) + F_min));
 double exp_3 = b * (sqrt(abs(F) + F_min) / t);
 dmudF = temp_dmudF * exp(exp_3);
}

class PMI_e_Morozzi_Polzoni : public PMI_Morozzi_Polzoni {

public:
  PMI_e_Morozzi_Polzoni (const PMI_Environment& env,
                   const PMI_HighFieldDrivingForce force,
                   const PMI_AnisotropyType anisotype);

  ~PMI_e_Morozzi_Polzoni ()
  {
  }

};

PMI_e_Morozzi_Polzoni::
PMI_e_Morozzi_Polzoni (const PMI_Environment& env,
                       const PMI_HighFieldDrivingForce force,
                       const PMI_AnisotropyType anisotype) :

PMI_Morozzi_Polzoni (env, force, anisotype)

{ // default values
  A_star = InitParameter("A_star_e", 2.89e-14);
  b = InitParameter("b_e", 3.86e-04);
}

class PMI_h_Morozzi_Polzoni : public PMI_Morozzi_Polzoni {

public:
  PMI_h_Morozzi_Polzoni (const PMI_Environment& env,
                   const PMI_HighFieldDrivingForce force,
                   const PMI_AnisotropyType anisotype);

  ~PMI_h_Morozzi_Polzoni ()
  {
  }

};

PMI_h_Morozzi_Polzoni::
PMI_h_Morozzi_Polzoni (const PMI_Environment& env,
                 const PMI_HighFieldDrivingForce force,
                 const PMI_AnisotropyType anisotype) :

PMI_Morozzi_Polzoni (env, force, anisotype)

{ // default values
  A_star = InitParameter("A_star_h", 2.89e-14);
  b = InitParameter("b_h", 3.86e-04); 
}

extern "C"
PMI_HighFieldMobility* new_PMI_HighField_e_Mobility
  (const PMI_Environment& env, const PMI_HighFieldDrivingForce force,
   const PMI_AnisotropyType anisotype)

{ return new PMI_e_Morozzi_Polzoni (env, force, anisotype);
}

extern "C"
PMI_HighFieldMobility* new_PMI_HighField_h_Mobility
  (const PMI_Environment& env, const PMI_HighFieldDrivingForce force,
   const PMI_AnisotropyType anisotype)

{ return new PMI_h_Morozzi_Polzoni (env, force, anisotype);
}
