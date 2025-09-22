#include <math.h>
#include "PMIModels.h"


// Implementation of the Poole-Frenkel mobility model (Wollongong) using the PMI interface

class PF_Wollongong : public PMI_HighFieldMobility {

protected:
    const double q; // elementary charge
    const double F_min; // minimum electric field to avoid convergence problems
    double A_star, b; // A_star and b will be specified as the average of values obtained by simulations

public:
  PF_Wollongong (const PMI_Environment& env,
                 const PMI_HighFieldDrivingForce force,
                 const PMI_AnisotropyType anisotype);

  ~PF_Wollongong ();
 

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

PF_Wollongong::
PF_Wollongong (const PMI_Environment& env,
               const PMI_HighFieldDrivingForce force,
               const PMI_AnisotropyType anisotype) :

PMI_HighFieldMobility (env, force, anisotype),

q (1.6021e-19), // [C]
F_min (1.0e-06) // [V/cm]
{
}

PF_Wollongong::
~PF_Wollongong ()
{
}

void PF_Wollongong::
Compute_mu (const double pot, const double n,
            const double p, const double t, const double ct,
            const double mulow, const double F, double& mu)

{
  double exponent = (b * t) * (sqrt(abs(F) + F_min));
  double temp_mu_h = (A_star * pot * pow(t, 2)) / (p * q * (abs(F) + F_min)); // no simplifications
  mu = temp_mu_h * (exp(exponent) - 1);
}

void PF_Wollongong::
Compute_dmudpot (const double pot, const double n,
                 const double p, const double t, const double ct,
                 const double mulow, const double F, double& dmudpot)

{
 double temp_dmudpot = (A_star * pow(t, 2)) / (p * q * (abs(F) + F_min));
 double exp_1 = (b * t) * (sqrt(abs(F) + F_min));
 dmudpot = temp_dmudpot * (exp(exp_1) - 1);
}

void PF_Wollongong::
Compute_dmudn (const double pot, const double n,
               const double p, const double t, const double ct,
               const double mulow, const double F, double& dmudn)

{
 dmudn = 0.0;
}

void PF_Wollongong::
Compute_dmudp (const double pot, const double n,
               const double p, const double t, const double ct,
               const double mulow, const double F, double& dmudp)

{
 double temp_dmudp = (A_star * pot * pow(t, 2)) / (pow(p, 2) * q * (abs(F) + F_min));
 double exp_2 = (b * t) * (sqrt(abs(F) + F_min));
 dmudp = (-temp_dmudp) * (exp(exp_2) - 1);
}

void PF_Wollongong::
Compute_dmudt (const double pot, const double n,
               const double p, const double t, const double ct,
               const double mulow, const double F, double& dmudt)

{
 double temp_dmudt_1 = (A_star * pot * pow(t, 2) * b) / (sqrt(abs(F) + F_min) * p * q);
 double temp_dmudt_2 = 2 * ((A_star * pot * t) / (p * q * (abs(F) + F_min)));
 double exp_3 = (b * t) * (sqrt(abs(F) + F_min));
 dmudt = (temp_dmudt_1 * exp(exp_3)) + (temp_dmudt_2 * (exp(exp_3) - 1));
}

void PF_Wollongong::
Compute_dmudct (const double pot, const double n,
                const double p, const double t, const double ct,
                const double mulow, const double F, double& dmudct)

{ dmudct = 0.0;
}

void PF_Wollongong::
Compute_dmudmulow (const double pot, const double n,
                   const double p, const double t, const double ct,
                   const double mulow, const double F, double& dmudmulow)
{ dmudmulow = 0.0;
}

void PF_Wollongong::
Compute_dmudF (const double pot, const double n,
               const double p, const double t, const double ct,
               const double mulow, const double F, double& dmudF)

{
 double temp_dmudF_1 = (A_star * pot * pow(t, 3) * b) / (2 * p * q * ((abs(F) + F_min) * sqrt(abs(F) + F_min)));
 double temp_dmudF_2 = (A_star * pot * pow(t, 2)) / (p * q * pow((abs(F) + F_min), 2));
 double exp_4 = (b * t) * (sqrt(abs(F) + F_min));
 dmudF = (temp_dmudF_1 * exp(exp_4)) - (temp_dmudF_2 * (exp(exp_4) - 1));
}

class PF_e_Wollongong : public PF_Wollongong {

public:
  PF_e_Wollongong (const PMI_Environment& env,
                   const PMI_HighFieldDrivingForce force,
                   const PMI_AnisotropyType anisotype);

  ~PF_e_Wollongong ()
  {
  }

};

PF_e_Wollongong::
PF_e_Wollongong (const PMI_Environment& env,
                 const PMI_HighFieldDrivingForce force,
                 const PMI_AnisotropyType anisotype) :

PF_Wollongong (env, force, anisotype)

{ // default values
  A_star = InitParameter("A_star_e", 2.89e-14);
  b = InitParameter("b_e", 3.86e-04);
}

class PF_h_Wollongong : public PF_Wollongong {

public:
  PF_h_Wollongong (const PMI_Environment& env,
                   const PMI_HighFieldDrivingForce force,
                   const PMI_AnisotropyType anisotype);

  ~PF_h_Wollongong ()
  {
  }

};

PF_h_Wollongong::
PF_h_Wollongong (const PMI_Environment& env,
                 const PMI_HighFieldDrivingForce force,
                 const PMI_AnisotropyType anisotype) :

PF_Wollongong (env, force, anisotype)

{ // default values
  A_star = InitParameter("A_star_h", 2.89e-14);
  b = InitParameter("b_h", 3.86e-04); 
}

extern "C"
PMI_HighFieldMobility* new_PMI_HighField_e_Mobility
  (const PMI_Environment& env, const PMI_HighFieldDrivingForce force,
   const PMI_AnisotropyType anisotype)

{ return new PF_e_Wollongong (env, force, anisotype);
}

extern "C"
PMI_HighFieldMobility* new_PMI_HighField_h_Mobility
  (const PMI_Environment& env, const PMI_HighFieldDrivingForce force,
   const PMI_AnisotropyType anisotype)

{ return new PF_h_Wollongong (env, force, anisotype);
}
