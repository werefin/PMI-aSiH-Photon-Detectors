## a-Si:H carrier mobility models in TCAD

This repository provides **Physical Model Interface (PMI)** implementations in C++ for **carrier mobility models** in hydrogenated amorphous silicon (a-Si:H). This work was developed for the **Synopsys Sentaurus TCAD** simulation environment, in collaboration with **INFN (Istituto Nazionale di Fisica Nucleare)** and the **University of Wollongong**.

### Overview

Hydrogenated amorphous silicon (a-Si:H) is promising for **particle detection** thanks to:
- Large-area coverage,
- Radiation damage resistance,
- Open challenges in electrical transport characterization.

Since standard TCAD lacks suitable models, we developed custom **PMI-based mobility models** to reproduce experimental I-V data and transient detector response.

### Implemented models

- **Poole–Frenkel mobility model** (baseline, adapted from organic semiconductors).  
- **Wollongong PF variant** (drift-based, modified for stability).  
- **Morozzi–Polzoni model** (best agreement with experimental data).  

Additional: **SRH DopingDependence lifetime model**, tuned for a-Si:H.

### Results

- Best reproduction of experimental I-V curves using **Morozzi–Polzoni + DopingDependence**.  
- **3D transient simulations** with Minimum Ionizing Particles show near-ideal **Charge Collection Efficiency** (i.e., $\text{CCE} \approx 100\%$).  
- Validated for **detector design studies** in a-Si:H.

All PMI models are collected in the `pmi` folder, while the final implementation that provides the best agreement with experimental data (Morozzi-Polzoni model) is provided separately in `pmi_PF_e_v2_1.C`.
