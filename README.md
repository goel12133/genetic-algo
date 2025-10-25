This repository contains a **multi-physics simulation framework** for tokamak fusion reactors, combining **neutron transport modeling (OpenMC)** and **magnetic field optimization (MATLAB)**.  
The suite enables joint study of **magnetic confinement** and **neutronics performance**, supporting design of optimized submersion-style tokamak geometries with high tritium breeding ratios (TBR) and stable plasma confinement.


# Part 1: OpenMC Tokamak Neutronics Simulation


## Requirements
- Python 3.6+
- OpenMC
- Paramak
- openmc-data-downloader
- CAD-to-OpenMC converter
- Nuclear cross section data 

## Installation
```bash
git clone https://github.com/goel12133/openmc-paramak-opt.git
cd openmc-paramak-opt
pip install openmc paramak openmc-data-downloader
```

## Geometry Explanation

The simulation uses a **Submersion Tokamak** reactor geometry, created using the Paramak library. The reactor is defined by various radial thicknesses, coil placements, and other key parameters.

### Tokamak Geometry Parameters

| Parameter                            | Description                                      | Value  |
|--------------------------------------|--------------------------------------------------|--------|
| `inner_bore_radial_thickness`       | Radial thickness of the central bore            | 30 mm  |
| `inboard_tf_leg_radial_thickness`   | Thickness of the inboard toroidal field coil leg | 30 mm  |
| `center_column_shield_radial_thickness` | Thickness of the shielding around the center column | 30 mm  |
| `divertor_radial_thickness`         | Radial thickness of the divertor                | 80 mm  |
| `inner_plasma_gap_radial_thickness` | Gap between the plasma and center column        | 50 mm  |
| `plasma_radial_thickness`           | Radial thickness of the plasma                  | 200 mm |
| `outer_plasma_gap_radial_thickness` | Gap between the plasma and first wall           | 50 mm  |
| `firstwall_radial_thickness`        | Thickness of the first wall                     | 30 mm  |
| `blanket_rear_wall_radial_thickness` | Thickness of the blanket rear wall              | 30 mm  |
| `number_of_tf_coils`                | Number of toroidal field coils                  | 16     |
| `rotation_angle`                     | Rotation angle of the reactor geometry          | 360°   |
| `support_radial_thickness`          | Radial thickness of structural supports         | 90 mm  |

The toroidal and poloidal field coils are positioned strategically to ensure optimal plasma confinement. The reactor includes a vacuum boundary at a radius of **10 m**.

---

## Neutron Source Explanation

A **fixed neutron source** is defined using a Muir energy distribution to simulate **fusion-generated neutrons**. The source has the following properties:

### Neutron Source Parameters

| Parameter       | Description                                         | Value        |
|----------------|-----------------------------------------------------|-------------|
| `radius`       | Radial position of neutron emission                 | 293 mm      |
| `z_values`     | Vertical position of neutron emission               | 0 mm        |
| `angle`        | Angular distribution of emitted neutrons            | 0° to 90°   |
| `energy`       | Muir distribution peak energy                       | 14.1 MeV    |
| `m_rat`        | Mass ratio for energy distribution                   | 5.0         |
| `kt`          | Thermal spread of neutron energy distribution        | 20 keV      |

The neutron source is **isotropic** and **cylindrically distributed**, ensuring an even spread of fusion neutrons. These neutrons interact with the reactor components, and their behavior is analyzed using OpenMC tallies.

## Customizing

### Breeder Material Change

The breeder material in a fusion reactor is critical for neutron absorption and the production of tritium. In this model, the breeder material is initially set to **Li17Pb83**, a common lithium-lead alloy that is used for its ability to absorb neutrons and breed tritium. 

To customize the breeder material in this model, you can modify the material properties 



```python
mat_blanket = openmc.Material(name="blanket")
mat_blanket.add_elements_from_formula("Li4SiO4")  # New breeder material formula
mat_blanket.set_density("g/cm3", 2.0)  # Adjusted density for Li4SiO4
```

## Outputs

| Output               | Description                                                   |
|----------------------|---------------------------------------------------------------|
| `TBR tally`         | Tritium Breeding Ratio tally from neutron interactions        |
| `Neutron transport` | Simulation of neutron behavior and interactions in the tokamak |
| `Particle tracks`   | Monte Carlo tracking of neutron paths                         |
| `Energy deposition` | Heat deposition from neutron interactions                     |
| `Material reactions`| Nuclear reactions within materials (e.g., (n,Xt) reactions)  |


# Part 2: MATLAB Tokamak Magnetic Field Optimization with genetic algorithms



![Field Visualization](field.png)

This module extends the OpenMC-based neutronics model by introducing a **MATLAB-based magnetic confinement optimization**.  
It determines the **optimal toroidal and poloidal coil currents** that reproduce desired magnetic field profiles while maintaining plasma equilibrium and wall integrity — forming the electromagnetic half of the coupled tokamak simulation.

---

##  Overview

A tokamak’s performance depends on carefully balanced magnetic fields:
- **Toroidal field (Bₜₒᵣ):** produced by large circular coils surrounding the plasma; responsible for primary confinement.
- **Poloidal field (Bₚₒₗ):** generated by vertical field coils; stabilizes plasma position and shape.

This MATLAB script uses the **Biot–Savart law** to compute 3D magnetic fields from coil geometries and applies a **Genetic Algorithm (GA)** to find current values that achieve target field strengths and stability constraints.

---

##  Key Features

- Accurate **Biot–Savart magnetic field solver**
- **Toroidal and poloidal coil generation** for configurable geometry  
- **Genetic Algorithm (GA)** optimization of coil currents  
- Physically grounded **multi-term cost function** balancing:
  - Target field matching  
  - Wall confinement (B ≥ 1 T)  
  - Center field stability (B ≈ 0.2 T)
- Contour-map visualizations of optimized **toroidal** and **poloidal** magnetic fields  

---

##  Algorithm Description

### 1. Coil Geometry Generation
The geometry is defined by:
- Major radius `R₀ = 2.5 m`
- Minor radius `a = 1.0 m`
- `n_coils = 12` toroidal coils
- `n_poloidal = 6` poloidal coils per ring  

Toroidal coils are distributed evenly in φ, while poloidal coils are placed vertically along ±a.  
Positions and orientations are generated in `generate_coils()`.

### 2. Magnetic Field Formulation

The magnetic environment of a tokamak arises from the superposition of fields produced by both toroidal and poloidal coils.  
Each coil segment contributes an infinitesimal magnetic field determined by the **Biot–Savart law**, derived from Maxwell’s equations in magnetostatics:

$$
\mathbf{B}(\mathbf{r}) = \frac{\mu_0}{4\pi} 
\int_{\text{coil}} \frac{I\, d\mathbf{l} \times (\mathbf{r}-\mathbf{r}')}{\lVert \mathbf{r}-\mathbf{r}' \rVert^3}
$$

Here:
- $\mathbf{r}'$ denotes the position of an infinitesimal coil element $d\mathbf{l}$,
- $\mathbf{r}$ is the observation point,
- $I$ is the coil current, and  
- $\mu_0$ is the magnetic permeability of free space.

This integral encapsulates the **inverse-cube decay** of magnetic influence with distance and the **vector cross-product geometry** governing field directionality.  
In the simulation, the net field is evaluated as the linear superposition of all coil contributions at each point on a two-dimensional $(R, Z)$ cross-section of the tokamak.  
The resulting field map reveals regions of confinement, shear, and symmetry critical for plasma stability.

---

### 3. Optimization Functional

To achieve equilibrium confinement, the magnetic configuration must meet geometric, physical, and operational requirements at once.  
An **objective functional** is therefore used to measure how far the computed magnetic fields deviate from the desired configuration.

In simple terms, it combines several goals:
- Match the **toroidal** and **poloidal** magnetic fields to their target values.  
- Ensure the magnetic field at the reactor **walls** is strong enough to prevent plasma escape.  
- Keep the **central field** moderate to avoid over-compressing the plasma.

Each of these goals contributes to the overall optimization cost, with weighting factors that balance their importance.  
By minimizing this cost, the algorithm finds a set of coil currents that produces a stable, self-consistent magnetic equilibrium while respecting the reactor’s physical constraints.


Physically, these penalty terms prevent **magnetic leakage at the plasma boundary** and ensure the **central field remains weak enough** to avoid over-compression.  
The weights $\lambda_{\text{wall}} = 10^4$ and $\lambda_{\text{core}} = 500$ balance numerical sensitivity between global field alignment and local constraint enforcement.

Through minimization of $\mathcal{J}$, the optimization process seeks a **self-consistent current distribution** that yields magnetostatic equilibrium under engineering constraints.

---

### 4. Genetic Algorithm Optimization

Instead of relying on local gradient descent, the current configuration is optimized via a **Genetic Algorithm (GA)** — a population-based stochastic search inspired by evolutionary principles.  
Each candidate solution represents a vector of toroidal and poloidal coil currents:

$$
\mathbf{I} = [I_{\text{tor,1}},\,I_{\text{tor,2}},\ldots,I_{\text{pol,N}}]
$$

The GA iteratively evolves this population through selection, crossover, and mutation, favoring individuals that minimize the objective functional $\mathcal{J}$.  
This approach is robust to the **non-convex, multi-modal nature** of the magnetic equilibrium landscape and can capture subtle inter-coil couplings that affect field topology.

Typical hyperparameters are chosen to balance exploration and convergence:

```matlab
options = optimoptions('ga', ...
    'PopulationSize', 100, ...
    'MaxGenerations', 200, ...
    'FunctionTolerance', 1e-6, ...
    'CrossoverFraction', 0.8, ...
    'MutationFcn', {@mutationadaptfeasible, 0.2}, ...
    'Display', 'iter');
```

## License

This project is licensed under the [MIT License](https://opensource.org/licenses/MIT).  
© 2025 goel12133.


