# Aeroelastic Panel Flutter Model

Numerical model for the aeroelastic analysis of thin plates subjected to supersonic flow, combining **Kirchhoff Plate Theory**, the **Finite Element Method (FEM)** and aerodynamic loading based on **Piston Theory**.

The project focuses on the study of **panel flutter**, structural dynamic stability and fluid-structure interaction.

---

## Overview

This repository contains a computational framework for the aeroelastic analysis of thin plates exposed to high-speed aerodynamic flow.

The structural model is based on Kirchhoff plate theory and discretized using the Finite Element Method.

The aerodynamic contribution is introduced through Piston Theory, allowing the investigation of the interaction between structural vibration and aerodynamic loading.

The model can be used to evaluate:

- Natural frequencies;
- Mode shapes;
- Aeroelastic eigenvalues;
- Flutter velocity;
- Modal damping;
- Frequency evolution with flow velocity;
- Structural instability;
- Frequency Response Functions;
- Influence of aerodynamic parameters;
- Sensitivity to structural properties.

---

## Physical Model

The structure is modeled as a thin plate governed by Kirchhoff Plate Theory.

Main assumptions include:

- Thin plate behavior;
- Small transverse displacements;
- Linear elastic structural response;
- Negligible transverse shear deformation;
- Supersonic aerodynamic flow;
- Aerodynamic loading modeled using Piston Theory.

Each node contains three mechanical degrees of freedom:

- Transverse displacement \(w\);
- Rotation \(\theta_x\);
- Rotation \(\theta_y\).

---

## Structural Dynamic Equation

The structural system can be expressed as:

\[
\mathbf{M}\ddot{\mathbf{q}}
+
\mathbf{C}\dot{\mathbf{q}}
+
\mathbf{K}\mathbf{q}
=
\mathbf{F}_{aero}
\]

where:

- \(\mathbf{M}\) is the global mass matrix;
- \(\mathbf{C}\) is the structural damping matrix;
- \(\mathbf{K}\) is the structural stiffness matrix;
- \(\mathbf{q}\) contains the structural degrees of freedom;
- \(\mathbf{F}_{aero}\) represents the aerodynamic forces.

When aerodynamic contributions are introduced, the aeroelastic system can be written in a generalized form as:

\[
\mathbf{M}\ddot{\mathbf{q}}
+
\left(\mathbf{C}+\mathbf{C}_{aero}\right)\dot{\mathbf{q}}
+
\left(\mathbf{K}+\mathbf{K}_{aero}\right)\mathbf{q}
=
0
\]

where:

- \(\mathbf{C}_{aero}\) represents aerodynamic damping effects;
- \(\mathbf{K}_{aero}\) represents aerodynamic stiffness effects.

The resulting system depends on the flow velocity and is used to evaluate aeroelastic stability.

---

## Aerodynamic Model

The aerodynamic loading is modeled using **Piston Theory**, which provides an approximation for aerodynamic pressure in supersonic flow.

The aerodynamic pressure depends on the structural deformation and its time variation.

A generic linearized form can be represented as:

\[
\Delta p
=
A(U)\frac{\partial w}{\partial x}
+
B(U)\frac{\partial w}{\partial t}
\]

where:

- \(w\) is the transverse displacement;
- \(U\) is the flow velocity;
- \(A(U)\) represents the aerodynamic stiffness contribution;
- \(B(U)\) represents the aerodynamic damping contribution.

These aerodynamic terms are converted into finite element matrices and coupled with the structural system.

---

## Finite Element Formulation

The plate is discretized using rectangular Kirchhoff plate elements.

The numerical procedure includes:

1. Definition of plate geometry;
2. Definition of material properties;
3. Mesh generation;
4. Calculation of element mass matrices;
5. Calculation of structural stiffness matrices;
6. Global matrix assembly;
7. Application of boundary conditions;
8. Modal analysis;
9. Modal reduction;
10. Assembly of aerodynamic matrices;
11. Aeroelastic system formulation;
12. Flow velocity sweep;
13. Eigenvalue calculation;
14. Flutter detection;
15. Post-processing of frequencies and damping.

---

## Aeroelastic Stability

Flutter is identified through the evolution of the system eigenvalues as the flow velocity increases.

The aeroelastic eigenvalues can be expressed as:

\[
\lambda = \sigma + i\omega
\]

where:

- \(\sigma\) represents the real part of the eigenvalue;
- \(\omega\) represents the oscillation frequency.

The system is stable when:

\[
\sigma < 0
\]

Flutter onset occurs when one pair of eigenvalues crosses the stability boundary:

\[
\sigma = 0
\]

and becomes unstable for:

\[
\sigma > 0
\]

The corresponding flow velocity is defined as the **critical flutter velocity**.

---

## Modal Reduction

To reduce computational cost, the structural system can be projected into modal coordinates.

The physical displacement vector is approximated as:

\[
\mathbf{q} = \mathbf{\Phi}\mathbf{\eta}
\]

where:

- \(\mathbf{\Phi}\) contains selected structural mode shapes;
- \(\mathbf{\eta}\) contains the modal coordinates.

The reduced system becomes:

\[
\mathbf{M}_r\ddot{\mathbf{\eta}}
+
\mathbf{C}_r\dot{\mathbf{\eta}}
+
\mathbf{K}_r\mathbf{\eta}
=
\mathbf{F}_{aero,r}
\]

This approach significantly reduces the size of the aeroelastic eigenvalue problem.

---

## Main Features

### Structural Analysis

- Finite Element Method implementation;
- Kirchhoff thin plate formulation;
- Natural frequency calculation;
- Mode shape extraction;
- Structural stiffness and mass matrix assembly;
- Modal reduction.

### Aeroelastic Analysis

- Supersonic aerodynamic loading;
- Piston Theory implementation;
- Aerodynamic stiffness matrix;
- Aerodynamic damping matrix;
- Aeroelastic eigenvalue analysis;
- Flutter velocity calculation;
- Stability boundary identification.

### Post-Processing

- Frequency vs. flow velocity plots;
- Damping vs. flow velocity plots;
- Eigenvalue tracking;
- Flutter mode identification;
- Mode shape visualization;
- Parametric studies;
- Sensitivity analysis.

---

## Typical Results

The model can generate:

- Structural natural frequencies;
- Structural mode shapes;
- Aeroelastic frequencies;
- Aeroelastic damping;
- Real and imaginary parts of eigenvalues;
- Critical flutter velocity;
- Flutter mode identification;
- Frequency coalescence behavior;
- Velocity-dependent stability diagrams.

Example outputs:

```text
results/
├── natural_frequencies/
├── mode_shapes/
├── frequency_velocity/
├── damping_velocity/
├── eigenvalues/
└── flutter_boundary/
