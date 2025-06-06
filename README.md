# Coupled-Mode Theory in DFB Lasers – Numerical Simulation

This project implements numerical methods to solve coupled-mode equations modeling optical wave propagation in Distributed Feedback (DFB) lasers.

## 🔬 Project Overview

Coupled-mode theory describes the interaction between forward- and backward-propagating waves in periodic optical structures. In this project, we analyze such systems using:

- Finite difference methods for ODE and PDE models
- Numerical integration via Euler and Runge-Kutta schemes
- MATLAB simulations of field amplitude evolution and transmission spectra

## 🧪 Key Features

- Implementation of time-dependent and time-independent models
- Validation against analytical solutions for uniform gratings
- Transmission analysis as a function of detuning
- Convergence and stability checks

## 📁 Code Structure

- `methods_for_ODEs_space.m`: Solves the coupled ODE system
- `Transmission_ODEs.m`: Computes transmission and compares with theory
- `Transmission_PDEs.m`: Solves the full time-dependent PDE system
- `convergence_eul_rk4.m`: Checks convergence order of methods

## 📚 Tools & Environment

- MATLAB
- Based on coursework in *Engineering Physics* at Politecnico di Milano

## 👥 Author

- Marco Trento

---

