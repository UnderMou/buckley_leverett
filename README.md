# Buckley-Leverett Equation Solver

This repository contains implementations and comparative analyses of the Buckley-Leverett equation solution using both analytical and numerical methods. The code is organized into four main folders:

## 1. buckley_leverett_analytical

This folder contains the full implementation of the analytical construction solution of the Buckley-Leverett equation. Each module focuses on different aspects:
- `fractional_flow.py`: Module for fractional flow functions.
- `relative_permeability.py`: Module for relative permeability curves.
- `bl_solution.py`: Module for the analytical solution construction procedure.
- `recovery_calc.py`: Script for calculating recovery using the analytical solution.
- `main.py`: Main application code that calls the modules in an example.

## 2. buckley_leverett_upwind

In this folder, the Buckley-Leverett equation solution is implemented in a numerical way using finite differences. The numerical approach includes upwind handling of the advective term. It utilizes the modules of fractional flow and relative permeability described in the analytical section.

## 3. compare_BL_analytical_numerical

This folder is a duplicate of the previous two and is dedicated to a comparative analysis of solutions obtained by both methods (analytical and numerical). It provides insights into the differences and similarities between the two approaches.

## 4. others

The `others` folder serves as a backup for additional codes and tests developed throughout the entire implementation. It contains miscellaneous files and experimental code that may not be directly related to the main Buckley-Leverett equation solutions but are valuable for reference and future development.

## Usage

Follow the instructions in each folder's respective README or code comments for details on how to use and run the provided implementations.

