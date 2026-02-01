# Binary Polymer-Solvent Thermodynamic Modeling

A comprehensive Python framework for modeling and simulating thermodynamic properties of binary polymer-gas (solvent) systems using the Sanchez-Lacombe Equation of State (EOS). This toolkit enables researchers to predict solubility, swelling behavior, and glass transition temperatures for polymer-CO₂ systems.

## Overview

This project implements advanced thermodynamic models for predicting phase equilibrium and transport properties in binary polymer-solvent mixtures. The framework is particularly focused on polymer-CO₂ systems and uses the Sanchez-Lacombe EOS combined with various mixing rules to calculate:

- Gas solubility in polymers (mass fraction absorbed)
- Polymer swelling behavior (volume changes upon gas absorption)
- Pressure-dependent glass transition temperature Tg(P)
- Phase equilibrium using chemical potential equality
- Entropy of mixing for polymer-solvent systems

The code supports multiple polymer types and includes utilities for fitting theoretical models to experimental data, with built-in parameters from various published research papers.

## Features

### Core Capabilities
- **Multiple Polymer Support**: PMMA, PS (Polystyrene), PC (Polycarbonate), PLA (Polylactic acid), and more
- **Sanchez-Lacombe EOS Implementation**: Classical lattice-fluid theory for polymer solutions
- **Multiple Mixing Rules**: Hassan, Kier, and Condo variations for different modeling scenarios
- **Bisection-Based Solving**: Robust numerical methods for phase equilibrium calculations
- **Reference Volume Correction**: Tait equation implementation for pressure-volume-temperature (PVT) corrections

### Analysis Tools
- **Glass Transition Temperature Modeling**: Calculate Tg as a function of pressure and solvent concentration
- **Solubility and Swelling Predictions**: Determine gas absorption and resulting polymer volume changes
- **Entropy Calculations**: Compute mixing entropy for thermodynamic analysis
- **Experimental Data Fitting**: Least-squares optimization to fit model parameters to experimental data
- **Isotherm/Isobar Analysis**: Split and analyze data at constant temperature or pressure
- **Discontinuity Detection**: Identify phase transitions in the equation of state

### Visualization
- **2D Plotting**: Generate publication-quality plots of Tg(P), solubility, and swelling curves
- **3D Surface Visualization**: Create 3D plots of thermodynamic surfaces
- **Comparative Analysis**: Plot multiple isotherms and polymer types for comparison

### Data Management
- **Experimental Data Loading**: Import and process experimental solubility/swelling data
- **CO₂ PVT Interpolation**: High-accuracy CO₂ property calculations
- **Data Splitting Utilities**: Separate data above/below Tg and into isothermal sets
- **Colored Console Output**: Enhanced terminal feedback for debugging

## Installation

### Prerequisites
- Python 2.7 (legacy code base - **Note**: Python 2.7 reached end-of-life in 2020. Migration to Python 3.x is recommended for security and compatibility)
- NumPy
- SciPy
- Matplotlib
- SymPy

### Setup

1. Clone the repository:
```bash
git clone https://github.com/hassanalam014/binary-all-polymer-types-reference-volume-changing.git
cd binary-all-polymer-types-reference-volume-changing
```

2. Install required dependencies:
```bash
pip install numpy scipy matplotlib sympy
```

3. Verify installation by running a test script:
```bash
python plot_Tg_nsolve_POST_THESIS.py
```

## Usage

### Basic Example: Solubility and Swelling Calculations

```python
from Parameters_of_Different_Polymers import Parameters_of_Different_Polymers
from Parameters_for_Mixtures_and_Tg import Parameters_for_Mixtures_and_Tg
from All_Functions import *

# Configure system
Polymer_Type = 'PS'  # Polystyrene
Solvent = 'CO2'
Parameters_Paper = 'Self_Grassia'
Paper_Number = 'Paper4_11_12'
Cp_Polymer_Weight = '02kilo_POST_THESIS'

kwargs = {
    'Polymer_Type': Polymer_Type,
    'Solvent': Solvent,
    'Parameters_Paper': Parameters_Paper,
    'Paper_Number': Paper_Number,
    'Cp_Polymer_Weight': Cp_Polymer_Weight
}

# Get polymer and solvent parameters
Ppstar, Tpstar, Rpstar, Mp, Psstar, Tsstar, Rsstar, Ms, P_exp, Tg_exp = \
    Parameters_of_Different_Polymers(**kwargs)

# Calculate mixture properties
# (see specific scripts for detailed examples)
```

### Running Main Analysis Scripts

#### 1. Glass Transition Temperature Analysis
```bash
python plot_Tg_nsolve_POST_THESIS.py
```
This script calculates and plots the glass transition temperature as a function of pressure for polymer-CO₂ mixtures.

#### 2. Solubility and Swelling Predictions
```bash
python plot_solubility_swelling_phi_through_bisect_improved_modern.py
```
Computes solubility (mass fraction) and swelling (volume change) across temperature and pressure ranges.

#### 3. Entropy Analysis
```bash
python plot_entropy_nsolve_POST_THESIS.py
```
Calculates mixing entropy for thermodynamic analysis of polymer-solvent systems.

#### 4. Parameter Fitting to Experimental Data
```bash
python fit_mixtureDHV_bisect_method.py
```
Fits mixture interaction parameters (delta, zeta) to experimental solubility/swelling data.

#### 5. 3D Visualization
```bash
python 3Dplot.py
```
Creates 3D surface plots of thermodynamic properties.

## Configuration

### Polymer Selection
Edit the main script to select polymer type:
```python
Polymer_Type = 'PS'      # Options: 'PS', 'PMMA', 'PC', 'PLA'
Solvent = 'CO2'          # Currently optimized for CO2
```

### Parameter Source Selection
Choose which research paper's parameters to use:
```python
Parameters_Paper = 'Self_Grassia'  # Options: 'Self_Grassia', 'Self_Park', 'Self_Schmidt', etc.
```

### Data Source Configuration
Specify experimental data reference:
```python
Paper_Number = 'Paper4_11_12'  # For PS
# Paper_Number = 'Paper15'     # For PMMA
```

### Mixing Rule Selection
In the calculation functions, select the mixing rule:
```python
Kier = True          # Use Kier's mixing rule
Hassan = False       # Hassan's modified approach
Condo = False        # Condo's approach
Hassan_Var_Vol = False  # Hassan with variable volume
```

### Adjustable Parameters
Key parameters that can be tuned:
- **delta**: Hole volume ratio in mixture vs pure solvent (typically 0.95-1.05)
- **zeta**: Cross-interaction energy parameter (typically 0.95-1.05)
- **g, epsilon_p, x**: Heat capacity parameters for Tg calculations

## Folder Structure

```
binary-all-polymer-types-reference-volume-changing/
├── All_Functions.py                      # Core EOS and thermodynamic functions
├── Parameters_of_Different_Polymers.py   # Sanchez-Lacombe parameters for polymers
├── Parameters_for_Mixtures_and_Tg.py     # Mixture parameters and Tg-related values
├── Tait_Parameters_of_Different_Polymers.py  # PVT correction parameters
│
├── loadExperimentalData.py               # Load solubility/swelling data
├── loadExperimentalDataCO2.py            # Load CO2-specific experimental data
├── CO2PVT_interpolation.py               # CO2 property interpolation
│
├── plot_Tg_nsolve_POST_THESIS.py         # Main Tg(P) calculation and plotting
├── plot_entropy_nsolve_POST_THESIS.py    # Entropy calculation and analysis
├── plot_solubility_swelling_phi_through_bisect_improved_modern.py  # Solubility/swelling predictions
├── plot_mixture.py                       # General mixture property visualization
├── plot_glass_temperature_Retro_improved.py  # Alternative Tg calculation method
├── 3Dplot.py                             # 3D thermodynamic surface visualization
│
├── fit_mixtureDHV_bisect_method.py       # Parameter fitting (bisection method)
├── fit_mixtureDHV_Kier_Original_Program_nsolve_method.py  # Parameter fitting (nsolve method)
│
├── SplittingExperimentalDataOf_X_Sw_aboveANDbelowTgANDIsotherms.py  # Data splitting utilities
├── Split_Exp_Data_in_Isotherms.py        # Isotherm extraction
├── find_discontinuity.py                 # Phase transition detection
│
├── To_get_colored_print.py               # Colored terminal output utilities
├── Question1.py                          # Specific analysis question/test
├── Random Testing.py                     # Development/debugging script
├── Tg(P) For Direct Plotting from Data Points.py  # Direct Tg data visualization
│
├── LICENSE                               # MIT License
└── README.md                             # This file
```

### Key File Descriptions

**Core Modules:**
- `All_Functions.py`: Contains fundamental equation of state implementations, chemical potential calculations, and phase equilibrium solvers
- `Parameters_of_Different_Polymers.py`: Database of Sanchez-Lacombe characteristic parameters (P*, T*, R*) for various polymers from literature
- `Parameters_for_Mixtures_and_Tg.py`: Mixture interaction parameters and glass transition-related thermodynamic properties

**Data Handling:**
- `loadExperimentalData.py`: Imports experimental solubility and swelling measurements
- `CO2PVT_interpolation.py`: Provides high-accuracy CO₂ density and volume calculations

**Analysis Scripts:**
- `plot_Tg_nsolve_POST_THESIS.py`: Primary tool for calculating pressure-dependent glass transition temperatures
- `plot_solubility_swelling_phi_through_bisect_improved_modern.py`: Main solubility/swelling prediction engine
- `fit_mixtureDHV_bisect_method.py`: Optimization routine for fitting model parameters to data

**Utilities:**
- `Split_Exp_Data_in_Isotherms.py`: Separates experimental data into constant temperature sets
- `find_discontinuity.py`: Detects liquid-vapor or glass-rubber transitions
- `To_get_colored_print.py`: Provides colored console output for enhanced debugging

## Contributing

Contributions are welcome! This is a research-oriented codebase, so please consider the following:

### How to Contribute

1. **Fork the repository**
2. **Create a feature branch**: `git checkout -b feature/your-feature-name`
3. **Make your changes**: Ensure code follows existing style conventions
4. **Test thoroughly**: Verify calculations against known results
5. **Document changes**: Update docstrings and comments
6. **Submit a pull request**: Include a clear description of changes and their scientific motivation

### Contribution Guidelines

- **Maintain scientific accuracy**: All changes should be theoretically sound and well-referenced
- **Preserve backward compatibility**: This code may be used in ongoing research projects
- **Add appropriate comments**: Explain the physical/mathematical meaning of calculations
- **Include references**: Cite papers for new models or parameter sets
- **Test with multiple polymers**: Verify changes work across different polymer types

### Areas for Contribution

- **Python 3 migration**: Update codebase to Python 3.x
- **Additional polymers**: Add parameters for new polymer types
- **Performance optimization**: Improve numerical solver efficiency
- **Unit tests**: Add comprehensive test coverage
- **Documentation**: Improve inline documentation and add examples
- **New EOS models**: Implement alternative equations of state
- **GUI development**: Create graphical interface for easier parameter selection

### Reporting Issues

If you encounter bugs or have feature requests:
- Open an issue on GitHub with a clear description
- Include error messages and relevant code snippets
- Specify polymer type, parameters, and conditions that trigger the issue
- Provide expected vs. actual behavior

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

```
MIT License

Copyright (c) 2018 Hassan Alam

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## Citation

If you use this code in your research, please cite:

```
Alam, H. (2018). Binary Polymer-Solvent Thermodynamic Modeling Framework.
GitHub repository: https://github.com/hassanalam014/binary-all-polymer-types-reference-volume-changing
```

## Acknowledgments

This project implements and extends thermodynamic models from various research papers, including work by:
- Sanchez and Lacombe (original equation of state)
- Kier, von Konigslow, and contributors
- Park, Schmidt, Grassia (parameter development)

Special thanks to the polymer thermodynamics research community for experimental data and theoretical foundations.

## Contact

For questions, suggestions, or collaboration opportunities:
- **Author**: Hassan Alam
- **Repository**: https://github.com/hassanalam014/binary-all-polymer-types-reference-volume-changing
- **Issues**: https://github.com/hassanalam014/binary-all-polymer-types-reference-volume-changing/issues

---

**Note**: This is research code developed for academic purposes. While efforts have been made to ensure accuracy, users should validate results against experimental data for their specific applications.
