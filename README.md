# MagneticNanorod-Simulation-and-Analysis
- A collection of MATLAB codes for a capstone design project by Team Sensor Manjeom. This toolkit provides tools for both the theoretical simulation and experimental analysis of magnetic nanorod filler alignment in composites.
- Modeling &lt;nanorod alignment with Langevin Dynamics/EMT> and &lt;quantifying orientation from microscope images>

## 📝 Overview

This repository contains two main components that work together to characterize magnetic filler alignment:

1. **`simulation/`**: A MATLAB script that simulates the alignment dynamics of magnetic nanorods in a viscous fluid using the **Langevin equation** and predicts the composite's permittivity change with **Effective Medium Theory (EMT)**.
    
2. **`image_analysis/`**: A MATLAB script that automates the analysis of microscope images to detect fillers and **quantitatively measure their real-world orientation distribution**.
    

This allows for a direct comparison between theoretical predictions and experimental results.

---

## Part 1: Nanorod Dynamics Simulation (`simulation/`)

This script models the rotational behavior of magnetic nanorods under an external magnetic field, considering factors like magnetic torque, viscous drag, and Brownian motion.

!

### How to Use (Simulation)

1. Open `simulation/magnetic_nanorod_dynamics.m` in MATLAB.
    
2. Modify the variables in the **`USER INPUT PARAMETERS`** section to match your material properties and simulation conditions.
    
3. Run the script. An Excel file (`simulation_results.xlsx`) and two plot images will be generated inside the `simulation/` folder.
    

---

## Part 2: Microscope Image Analysis (`image_analysis/`)

This script analyzes a series of images to automatically measure the orientation of fillers, providing experimental data on their alignment.

!

### How to Use (Image Analysis)

1. Place all your microscope images inside the `image_analysis/Sample_image/` folder.
    
2. Open `image_analysis/OMimage_angle_analysis.m` in MATLAB.
    
3. In the **`USER PARAMETERS`** section, update the `params.image_files` list with the names of your image files.
    
4. Run the script. An intermediate plot showing the detected outlines and a final plot of the angle distribution will be generated.
    

---

## 🔧 Dependencies

- MATLAB (tested on R2021b and later)
    
- **Image Processing Toolbox** (required for the `image_analysis` script)
    

## 📜 License

This project is distributed under the [MIT License](https://www.google.com/search?q=LICENSE).
