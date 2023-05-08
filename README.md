# Soft_Tissue_Meterial_Charectrization

The project's main objective is to determine the material constants ci using three distinct strain energy models (2-D Fung, 3-D Guccione et al., and 3-D Holzapfel, Gasser, and Ogden) based on provided experimental data. The input data comprises two separate files, each containing columns for P, ri, v, and Fz variables, with the initial line providing values for Ro, Ri, and B0.

To achieve this, the project first calculates force and pressure in theory and calculates the error between the theoretical and experimental values. The error is then minimized using gradient descent and Newton Raphson methods.

The project then generates plots comparing the experimental data with the results from the three different strain energy models. The objective is to determine which model produces the best result by analyzing the data and making comparisons.
