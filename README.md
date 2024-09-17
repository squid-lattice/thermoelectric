--------------TKB.py-------------
This code can be used to generate a Pisarenko (carrier concentration [carriers/cm^3] vs Seebeck coefficient [uV/K]) plot using a two band Kane model. 
Lines 25 - 37 will need to edited for your material. The E_gap, K, etc values are currently for SnTe.
Coded using math from: Wiendlocha, Candolfi et al, 2021, Residual resistivity as an independent indicator of resonant levels in semiconductors and Zhang, Ren 2013 PNAS, High thermoelectric performance by resonant dopant indium in nanostructured SnTe	9
--------------MOFC_SPB.py-------------
This code uses a single parabolic band model (SPB) to solve for the following four material parameters: reduced Fermi level, scattering time (tau), scattering prefactor (r), and carrier effective mass (m*/me). You must provide experimental input (measured at the same temperature): resistivity, Seebeck coefficient, Hall carrier concentration, and Nernst coefficient. Units are specified in the code. Output is generated as a .csv file saved to the directory where the .py file is stored on your computer.
