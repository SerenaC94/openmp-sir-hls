# HLS kernels for an SIR model solver

This repository contains source code, scripts, and evaluation results for solver kernels taken from an epidemiology SIR (Susceptible-Infected-Recovered) model and accelerated on FPGA through High-Level Synthesis.

The chosen HLS tool is Bambu, exploiting its SPARTA methodology for the synthesis of OpenMP applications.

# References

F. Ferrandi, V. G. Castellana, S. Curzel, P. Fezzardi, M. Fiorito, M. Lattuada, et al. “Invited:
Bambu: an Open-Source Research Framework for the High-Level Synthesis of Complex Applications”.
In: Proceedings of the 58th ACM/IEEE Design Automation Conference (DAC). 2021, pp. 1327–1330 

G. Gozzi, M. Fiorito, S.Curzel, C. Barone, V. G. Castellana, et al. “SPARTA: High-Level Synthesis
of Parallel Multi-Threaded Accelerators”. In: ACM Trans. Reconfigurable Technol. Syst. 18.1 (Dec. 2024)
