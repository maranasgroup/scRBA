# *sc*RBA: *S. cerevisiae* Resource Balance Analysis model
This repository provides resources and model files for the genome-scale resource balance analysis (RBA) model *sc*RBA that accompanies the following manuscript:<br>
**Evaluating proteome allocation of *Saccharomyces cerevisiae* phenotypes with resource balance analysis**<br>
Hoang V. Dinh and Costas D. Maranas<br>
Metabolic Engineering, 2023, 77:242-255; doi: https://doi.org/10.1016/j.ymben.2023.04.009<br>
<br>
Formulation and software usage guide are available at /scRBA/suppMat/scRBA_suppMethods_2022-09-12.docx<br>
Program and packages requirements for the provided software implementation: Python 3 (plus packages: cobra, pandas, numpy, matplotlib, jupyter, scikit-learn, cplex (compiled from installed CPLEX linear programming solver), and all automatically-installed associated dependencies to those already listed) and GAMS with ("built-in") Soplex solver. (Tested on Python 3.6 and GAMS 39.1.0. Python package installation could be done through "pip" and could take up to 1 hour)
