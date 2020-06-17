# Sensitivity-based Model Reduction for In-Process Identification of Industrial Robots Inverse Dynamics (Preprint)

Authors:

Björn Volkmann, Daniel Kaczor, Mathias Tantau, Moritz Schappler and Tobias Ortmaier

This repository provides supplemental material to the following submission to the IEEE ICMA 2020 conference:

```bib
@Conference{SchapplerOrt2020,
  author    = {Volkmann, Björn and Kaczor, Daniel and Tantau, Mathias and Schappler, Moritz and Ortmaier, Tobias},
  title     = {Sensitivity-based Model Reduction for In-Process Identification of Industrial Robots Inverse Dynamics},
  booktitle = {IEEE International Conference on Mechatronics and Automation},
  year      = {2020},
  month     = {October},
  note      = {submitted},
}
```

## Abstract

This paper presents a sensitivity-based approach for optimal model design and identification of the dynamics of a state-of-the-art industrial robot considering process-related restrictions.
The possibility of parameter excitation for subsequent identification of the model parameters is severely limited due to restrictions imposed by the process environment, especially the limited available workspace. 
Without sufficient parameter excitation, a satisfactory quality of the full model identification cannot be achieved, since non-excited parameters cannot be identified correctly.  
Furthermore, optimal excitation requires time-consuming calculations and distinct experiments during which the robot is not available for daily operation. 
It is therefore of interest to use process-related trajectories instead of dedicated excitation trajectories, which is expected to deteriorate the identifiability of the model parameters.
For this reason, the presented method uses a sensitivity-based approach allowing model order reduction in the identification process. The resulting model contains only those parameters excited by the excitation trajectory. 
For process-related trajectories this implies the model being limited to parameters relevant for the process. 
In experiments with a standard serial-link industrial robot controlled by standard industrial programmable logic control and servo inverters it is shown that the method produces significantly reduced models with a good measure of identifiability and quality.

## Contents

This repository contains Matlab scripts and Latex code to reproduce all figures and tables of the paper.
