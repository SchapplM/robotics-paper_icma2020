# Sensitivity-based Model Reduction for In-Process Identification of Industrial Robots Inverse Dynamics

Authors:
Björn Volkmann, Daniel Kaczor, Mathias Tantau, Moritz Schappler, and Tobias Ortmaier

This repository provides supplemental material to the following publication at the IEEE ICMA 2020 conference:

```bib
@Conference{VolkmannKacTanSch2020,
  author    = {Volkmann, Björn and Kaczor, Daniel and Tantau, Mathias and Schappler, Moritz and Ortmaier, Tobias},
  title     = {Sensitivity-based Model Reduction for In-Process Identification of Industrial Robots Inverse Dynamics},
  booktitle = {2020 IEEE International Conference on Mechatronics and Automation (ICMA)}, 
  year      = {2020},
  month     = {October},
  pages     = {912-919},
  doi       = {10.1109/ICMA49215.2020.9233709},
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

### Prerequisites

The following scripts have to be run in Matlab to be able to re-create the tables:

* `identification/Model_Reduction.m`
* `identification/Optimal_Excitation.m`

### Tables and Figures

Matlab scripts for re-creating the data of the tables are located in directory `figures_tables`:

* Tab. I: `Create_Table_Parameters.m`
* Tab. II: `Create_Table_errorModelA.m`
* Tab. III: `Create_Table_errorModelB.m`
* Tab. IV: `Create_Table_ParametersProcess.m`
* Tab. V: `Create_Table_errorProcessOnProcess.m`
* Tab. VI: `Create_Table_errorModelAonProcess.m` (and not presented in the paper for model B `Create_Table_errorModelBonProcess.m`)
* Tab. VII: `Create_table_ControllError_comparison_Process1.m` (and not presented in the paper for process 2 `Create_table_ControllError_comparison_Process2.m`
* Fig. 4: `Create_Figure_squareEE.m`

### Robot Model

The directory `robot_model` contains the model of the Cloos QRC 350 robot used to obtain the results of the paper.

* Functions for creating the papers results:
  * `RegressorMatrix.m` provides the dynamics base parameter regressor X from equ. 3 of the paper
  * `WLS_method.m` implements equ. 5 of the paper
* The model is generated with the stand-alone Maple worksheet from `generation/main.mw`. The subdirectories `kinematics` and `dynamics` contain the generated output and intermediate results for validation. The directory `output` stores the automatically generated files from the worksheet.

Second model implementation as reference:
* The directory `codegen_reference` contains files to generate the robot model with another toolbox. Running the generation is not necessary at this point.
* `matlabfcn_CloosQRC350DE` contains the second implementation of the robot model with consideration of the differential gear mechanism in the hand joints. `matlabfcn_CloosQRC350OL` contains a simplified model only regarding the joints.
* `testfcn_CloosQRC350DE` and `testfcn_CloosQRC350OL` contain a unit test suite for the second implementation with the scripts `CloosQRC350DE_test_everything_fixbase.m` and `CloosQRC350OL_test_everything_fixbase.m`
* Both model implementations are compared against each other with the Matlab script `validation/validate_hybrdyn.m`, proving the validity of the model.
* The validation with the second implementation has the following dependencies:
   * [Matlab toolbox](https://github.com/SchapplM/matlab_toolbox) (tested with rev. 59c70af , 2020-05-25)
   * [Various external toolboxes](https://github.com/SchapplM/robotics-dep-ext) (tested with rev. 5a5fa3a, 2020-05-25)
   * [Robotics toolbox](https://github.com/SchapplM/robotics-toolbox) (tested with rev. 0e84d24, 2020-06-23)
* To regenerate the Matlab functions for the second implementation, the files in the directory `codegen_reference` can be used in connection to the [robot dynamics toolbox](https://github.com/SchapplM/robsynth-modelgen) from the co-author Moritz Schappler. The following steps are necessary
   * Download the toolbox from GitHub. Tested with rev. f50aeda, 2020-06-24. Set it up according to it's README file.
   * The code generation runs better on Linux operating system. For Windows an appropriate Bash terminal is necessary, see [toolbox Readme](https://github.com/SchapplM/robsynth-modelgen#unter-windows-windows-linux-subsystem).
   * The link to the toolbox in the file `hybrdyn_repo_path` has to be created (see template file with the same name)
   * Run `generate_models.sh` in a Bash terminal to create the second implementation of the model.
