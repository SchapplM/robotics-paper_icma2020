#!/bin/bash -ex
# Kopiere alle benötigten Dateien für die Code-Generierung ins HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
#
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-06
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
this_path=$(pwd)

if [ "$maplerepopath" == "" ]; then
  echo "Fehlendes Eingabeargument"
  exit 2
fi;
defpath=$maplerepopath/robot_codegen_definitions
constrpath=$maplerepopath/robot_codegen_constraints

## Definitionen kopieren
cp $this_path/robot_env_CloosQRC350DE $defpath/robot_env_CloosQRC350DE
cp $this_path/robot_env_CloosQRC350OL $defpath/robot_env_CloosQRC350OL

# Maple-Skripte (Kinematische Zwangsbedingungen)
cp $this_path/CloosQRC350DE_kinematic_constraints.mpl $constrpath/CloosQRC350DE_kinematic_constraints.mpl
cp $this_path/CloosQRC350DE_kinematic_constraints.mw  $constrpath/CloosQRC350DE_kinematic_constraints.mw

# Werte für Kinematikparameter (für Modultests)
cp $this_path/CloosQRC350_kinematic_parameter_values.m $constrpath/CloosQRC350OL_kinematic_parameter_values.m
cp $this_path/CloosQRC350_kinematic_parameter_values.m $constrpath/CloosQRC350DE_kinematic_parameter_values.m
