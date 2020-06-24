#!/bin/bash
# Kopiere generierten Code aus dem HybrDyn-Repo
# Dieses Skript im Ordner ausführen, in dem es liegt.
# Argument 1: Pfad zum HybrDyn-Repo

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-06
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

maplerepopath=$1
mdlext="DE OL"
for ext in $mdlext; do
  echo "Kopiere Code für CloosQRC350$ext"
  mkdir -p ../matlabfcn_CloosQRC350${ext}
  mkdir -p ../testfcn_CloosQRC350${ext}
  cp -u $maplerepopath/codeexport/CloosQRC350${ext}/testfcn/*.* ../testfcn_CloosQRC350${ext}
  cp -u $maplerepopath/codeexport/CloosQRC350${ext}/matlabfcn/*.* ../matlabfcn_CloosQRC350${ext}
  cp -u $maplerepopath/codeexport/CloosQRC350${ext}/tmp/robot*.sh ../matlabfcn_CloosQRC350${ext}
done
