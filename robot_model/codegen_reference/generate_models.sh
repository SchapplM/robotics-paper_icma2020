#!/bin/bash -ex
# Generate the robot models
# execute in the containing folder

# Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-06
# (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

this_path=$(pwd)
hybrdyn_repo_path=`cat hybrdyn_repo_path`

if [ "$hybrdyn_repo_path" == "" ]; then
  echo "you need to place a file named hybrdyn_repo_path at this location; see template file"
fi;

cd $this_path
./prepare_maple_repo.sh $hybrdyn_repo_path

cd $this_path
./generate_maple_code.sh $hybrdyn_repo_path

cd $this_path
./copy_generated_code.sh $hybrdyn_repo_path

