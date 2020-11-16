path="/home/alexis/data/graph_data/docking_power/AGL_poses/1apw_training_poses/superposed"

for file in $(ls $path);
  do
    python graph_MLAI.py -i $path/$file -l LIG -gc -r 20.0 -ex
  done
