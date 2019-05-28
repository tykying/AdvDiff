#!/bin/bash
#SIV=32
#echo \"${SIV}\"
#nohup matlab -nodisplay -nodesktop -r "SIv_in=${SIV}; Diffusivity_inference_from_MCMC" &> ./zkappa_$SIV.txt &

nohup matlab -nodisplay -nodesktop -r "Script_generate_particle_trajectories" &> ./all_traj2.txt &

