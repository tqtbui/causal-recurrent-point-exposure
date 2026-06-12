This repository contains code for "Causal inference for the expected number of recurrent events in the presence of a terminal event"

The simulation results presented in the paper were obtained by running the code on a cluster server. To run the code on your own server, adjust the configurations in the 'sim_run.sh' file, create the appropriate directories for logs and outputs, and submit the job from the terminal using the command 

cd YOUR_DIRECTORY
sbatch sim_run.sh SCENARIO

where YOUR_DIRECTORY is where you store the code files and SCENARIO is the scenario you want to run, which is either 1, 2 or 3. 

You can also try running the simulation on a personal computer using the file 'small simulation example.R'. A data analysis example is also given in file 'data analysis example.R'. The comments in the two files make them self-explanatory.

