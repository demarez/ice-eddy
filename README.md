# ice-eddy
Ice-Eddy interaction project in collaboration with Josue

The github comprises the directories to put on the home directory

# Architecture of the project

## on the home dir, 

3 folders: INIT, RUN, DIAGS.

### in the INIT folder, 

create initial conditions using the *Initiate_vortex.ipynb* notebook: Choose parameters and name of the simulation

--> It creates a folder on datawork,  in the ...../init_cond/  repository, that contains *init_julia.nc*, the init file, and *parameters.txt* with all parameters of the simulation (eddy, grid, etc.)

### in the RUN folder, 

launch the julia simulation 

- adjust `path_root` in *run.jl*, and other parameterizations if needed.
- adjust name of simulation in *submit2GPU.pbs*
- $qsub submit2GPU.pbs

--> it will run the simulation on GPU, with nothing more needed because everything is written in the *parameters.txt* file
--> all will be stored on the datawork dir ....../RUN/ in the correct folder

note that there is a readme with some tricks, update it with new infos

### in the DIAGS folder, 

put notebooks to do diagnostics/plots or whatever wanted

  

## on the datawork dir, 

create 2 directories, init_cond and RUN. All the stuff will be stored there by the aforementioned actions. 
