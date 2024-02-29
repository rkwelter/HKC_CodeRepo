Code implementing the HKC hierarchy of rotating Rayleigh-Benard convection models, as described in the paper 
"Rotating Rayleigh-Benard convection: Attractors, bifurcations and heat transport via a Galerkin hierarchy"

For details on the derivation of the model, see the paper: https://arxiv.org/abs/2402.14724

The code is organized into several directories, as follows:

1.) ModelConstruction: This is the natural starting place for working with the HKC hierarchy.  In this folder,
                       one finds the script "ModelConstructor.m".  By setting the model number M and running 
                       the script, a user can generate the time stepper for the Mth HKC model, which is then
                       written in the folder "FluidSolver/TimeSteppers".  The time steppers for several HKC
                       models have already been generated, but if the user desires to study others then they
                       can use this function.  The folder "ModelConstruction" contains several other similar
                       scripts, which write the same time stepper in different formats, such as for LaTeX and
                       Mathematica.
                       
2.) FluidSolver: In this folder, one finds the function "FluidSolver.m", which uses ode45 to call the time 
                 stepper written by ModelConstructor to solve a HKC model.  This is primarily used in via 
                 the functions in the HeatTransport folder, but can be used directly by first running the
                 script "AddPaths.m".  "FluidSolver.m" takes a set of parameter values, a previously
                 computed trajectory (can be empty), an "addTime" arugment, and a "save" boolean, solves the
                 corresponding HKC model starting from the last point on the previously computed trajectory
                 for the length of time "addTime", adds this to the previously computed trajectory, and if
                 the "save" boolean is true, then it saves the result in the folder "SavedData".  If the
                 previously computed trajectory is empty, then the function "InitialCondition.m" is called
                 to provide an initial condition.

3.) HeatTransport: This folder contains a number of functions for generating heat transport data from the
                   solutions to the HKC models.  More description to come.

4.) FluidVisualizer: In this folder, one finds the functions "FluidSimulator_TemperatureField.m",
                     "FluidSimulator_VelocityField.m", "FluidSimulator_PhaseSpace.m" and
                     "FluidSimulator_all.m".  

5.) SavedData: This folder contains the heat transport data in the subfolder "HeatTransport", as well as
               all trajectories generated in the subfolder "Trajectories".  

6.) SavedImages: This folder contains the images generated from the functions in "FluidVisualizer", as well
                 as several figures generated for the RRBC paper.
