%  ========================================================================
%> @file  matlab_example.m
%>
%> @brief Brief example illustration how to run the cellular automaton,
%>        which is implemented in C++, in MATLAB.
%  ========================================================================
%>
%> @brief Brief example illustration how to run the cellular automaton,
%>        which is implemented in C++, in MATLAB.
%>
%> It allows to set the size of the domain, the porosity, the
%> jump_parameter, indicating how far cells are allowed to jump, and the
%> number of steps the cellular automaton should run. It creates the
%> resulting data and images visualizing the results.
%> 
%> This file is part of the GitHub repository
%>   https://github.com/AndreasRupp/cellular-automaton
%> Copyright and license conditions can be found there.

%% First step is including the directory 'mex/CAM' to the MATLAB path.
%  This allows to run all MATLAB functions of CAM. These MATLAB functions
%  take care of compiling and running respective C++ functions, in which
%  the actual work is done.
addpath('mex/CAM');
addpath('matlab_functions');

%% Define the arguments of the run_cam function.
%  This example illustrates how to run the run_cam function, see the file
%  mex/CAM/run_cam.m. It runs a cellular automaton method on a grid of size
%  'nx(1) * nx(2) * ... * nx(n)' for 'num_steps' steps. Moreover, we set
%  the porosity (i.e.the percentage of void space---not occupied by solid)
%  and the jump parameter (telling how far individual cells are allowd to
%  jump) as mandatory parameters.
nx             = [10 10];
num_steps      = 10;
porosity       = 0.9;
jump_parameter = 1;

%  Optionally, we can set a frame rate (saying after how many steps an
%  output shall be returned), a boolean parameter indicating if the output
%  of the results shall be returned or not, a boolean parameter indicating
%  if the calculated measures shall be returned or not and a boolean
%  parameter indicating if the number of the random seed shall be printed
%  or not. The parameter num_random_seed allows to choose between a random
%  seed, depending on the clock time (if set to 0), or a given seed.
frame_rate         = 1;
output_results     = true;
output_measures    = true;
output_random_seed = true;
num_random_seed    = 0;

%% Call the run_cam function and save the output data.
%  This calls the run_cam function with the chosen parameters. The output
%  consists of the matrix domain_data with nx(1)*nx(2)*...*nx(n) rows and
%  a column for every step with output. The entries of domain_data are
%  either 0, corresponding to void cells, i.e. non-solid cells, or
%  positive integers corresponding to solid cells. The matrix measures
%  contains 6 rows and a column for every step with output. Every row
%  corresponds to a certain geometric measure of the domain_data at the
%  given step: 
%  - row 1: Number of single solid pixels without solid neighbours.
%  - row 2: Number of solid particles, including single solid pixels and
%  agglomorates of solid pixels.
%  - row 3: Total number of solid pixels. This is constant over the run of
%  one simulation.
%  - row 4: Total solid surface, i.e. the sum over all solid pixels of
%  neighboring void pixels.
%  - row 5: Mean particle size, i.e. the total number of solid pixels
%  divided by the number of solid particles.
%  - row 6: Number of connected fluid, i.e. void, components.
[domain_data, measures] = run_cam(nx, num_steps,porosity,...
    jump_parameter,output_rate=frame_rate, print_results=output_results,...
    print_measures=output_measures,print_random_seed=output_random_seed,...
    random_seed=num_random_seed);

measures

%% Call the print function to illustrate CAM.
print(nx, domain_data, frame_rate, num_steps)
