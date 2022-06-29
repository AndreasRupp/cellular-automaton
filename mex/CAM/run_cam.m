% TODO: Joona and Simon!

%  ========================================================================
%> @file  matlab_example.m
%>
%> @brief Brief example illustration how to run the cellular automaton,
%         which is implemented in C++, in MATLAB.
%  ========================================================================
%>
%> @brief Brief example illustration how to run the cellular automaton,
%         which is implemented in C++, in MATLAB.
%>
%> TODO: Joona and Simon: Fill this.
%> 
%> The main loop repeatedly executes four steps until the parameter
%> <code>problemData.isFinished</code> becomes <code>true</code>.
%> These four steps are:
%>
%>  1. preprocessStep()
%>  2. solveStep()
%>  3. postprocessStep()
%>  4. outputStep()
%> 
%> This routine is executed second in each loop iteration.
%> It assembles the global system and solves it using the backslash-operator.
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%>
%> @retval problemData  The input struct enriched with solution data at
%>                      the next step. @f$[\text{struct}]@f$
%>
%> This file is part of the GitHub repository
%>   https://github.com/AndreasRupp/cellular-automaton
%> Copyright and license conditions can be found there.

function [ outputData, measures ] = run_cam( ...
    nx, ny, num_steps, porosity, jump_parameter, options )

arguments
    nx int32
    ny int32
    num_steps int32
    porosity double
    jump_parameter double
    options.output_rate int32 = 1
    options.print_results int8  = 1
    options.print_measures int8 = 1
    options.print_random_seed int8 = 0
    options.random_seed int64 = 0
end

current_folder = pwd;

[path,~,~] = fileparts(mfilename('fullpath'));
path = strcat(path, '/../..');
cd(path)

if not(isfolder('build'))
    mkdir('build')
end  % not is folder
cd build

file_name = strcat('m_run_cam_', string(nx), '_', string(ny));
if ~isfile(strcat(file_name, '.mexa64'))
    fid  = fopen('../mex/CAM/run_cam.cxx.in','r');
    text = fread(fid,'*char')';
    fclose(fid);

    text = strrep(text, 'NX_MATLAB_VAL', string(nx));
    text = strrep(text, 'NY_MATLAB_VAL', string(ny));
    
    fid  = fopen(strcat(file_name, '.cxx'),'w');
    fprintf(fid,'%s',text);
    fclose(fid);
    
    mex(strcat(file_name, '.cxx'), '-I../include', ...
        'COMPFLAGS=$COMPFLAGS -O3')
end  % not exists file

n_outputs = 0;
if output_rate ~= 0
    n_outputs = 1 + floor(num_steps / output_rate);
end  % if output rate not 0


if options.print_results
    results_matrix = zeros(nx * ny, n_outputs);
else
    results_matrix = [];
end  % if print_results

if options.print_measures
    measures_matrix = zeros(6, n_outputs);
else
    measures_matrix = [];
end  % if print measures

command = strcat(file_name, '(num_steps, porosity, jump_parameter, ', ...
    'options.output_rate, results_matrix, measures_matrix, ', ...
    'options.print_random_seed, options.random_seed)');
[outputData, measures] = eval(command);

cd(current_folder)

end  % function run_cam