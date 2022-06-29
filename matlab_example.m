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



%% First step is including the directory 'mex/CAM' to the MATLAB path.
%  This allows to run all MATLAB functions of CAM. These MATLAB functions
%  take care of compiling and running respective C++ functions, in which
%  the actual work is done.
addpath('mex/CAM')

%% Define the arguments of the run_cam function.
%  This example illustrates how to run the run_cam function, see the file
%  mex/CAM/run_cam.m. It runs a celular automaton method on a grid of size
%  'nx * ny' for 'num_steps' steps. Moreover, we set the porosity (i.e.the
%  percentage of void space---not occupied by solid), the jump parameter
%  (telling how far individual cells are allowd to jump), and an output
%  rate (saying after how many steps an output shall be returned).
nx             = 10;
ny             = 10;
num_steps      = 50;
porosity       = 0.9;
jump_parameter = 1;
output_rate    = 1;


outputData = run_cam(nx, ny, num_steps,porosity,jump_parameter,output_rate);

if ~exist('output', 'dir')
    mkdir output;
end

numSolidPixels = sum(outputData(:,1)>0);

for i = 1 : num_steps + 1
    if sum(outputData(:,i)>0) ~= numSolidPixels
        error('Error. Number of solid pixels changed.')
    end
    outputPic = reshape(outputData(:,i), [nx ny]);
    visualizeBinaryMatrix(outputPic, strcat(strcat('output/fig.',num2str(i-1)),'.png')) 
end


function visualizeBinaryMatrix(inputMatrix, figPath) 
%     inputMatrix = inputMatrix > 0;
%     [row,col] = find(inputMatrix > 0);
%     outputPixels = [col-1 size(inputMatrix,2)-row];
%     axh = axes('XTick',0:1:10,'YTick',0:1:10); % Assumes square axes
%     grid on
%     axis equal
%     xlim([0,10])
%     ylim([0,10])
%     hold on
%     % Plot rectangles 
%     rh = arrayfun(@(i) rectangle('position',[outputPixels(i,:), 1, 1], 'FaceColor','k','EdgeColor', 'none'), 1:size(outputPixels,1)); 
%     exportgraphics(axh,figPath)
    inputMatrix = repelem(inputMatrix, floor(1000/size(inputMatrix,1)), floor(1000/size(inputMatrix,2)));
    imwrite(~inputMatrix, figPath);
end