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

%% Define the arguments of the run_cam function.
%  This example illustrates how to run the run_cam function, see the file
%  mex/CAM/run_cam.m. It runs a cellular automaton method on a grid of size
%  'nx(1) * nx(2) * ... * nx(n)' for 'num_steps' steps. Moreover, we set
%  the porosity (i.e.the percentage of void space---not occupied by solid)
%  and the jump parameter (telling how far individual cells are allowd to
%  jump) as mandatory parameters.
nx             = [10 10 10];
num_steps      = 3;
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
num_random_seed    = 520646584;

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

%  For every step with output data, we create an image file visualizing the
%  domain with solid pixels in black and void pixels in white. The
%  images are saved as output/fig.i.png. 
% for i = 1 : num_steps + 1
%     if (frame_rate ~= 0 && mod(i , frame_rate) == 0)
%         %  The domain vector at the given step is reshaped to a matrix of
%         %  size [nx(1) nx(2) ... nx(n)].
%         domain_geo = reshape(domain_data(:,floor(i / frame_rate)), ...
%             [nx(1) nx(2)]);
%         visualize_binary_matrix(domain_geo, strcat(strcat('output/fig.',...
%             num2str(floor(i / frame_rate) -1)),'.png'))
%     end
% end

%  For every step with output data, we create a figure of the
%  domain with solid pixels as black cubes and void pixels as empty space.
%  Press one (1) to move forward and two (2) to move backward in figures.
%  Press any other button and the program ends.
%  Remember to press enter after every option to complete the action.

i = 0;

if (frame_rate ~= 0)
    while 1  
        eventdata= input(['Press "1" (next) or "2" (previous) to move ' ...
           'between time steps.']);
        switch eventdata
            case 1
                if i > num_steps
                    i = 0;
                end
                i = i + 1;
                plot_3D_with_grid(nx, domain_data(:,i));
            case 2
                i = i - 1;
                if i < 1
                    i = num_steps + 1;
                end
                plot_3D_with_grid(nx, domain_data(:,i));
            otherwise
                return
        end
     end
end

%% This routine visualizes the given matrix by a black and white image.
%
%> @brief Visualize output of cellular_automaton, in .png files.
%> 
%> @param geo       A matrix with entries of 0 (void pixels) and
%>                  non-zero integers, corresponding to solid pixels.
%> @param fig_path  The name of the file where the image is written to.
%>
%>  The function resizes the matrix geo, s.t. the created image is at least
%>  of size 500 pixels per direction. In the created image, black pixels 
%>  correspond to the non-zero entries of geo, while white pixels
%>  correspond to the zero entries of geo.
% function visualize_binary_matrix(geo, fig_path)
%     if not(isfolder('output'))
%         mkdir('output')
%     end  
%     resize_factor = max(1, min(floor(1000/size(geo,1)), ...
%         floor(1000/size(geo,2))));
%     geo = repelem(geo, resize_factor, resize_factor);
%     imwrite(~geo, fig_path);
% end

%%
function plot_3D_with_grid(nx, data)

    NX = nx;
    [xdim, dim] = size(NX);
    v = data;
    x = 0:NX(1);
    if dim == 2 || dim == 3
        y = 0:NX(2);
    end
    if dim == 3
        z = 0:NX(3);
    end
    cla
    
    % Creates gridlines for the domain
    if dim == 3
        [X1, Y1, Z1] = meshgrid(x,y,z);
        X1 = permute(X1,[2 1 3]); Y1 = permute(Y1,[2 1 3]); Z1 = ...
            permute(Z1,[2 1 3]);
        X1(end+1,:,:) = NaN; Y1(end+1,:,:) = NaN; Z1(end+1,:,:) = NaN;
        [X2, Y2, Z2] = meshgrid(x,y([1 end]),z);
        X2(end+1,:,:) = NaN; Y2(end+1,:,:) = NaN; Z2(end+1,:,:) = NaN;
        [X3, Y3, Z3] = meshgrid(x,y,z([1 end]));
        X3 = permute(X3,[3 1 2]); Y3 = permute(Y3,[3 1 2]); Z3 = ...
            permute(Z3,[3 1 2]);
        X3(end+1,:,:) = NaN; Y3(end+1,:,:) = NaN; Z3(end+1,:,:) = NaN;
        
        h = line([X1(:);X2(:);X3(:)], [Y1(:);Y2(:);Y3(:)], ...
            [Z1(:);Z2(:);Z3(:)]);
        set(h, 'Color',[0.5 0.5 1], 'LineWidth',1, 'LineStyle','-');
    
        % Creates faces of the cube (particle)
        for i = 1:size(v)
            if v(i) ~= 0
                xi = mod(i,NX(1));
                yi = mod(floor(i/NX(1)),NX(2));
                zi = mod(floor(i/(NX(1)*NX(2))),NX(3));
                cube_xy = [0 1 1 0 0;
                           0 0 1 1 0;
                           0 0 0 0 0];
                patch('XData', xi + cube_xy(1,:), 'YData', yi + ...
                    cube_xy(2,:), 'ZData', zi + cube_xy(3,:))
                patch('XData', xi + cube_xy(1,:), 'YData', yi + ...
                    cube_xy(2,:), 'ZData', zi + cube_xy(3,:) + 1)
                
                cube_xz = [0 1 1 0 0;
                           0 0 0 0 0;
                           0 0 1 1 0];
                patch('XData', xi + cube_xz(1,:), 'YData', yi + ...
                    cube_xz(2,:), 'ZData', zi + cube_xz(3,:))
                patch('XData', xi + cube_xz(1,:), 'YData', yi + ...
                    cube_xz(2,:) + 1, 'ZData', zi + cube_xz(3,:))
                
                cube_yz = [0 0 0 0 0;
                           0 1 1 0 0;
                           0 0 1 1 0];
                patch('XData', xi + cube_yz(1,:), 'YData', yi + ...
                    cube_yz(2,:), 'ZData', zi + cube_yz(3,:))
                patch('XData', xi + cube_yz(1,:) + 1, 'YData', yi + ...
                    cube_yz(2,:), 'ZData', zi + cube_yz(3,:))
            end
        end
    end

    if dim == 2
        [X1, Y1] = meshgrid(x,y);
        X1 = permute(X1,[2 1 3]); Y1 = permute(Y1,[2 1 3]);
        X1(end+1,:,:) = NaN; Y1(end+1,:,:) = NaN;
        [X2, Y2] = meshgrid(x,y([1 end]));
        X2(end+1,:,:) = NaN; Y2(end+1,:,:) = NaN;
        [X3, Y3] = meshgrid(x,y);
        X3 = permute(X3,[3 1 2]); Y3 = permute(Y3,[3 1 2]);
        X3(end+1,:,:) = NaN; Y3(end+1,:,:) = NaN;
        
        h = line([X1(:);X2(:);X3(:)], [Y1(:);Y2(:);Y3(:)]);
        set(h, 'Color',[0.5 0.5 1], 'LineWidth',1, 'LineStyle','-');
    
        % Creates faces of the square (particle)
        for i = 1:size(v)
            if v(i) ~= 0
                xi = mod(i,NX(1));
                yi = mod(floor(i/NX(1)),NX(2));
                cube_xy = [0 1 1 0 0;
                           0 0 1 1 0];
                patch('XData', xi + cube_xy(1,:), 'YData', yi + ...
                    cube_xy(2,:))
            end
        end
    end

    if dim == 1
        y = [0 1];
        [X1, Y1] = meshgrid(x,y);
        X1 = permute(X1,[2 1 3]); Y1 = permute(Y1,[2 1 3]);
        X1(end+1,:,:) = NaN; Y1(end+1,:,:) = NaN;
        [X2, Y2] = meshgrid(x,y([1 end]));
        X2(end+1,:,:) = NaN; Y2(end+1,:,:) = NaN;
        [X3, Y3] = meshgrid(x,y);
        X3 = permute(X3,[3 1 2]); Y3 = permute(Y3,[3 1 2]);
        X3(end+1,:,:) = NaN; Y3(end+1,:,:) = NaN;
        
        h = line([X1(:);X2(:);X3(:)], [Y1(:);Y2(:);Y3(:)]);
        set(h, 'Color',[0.5 0.5 1], 'LineWidth',1, 'LineStyle','-');
    
        % Creates faces of the square (particle)
        for i = 1:size(v)
            if v(i) ~= 0
                xi = mod(i,NX(1));
                cube_xy = [0 1 1 0 0;
                           0 0 1 1 0];
                patch('XData', xi + cube_xy(1,:), 'YData', cube_xy(2,:))
            end
        end
    end
    axis equal
    axis off
    view;
    axis vis3d
    camproj perspective, rotate3d on
end