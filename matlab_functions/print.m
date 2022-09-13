%  ========================================================================
%> @file  matlab_functions/print.m
%>
%> @brief Use print_update function to move between time steps.
%  ========================================================================
%>
%> @brief Use print_update function to move between time steps.
%>
%> @param nx        The size of a row for each dimension of the domain.
%> @param data      A matrix with nx(1)*nx(2)*...*nx(n) rows and a column
%>                  for every step with output. The entries of domain_data
%>                  are either 0, corresponding to void cells, i.e.
%>                  non-solid cells, or positive integers corresponding to
%>                  solid cells.
%> @param rate      An integer indicating after how many steps an output
%>                  shall be returned.
%> @param steps     The number of iterations of the CAM.
%>
%>  For every step with output data, we create a figure of the domain.
%>  Press one (1) to move forward and two (2) to move backward in figures.
%>  Press any other button and the program ends.
%>  Remember to press enter after every option to complete the action.
%> 
%> This file is part of the GitHub repository
%>   https://github.com/AndreasRupp/cellular-automaton
%> Copyright and license conditions can be found there.

%% This routine controls the printing for time steps.
function print(nx, data, rate, steps)
    if rate == 0
        return
    end  % if rate == 0

    i = 0;
    while true
        title("Time step: " + i)
        print_update(nx, data(:,i+1));
        event= input(['Press "1" (next) or "2" (previous) to ' ...
                'move between time steps.']);
        switch event
            case 1
                i = mod(i + 1, steps);
            case 2
                i = mod(i - 1, steps);
            otherwise
                return
        end  % switch event
    end  % while true
end  % function print

%% This routine visualizes the given matrix.
%> @brief Visualize output of cellular_automaton in 3D, 2D or 1D for single
%>        time step.
%>
%> @param nx        The vector containing array sizes for each dimension.
%> @param data      The vector containing data for the current time step.
%>
%> We create a figure of the domain with solid pixels as black cubes and
%> void pixels as empty space.
%>
%> This file is part of the GitHub repository
%>   https://github.com/AndreasRupp/cellular-automaton
%> Copyright and license conditions can be found there.
function print_update(nx, data)
    [~, dim] = size(nx);
    if dim == 1
        nx(2) = 1;
        dim = 2;
    end

    x = 0:nx(1);
    y = 0:nx(2);
    if dim == 3
        z = 0:nx(3);
    end

    cla

    if dim == 2
        % Creates gridlines for the domain
        [X1, Y1] = meshgrid(x,y);
        X1 = permute(X1,[2 1]); Y1 = permute(Y1,[2 1]);
        X1(end+1,:) = NaN; Y1(end+1,:) = NaN;
        [X2, Y2] = meshgrid(x,y);
        X2(end+1,:) = NaN; Y2(end+1,:) = NaN;
        h = line([X1(:);X2(:)], [Y1(:);Y2(:)]);
        set(h, 'Color',[0.5 0.5 1], 'LineWidth',1, 'LineStyle','-');

        % Creates faces of the square (particle)
        cube_xy = [0 1 1 0 0; 0 0 1 1 0];
        x = repmat(mod(0:size(data)-1,nx(1))', 5);
        x = x((data ~= 0), :);
        y = repmat(mod(floor((0:size(data)-1)./nx(1)),nx(2))', 5);
        y = y((data ~= 0), :);
        patch('XData', (x + cube_xy(1,:))', 'YData', (y + cube_xy(2,:))')
    end

    if dim == 3
        % Creates gridlines for the domain
        [X1, Y1, Z1] = meshgrid(x,y,z);
        X1 = permute(X1,[2 1 3]); Y1 = permute(Y1,[2 1 3]);
        Z1 = permute(Z1,[2 1 3]);
        X1(end+1,:,:) = NaN; Y1(end+1,:,:) = NaN; Z1(end+1,:,:) = NaN;
        [X2, Y2, Z2] = meshgrid(x,y,z);
        X2(end+1,:,:) = NaN; Y2(end+1,:,:) = NaN; Z2(end+1,:,:) = NaN;
        [X3, Y3, Z3] = meshgrid(x,y,z);
        X3 = permute(X3,[3 1 2]); Y3 = permute(Y3,[3 1 2]);
        Z3 = permute(Z3,[3 1 2]);
        X3(end+1,:,:) = NaN; Y3(end+1,:,:) = NaN; Z3(end+1,:,:) = NaN;

        h = line([X1(:);X2(:);X3(:)], [Y1(:);Y2(:);Y3(:)], ...
            [Z1(:);Z2(:);Z3(:)]);
        set(h, 'Color',[0.5 0.5 1], 'LineWidth',1, 'LineStyle','-');

        % Creates faces of the cube (particle)
        x = repmat(mod(0:size(data)-1,nx(1))', 5);
        x = x((data ~= 0), :).';
        y = repmat(mod(floor((0:size(data)-1)./nx(1)),nx(2))', 5);
        y = y((data ~= 0), :).';
        z = repmat(mod(floor((0:size(data)-1)./(nx(1)*nx(2))),nx(3))', 5);
        z = z((data ~= 0), :).';

        cube_xy = [0 1 1 0 0; 0 0 1 1 0; 0 0 0 0 0];
        patch('XData', x + cube_xy(1,:)', 'YData', y + cube_xy(2,:)', ...
            'ZData', z + cube_xy(3,:)')
        patch('XData', x + cube_xy(1,:)', 'YData', y + cube_xy(2,:)', ...
            'ZData', z + cube_xy(3,:)' + 1)
        cube_xz = [0 1 1 0 0; 0 0 0 0 0; 0 0 1 1 0];
        patch('XData', x + cube_xz(1,:)', 'YData', y + cube_xz(2,:)', ...
            'ZData', z + cube_xz(3,:)')
        patch('XData', x + cube_xz(1,:)', 'YData', y + cube_xz(2,:)'+1, ...
            'ZData', z + cube_xz(3,:)');
        cube_yz = [0 0 0 0 0; 0 1 1 0 0; 0 0 1 1 0];
        patch('XData', x + cube_yz(1,:)', 'YData', y + cube_yz(2,:)', ...
            'ZData', z + cube_yz(3,:)')
        patch('XData', x + cube_yz(1,:)'+1, 'YData', y + cube_yz(2,:)', ...
            'ZData', z + cube_yz(3,:)')
    end

    axis equal
    axis off
    view;
    axis vis3d
    camproj perspective, rotate3d on
end
