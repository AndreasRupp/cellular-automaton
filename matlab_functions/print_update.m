%  ========================================================================
%> @file  matlab_functions/print_update.m
%>
%> @brief Visualize output of cellular_automaton in 3D, 2D or 1D for single
%>        time step.
%  ========================================================================
%>
%> @brief Visualize output of cellular_automaton in 3D, 2D or 1D for single
%>        time step.
%>
%> @param nx        The vector containing array sizes for each dimension.
%> @param data      The vector containing data for the current time step.
%>
%>  We create a figure of the domain with solid pixels as black cubes and
%>  void pixels as empty space.
%> 
%> This file is part of the GitHub repository
%>   https://github.com/AndreasRupp/cellular-automaton
%> Copyright and license conditions can be found there.

%% This routine visualizes the given matrix.
function print_update(nx, data)
    [~, dim] = size(nx);
    x = 0:nx(1);
    if dim == 2 || dim == 3
        y = 0:nx(2);
    end
    if dim == 3
        z = 0:nx(3);
    end
    cla
    
    if dim == 3
        % Creates gridlines for the domain
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
        for i = 1:size(data)
            if data(i) ~= 0
                xi = mod((i-1),nx(1));
                yi = mod(floor((i-1)/nx(1)),nx(2));
                zi = mod(floor((i-1)/(nx(1)*nx(2))),nx(3));
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
        % Creates gridlines for the domain
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
        for i = 1:size(data)
            if data(i) ~= 0
                xi = mod((i-1),nx(1));
                yi = mod(floor((i-1)/nx(1)),nx(2));
                cube_xy = [0 1 1 0 0;
                           0 0 1 1 0];
                patch('XData', xi + cube_xy(1,:), 'YData', yi + ...
                    cube_xy(2,:))
            end
        end
    end

    if dim == 1
        % Creates gridlines for the domain
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
        for i = 1:size(data)
            if data(i) ~= 0
                xi = mod((i-1),nx(1));
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