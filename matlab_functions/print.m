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
    i = 0;
    
    if (rate ~= 0)
        while 1  
            eventdata= input(['Press "1" (next) or "2" (previous) to ' ...
                'move between time steps.']);
            switch eventdata
                case 1
                    if i > steps
                        i = 0;
                    end
                    i = i + 1;
                    print_update(nx, data(:,i));
                case 2
                    i = i - 1;
                    if i < 1
                        i = steps + 1;
                    end
                    print_update(nx, data(:,i));
                otherwise
                    return
            end
            title("Time step: " + i)
         end
    end
end