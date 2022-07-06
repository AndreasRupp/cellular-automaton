function [ domain_data, measures ] = run_cam(...
    nx, num_steps, porosity, jump_parameter, options )

arguments
    nx (1,:) int32
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
cd(strcat(path, '/../..'))

if not(isfolder('build'))
    status = system('cmake -version');
    if status ~= 0
        cd (current_folder)
        error(strcat('Your MATLAB version does not support CMAKE. ', ...
            ' Please run ./setup.sh manually.'));
    end
    system(strcat('cmake -E make_directory build && cd build && ', ...
        'CXX=g++-10 cmake .. && cd ..'), '-echo');
end

cd('build')

[ domain_data, measures ] = compile_and_run( ...
    nx, num_steps, porosity, jump_parameter, ...
    output_rate       = options.output_rate, ...
    print_results     = options.print_results, ...
    print_measures    = options.print_measures, ...
    print_random_seed = options.print_random_seed, ...
    random_seed       = options.random_seed);

cd(current_folder);

end