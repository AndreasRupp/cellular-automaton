function [ outputData ] = run_cam( ...
    nx, ny, numSteps, porosity, jump_param, output_rate )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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

command = strcat(file_name, '(numSteps, porosity, jump_param, ', ...
    'output_rate, zeros(nx * ny, numSteps + 1))');
outputData = eval(command);

cd(current_folder)

end  % function run_cam