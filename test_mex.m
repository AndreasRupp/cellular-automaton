%mex('mex/CAM/run_cam.cxx', '-Iinclude', 'COMPFLAGS=$COMPFLAGS -O3')

addpath('mex/CAM')

numSteps = 50;
porosity = 0.9;
jump_param =1;
output_rate = 1;
nx = 10;
ny = 10;
domainSize = nx * ny;

outputData = run_cam(nx, ny, numSteps,porosity,jump_param,output_rate,zeros(domainSize,numSteps + 1));

if ~exist('output', 'dir')
    mkdir output;
end

numSolidPixels = sum(outputData(:,1)>0);

for i = 1 : numSteps + 1
    if sum(outputData(:,i)>0) ~= numSolidPixels
        error('Error. Number of solid pixels changed.')
    end
    outputPic = reshape(outputData(:,i), [nx ny]);
    visualizeBinaryMatrix(outputPic, strcat(strcat('output/fig.',num2str(i-1)),'.png')) 
end

