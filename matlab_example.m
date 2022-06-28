% TODO: Joona, and Simon: Please describe what is done here.

addpath('mex/CAM')

numSteps = 50;
porosity = 0.9;
jump_param =1;
output_rate = 1;
nx = 10;
ny = 10;
domainSize = nx * ny;

outputData = run_cam(nx, ny, numSteps,porosity,jump_param,output_rate);

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


function visualizeBinaryMatrix(inputMatrix, figPath) 
    inputMatrix = inputMatrix > 0;
    [row,col] = find(inputMatrix > 0);
    outputPixels = [col-1 size(inputMatrix,2)-row];
    axh = axes('XTick',0:1:10,'YTick',0:1:10); % Assumes square axes
    grid on
    axis equal
    xlim([0,10])
    ylim([0,10])
    hold on
    % Plot rectangles 
    rh = arrayfun(@(i) rectangle('position',[outputPixels(i,:), 1, 1], 'FaceColor','k','EdgeColor', 'none'), 1:size(outputPixels,1)); 
    exportgraphics(axh,figPath)
end