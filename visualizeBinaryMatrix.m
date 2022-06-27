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