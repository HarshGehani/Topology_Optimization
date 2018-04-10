function []=viewMatrix(densityMap)
    % function inputs a matrix(density structure) and views it as an image
    
    % set density values between 0 & 1:
    densityMap = 100*densityMap;
    
    % resizing image:
    [~,n] = size(densityMap);
    if n < 100
        densityMap2 = imresize(densityMap,100);
    
        % view resized image:
        figure;
        image(densityMap2);
        colorbar;
        title("Optimized structure: rescaled (x 100)");
        ylabel('length');
        xlabel('breadth');
    end
    
    % view original image:
    figure;
    image(densityMap);
    colorbar;
    title("Optimized structure");
    ylabel('length');
    xlabel('breadth');
end