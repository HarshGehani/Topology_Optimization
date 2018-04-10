function []=viewMatrix(densityMap)
    % function inputs a matrix(density structure) and views it as an image
    
    % check for 2D or 3D structure:
    [~,~,c] = size(densityMap,3);
    
    if c==1:
    % view 2D structure:
    
        % set density values between 0 & 100:
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
        
    else
    % view 3D structure:
        
    end
end
