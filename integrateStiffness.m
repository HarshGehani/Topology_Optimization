function [stiffnessMatrix] = integrateStiffness(equationMatrix, nVariable)
    % function integrates over the elemental stiffness matrix
    % and returns its numerical value
    
    % Input:
    %   equationMatrix - Matrix with equations of stiffness element
    %   nVariable - number of variables in the equation
    % Ouput:
    %   stiffnessMatrix - matrix with stiffness values for each element
    
    [m,n] = size(equationMatrix);
    
    % integration limits:
%     eta1Min = -1;
%     eta1Max = 1;
%     eta2Min = -1;
%     eta2Max = 1;
%     eta3Min = -1;
%     eta3Max = 1;

    eta1Min = 0;
    eta1Max = 1;
    eta2Min = 0;
    eta2Max = 1;
    eta3Min = 0;
    eta3Max = 1;

    % check for 2D or 3D structure:
    if nVariable==2
    % 2D structure iteration over each cell:
        for i=1:m
            for j=1:n
                
                % create matlabFunction for equation evaluation:
                fun = matlabFunction(equationMatrix(i,j));

                % check function for zero integral:
                if ~isequal(func2str(fun),func2str(@()0.0))
                    stiffnessMatrix(i,j) = 2*integral2(fun,eta1Min,eta1Max,eta2Min,eta2Max);
                else
                    stiffnessMatrix(i,j) = 0;
                end
            end
        end
        
    else
    % 3D structure iteration over each cell:
        for i=1:m
            for j=1:n
                
                % create matlabFunction for equation evaluation:
                fun = matlabFunction(equationMatrix(i,j));

                % check function for zero integral:
                if ~isequal(func2str(fun),func2str(@()0.0));
                    stiffnessMatrix(i,j) = integral3(fun,eta1Min,eta1Max,eta2Min,eta2Max,eta3Min,eta3Max);
                else
                    stiffnessMatrix(i,j) = 0;
                end
            end
        end      
    end
    
end
