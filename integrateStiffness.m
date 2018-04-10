function [stiffnessMatrix] = integrateStiffness(equationMatrix)
    % function integrates eta and zeta of the elemental stiffness matrix
    % and returns its numerical value
    
    [m,n] = size(equationMatrix);
    
    % integration limits:
    etaMin = -1;
    etaMax = 1;
    zetaMin = -1;
    zetaMax = 1;
    
    % iterate over each cell of elemental stiffness matrix:
    for i=1:m
        for j=1:n
            fun = matlabFunction(equationMatrix(i,j));
            
            % check function for zero integral:
            if ~isequal(func2str(fun),func2str(@()0.0))
                stiffnessMatrix(i,j) = integral2(fun,etaMin,etaMax,zetaMin,zetaMax);
            else
                stiffnessMatrix(i,j) = 0;
            end
        end
    end
    
end