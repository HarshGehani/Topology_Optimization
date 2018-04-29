% Assembling the global stiffness matrix
% using connectivity matrix, 'Conn' (pre-computed)
% 'x' is the matrix of elemental desities
% 'f' & 'k' are elemental stiffness matrices of size (8,8) for
% quadrilateral element (2 DOF per node) (dimensions are (8,8,nElem))
% 'F' & 'K' are global stifness matrices of size () initially all zeros
% 'nX' & 'nY' are the number of cells in the x and y-directions


function [K,K0] = globalAssembly(E_i,Conn,ke,dofPerNode,nX,nY,nZ)
    nNodes = (nX + 1)*(nY + 1);
    if nZ ~= 1
        nNodes = nNodes*(nZ + 1);
    end
    nElem = nX*nY*nZ;
    N = nNodes*dofPerNode;                           % Total DOFs of the system
%     K = zeros(N,N);                                  % Global stiffness matrix
    K = sparse(N,N);
    K0 = sparse(N,N);
    %K0 = cell(1,nElem);                             % Cell array; Each cell is an element to store sparse stiffness matrix
    
%     K0 = zeros(N,N,nElem);                           % Global version of the stiffness matrices for each element (3rd dimension corresponds to a particular element)
   
%     % Using sparse
%     K0 = sparse(K0);
%     K = sparse(K);
    
    % Forming the sparse global version of the elemental stiffness matrix
    element = struct();
    for e = 1:nElem
        if nZ == 1        
    %       2D
            nodeIndices = zeros(1,8);                   % 2D
            nodeIndices(1:2:numel(nodeIndices)) = 2*Conn(e,:) - 1;
            nodeIndices(2:2:numel(nodeIndices)) = 2*Conn(e,:);
%             K0(nodeIndices,nodeIndices,e) = k(:,:,e);

            % Sparse Matrix
            elementStiff = sparse(N,N);
            elementStiff(nodeIndices,nodeIndices) = ke(:,:,e);
            element(e).stiffness = elementStiff;
        else    
            % 3D
            nodeIndices = zeros(1,24);                   % 3D
            nodeIndices(1:3:numel(nodeIndices)) = 3*Conn(e,:) - 2;
            nodeIndices(2:3:numel(nodeIndices)) = 3*Conn(e,:) - 1;
            nodeIndices(3:3:numel(nodeIndices)) = 3*Conn(e,:);
%             K0(nodeIndices,nodeIndices,e) = k(:,:,e);
            

            % Sparse Matrix
            elementStiff = sparse(N,N);
            elementStiff(nodeIndices,nodeIndices) = ke(:,:,e);
            element(e).stiffness = elementStiff;
        end
    end
    
    % Forming the stiffness matrix
    i = 1; j = 1;k = 1;
    for e = 1:nElem
    
%         K = K + (E_i(i,j,k))*K0(:,:,e);
        K = K + (E_i(i,j,k))*element(e).stiffness;
        K0 = K0 + element(e).stiffness;
        j = j + 1;
        if j > nX
            j = 1;
            i = i + 1;
        end
        if i > nY
           k = k + 1;
           i = 1;
        end
    end
    
%     for e = 1:N
%        for i = 1:n
%            for j = 1:n
%                K(Conn(e,i),Conn(e,j)) = K(Conn(e,i),Conn(e,j)) + k(i,j);
%            end
%            F(Conn(e,i)) = F(Conn(e,i)) + f(i);
%        end
%     end
    
end
