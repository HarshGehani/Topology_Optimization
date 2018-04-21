% Assembling the global stiffness matrix
% using connectivity matrix, 'Conn' (pre-computed)
% 'x' is the matrix of elemental desities
% 'f' & 'k' are elemental stiffness matrices of size (8,8) for
% quadrilateral element (2 DOF per node) (dimensions are (8,8,nElem))
% 'F' & 'K' are global stifness matrices of size () initially all zeros
% 'nX' & 'nY' are the number of cells in the x and y-directions


function K = globalAssembly(E_i,Conn,k,dofPerNode,nX,nY)
    nNodes = (nX + 1)*(nY + 1)*(nZ + 1);
    nElem = nX*nY*nZ;
    N = nNodes*dofPerNode;                           % Total DOFs of the system
    K = zeros(N,N);                                  % Global stiffness matrix
    K0 = zeros(N,N,nElem);                           % Global version of the stiffness matrices for each element (3rd dimension corresponds to a particular element)
    
    
    % Forming the sparse global version of the elemental stiffness matrix
    for e = 1:nElem
%         nodeIndices = zeros(1,8);                   % 2D
        nodeIndices = zeros(1,24);                   % 3D
%         nodeIndices(1:2:numel(nodeIndices)) = 2*Conn(e,:) - 1;
%         nodeIndices(2:2:numel(nodeIndices)) = 2*Conn(e,:);

        % 3D
        nodeIndices(1:3:numel(nodeIndices)) = 3*Conn(e,:) - 2;
        nodeIndices(2:3:numel(nodeIndices)) = 3*Conn(e,:) - 1;
        nodeIndices(3:3:numel(nodeIndices)) = 3*Conn(e,:);
        K0(nodeIndices,nodeIndices,e) = k(:,:,e);
    end
    
    % Forming the stiffness matrix
    i = 1; j = 1;k = 1;
    for e = 1:nElem
        %K = (E_min + x(i,j)*(E_0 - E_min))*K0(:,:,e);
        K = K + (E_i(i,j,k))*K0(:,:,e);
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
