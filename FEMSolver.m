% FEM SOLVER

length = 1; breadth = 0.25;             % Dimensions of the beam
nX = 6; nY = 2;                      % Number of cells in x and y-directions
dx = length/nX; dy = breadth/nY;        % Dimensions of each element
nodesPerElement = 4;                    % Considering linear quadrilateral elementq
dofPerNode = 2;
nNodes = (nX + 1)*(nY + 1);
nElem = nX*nY;                          % Total number of elements
E_0 = 2e11;
E_min = 1e-7;

% Elemental densities
%x = sym('x',[nX nY]);
x = ones(nX,nY);
E_i = E_0*ones(nX,nY) + x*(E_0 - E_min);

% Elemental stiffness matrix
k = zeros(dofPerNode*nodesPerElement,dofPerNode*nodesPerElement,nElem);

% Global Load vector
F = zeros(dofPerNode*nNodes,1);

% Setting up the boundary conditions (Indices of the rows and columns to be deleted from K and F)
boundaryNodeIndices = zeros(1,2*(nY + 1));

% Forming the connectivity matrix
Conn = zeros(nElem,nodesPerElement);
i = 1; boundaryCount = 1;
for e = 1:nElem
    % Checking for boundary nodes
    if mod(e,nX) == 1
        boundaryNodeIndices(boundaryCount:boundaryCount + 1) = [2*i - 1,2*i];
        boundaryNodeIndices(boundaryCount + 2:boundaryCount + 3) = [2*(i + nX + 1) - 1,2*(i + nX + 1)];
        boundaryCount = boundaryCount + 4;
    end
    Conn(e,:) = [i,i + 1,i + nX + 2,i + nX + 1];
    i = i + 1;
    if i == nX + 1
        i = i + 1;
    end
end
boundaryNodeIndices = unique(boundaryNodeIndices);

% Calculating the elemental stiffness matrices
i = 1; j = 1;
for e = 1:nElem
    centerX = (j - 1)*dx + dx/2.0;
    centerY = (i - 1)*dy + dy/2.0;
    k(:,:,e) = k(:,:,e) + Elem_stifness(centerX, centerY, dx, dy);
    j = j + 1;
    if j > nX
        j = 1;
        i = i + 1;
    end
end

% Assembling the global stiffness matrix
K = globalAssembly(E_i,Conn,k,dofPerNode,nX,nY);

% Removing the rows + columns that correspond to nDOFs = 0 (Boundary nodes)
K(boundaryNodeIndices,:) = [];                      % Rows
K(:,boundaryNodeIndices) = [];                      % Columns

F(boundaryNodeIndices,:) = [];                      % Rows
F(:,boundaryNodeIndices) = [];                      % Columns

u = K\F;                                            % Nodal Displacement vector


