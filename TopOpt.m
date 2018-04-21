% FEM SOLVER
function TopOpt
volfrac = 0.8;

length = 1; breadth = 0.025;             % Dimensions of the beam
nX = 20; nY = 2; nZ = 2;                      % Number of cells in x and y-directions
dx = length/nX; dy = breadth/nY;      % Dimensions of each element
v = (dx*dy*ones(nY,nX));
% nodesPerElement = 4;                    % Considering linear quadrilateral elementq
nodesPerElement = 8;
dofPerNode = 3;
% dofPerNode = 2;
nNodes = (nX + 1)*(nY + 1)*(nZ + 1);
nElem = nX*nY*nZ;                          % Total number of elements

element_state = ones(nElem);            % Whether or not the element has been removed
x_threshold = 0.5;                      % Elements below this threshold are removed          

E_0 = 1;
E_min = 1e-7;
nele = nX*nY;

%x0 = volfrac*ones(nY,nX);
x0 = ones(nY,nX,nZ);

maxloop = 200;    % Maximum number of iterations
tolx = 0.01;   

f=@(x) FEMSolver(x);

disp(f(x0));

%A = -x0;
%B = volfrac;
A= [];
B= [];
Aeq = [];
Beq = [];
LB = zeros(size(x0));
UB = ones(size(x0));

options = optimoptions('fmincon','TolX',tolx, 'MaxIter',maxloop,'Display','iter','Algorithm','sqp');
z = fmincon(f, x0, A, B, Aeq, Beq, LB, UB, @(x) nonlcon(x),options);

disp (z);


% FEM Solver (Objective function)
function f=FEMSolver(x)
% Elemental densities
%x = sym('x',[nX nY]);
%x = ones(nY,nX);
E_i = E_0*ones(nY,nX,nZ) + x*(E_0 - E_min);

% Elemental stiffness matrix
k = zeros(dofPerNode*nodesPerElement,dofPerNode*nodesPerElement,nElem);

% Global Load vector
F = zeros(dofPerNode*nNodes,1);
F(2*(nX + 1)) = -1000;

% Setting up the boundary conditions (Indices of the rows and columns to be deleted from K and F)
boundaryNodeIndices = zeros(1,2*(nY + 1)*(nZ + 1));                 % 3D
% boundaryNodeIndices = zeros(1,2*(nY + 1));                          % 2D

% Forming the connectivity matrix (3D)
Conn = zeros(nElem,nodesPerElement);
i = 1; boundaryCount = 1; z_idx = 0;                % z_idx is basically the z-coordinate
nNodes_xy = (nX + 1)*(nY + 1);
for e = 1:nElem
    % Checking for boundary nodes
    if mod(e,nX) == 1
%         boundaryNodeIndices(boundaryCount:boundaryCount + 3) = [2*i - 1,2*i];
%         boundaryNodeIndices(boundaryCount + 4:boundaryCount + 7) = [2*(i + nX + 1) - 1,2*(i + nX + 1)];
%         boundaryCount = boundaryCount + 8;
    end
    lowerFace = [i + z_idx*nNodes_xy,i + 1 + z_idx*nNodes_xy,i + nX + 2 + z_idx*nNodes_xy,i + nX + 1 + z_idx*nNodes_xy];
    upperFace = [i + (z_idx +  1)*nNodes_xy,i + 1 + (z_idx +  1)*nNodes_xy, i + nX + 2 + (z_idx +  1)*nNodes_xy, i + nX + 1 + (z_idx +  1)*nNodes_xy];
    Conn(e,:) = [lowerFace,upperFace];
    i = i + 1;
    if i == nNodes_xy
       z_idx = z_idx + 1;
       i = 1;
    end
    if mod(i,nX + 1) == 0
        i = i + 1;
    end
end
boundaryNodeIndices = unique(boundaryNodeIndices);


% % Forming the connectivity matrix (2D)
% Conn = zeros(nElem,nodesPerElement);
% i = 1; boundaryCount = 1;
% for e = 1:nElem
%     % Checking for boundary nodes
%     if mod(e,nX) == 1
%         boundaryNodeIndices(boundaryCount:boundaryCount + 1) = [2*i - 1,2*i];
%         boundaryNodeIndices(boundaryCount + 2:boundaryCount + 3) = [2*(i + nX + 1) - 1,2*(i + nX + 1)];
%         boundaryCount = boundaryCount + 4;
%     end
%     Conn(e,:) = [i,i + 1,i + nX + 2,i + nX + 1];
%     i = i + 1;
%     if i == nX + 1
%         i = i + 1;
%     end
% end
% boundaryNodeIndices = unique(boundaryNodeIndices);

% Calculating the elemental stiffness matrices
i = 1; j = 1;
k_element = Elem_stiffness(x(i,j), dx, dy);             % Elemental stiffness matrix
for e = 1:nElem
    centerX = (j - 1)*dx + dx/2.0;
    centerY = (i - 1)*dy + dy/2.0;
   % k(:,:,e) = k(:,:,e) + Elem_stiffness(x(i,j), dx, dy);
    k(:,:,e) = k(:,:,e) + k_element;
    j = j + 1;
    if j > nX
        j = 1;
        i = i + 1;
    end
end

% Assembling the global stiffness matrix
K = globalAssembly(E_i,Conn,k,dofPerNode,nX,nY);

indices = 1:dofPerNode*nNodes;                  % Indices of all the nodes
indices(boundaryNodeIndices) = [];

% Removing the rows + columns that correspond to nDOFs = 0 (Boundary nodes)
K(boundaryNodeIndices,:) = [];                      % Rows
K(:,boundaryNodeIndices) = [];                      % Columns

F_calc = F;
F_calc(boundaryNodeIndices) = [];                      % Rows
%F(:,boundaryNodeIndices) = [];                      % Columns

u = K\F_calc;                                            % Nodal Displacement vector
% u_x = reshape(u(1:2:end),2,6);
% u_y = reshape(u(2:2:end),2,6);

% Re-inserting the boundary nodes into the displacement matrix
U = zeros(dofPerNode*nNodes,1);
U(indices) = u;
U(boundaryNodeIndices) = 0;

%f=F'*u;
f = F'*U;
%viewMatrix(u);
end


% Non-linear constraints (volume fraction)
function [c, ceq] = nonlcon(x)
        
        %c = sum((x'*v)) - volfrac*nX*nX;
        c = sum(sum((dy*dx)*(x))) - volfrac*length*breadth;
        ceq = [];
    end
end
