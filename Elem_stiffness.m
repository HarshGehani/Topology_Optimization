function [k0_i_final] = Elem_stiffness(dx, dy)
% % Function to determine element-wise stiffness matrix using coordinate 
% % values of center of every element and the dimensions of the element
    
% % Parameters
    p = 1;  %Penalty
    E_0 = 2e11; %Isotropic material Modulus of elasticity (Steel)
    E_min = 1e-7;   % Arbitrary min value of modulus for avoiding optimization algo breakdown
    nu = 0.3; % Poisson ratio (Steel)

    % syms x_i
    % Heuristic relation between (relative) element density, x_i, and element 
    % Young’s modulus, E_i for modified SIMP approach
    E_i = E_min + (x_i)^p*(E_0 - E_min);      % x_i is in range (0,1)

    % 2D Constitutive matrix
    coeff = 1/((1 + nu)*(1 - 2*nu));
    C0_i = coeff*[(1 - nu), nu, 0;
            nu, (1 - nu), 0;
            0, 0, (1 - 2*nu)/2];

    syms eta zeta

    % Linear Quadilateral element
    N1 = 0.25*(1-zeta)*(1-eta);
    N2 = 0.25*(1+zeta)*(1-eta);
    N3 = 0.25*(1+zeta)*(1+eta);
    N4 = 0.25*(1-zeta)*(1+eta);
    
    grad_N1 = grad(N1); grad_N2 = grad(N2);
    grad_N3 = grad(N3); grad_N4 = grad(N4);

    B = 1/(2*dx*dy)*[grad_N1(1),0,grad_N2(1),0,grad_N3(1),0,grad_N4(1),0;
        0,grad_N1(2),0,grad_N2(2),0,grad_N3(2),0,grad_N4(2);
        grad_N1(2),grad_N1(1),grad_N2(2),grad_N2(1),grad_N3(2),grad_N3(1),grad_N4(2),grad_N4(1)];

    unsolved_k0_i = (B' * C0_i * B);
    
    % Element stiffness matrix
    k0_i = integrateStiffness(unsolved_k0_i);
    k0_i_final = E_i * k0_i;