function [k0_i_final] = Elem_stiffness_3D(dx, dy, dz)
% % Function to determine element-wise stiffness matrix using coordinate 
% % values of center of every element and the dimensions of the element
    
% % Parameters
    p = 1;  %Penalty coefficient
    E_0 = 1; %Isotropic material Modulus of elasticity (Steel)
    % Arbitrary min value of modulus for avoiding optimization algo breakdown
    E_min = 1e-7;   
    nu = 0.3; % Poisson ratio (Steel)

    % syms x_i
    % Heuristic relation between (relative) element density, x_i, and element 
    % Young’s modulus, E_i for modified SIMP approach
    E_i = E_min + ((x_i)^p)*(E_0 - E_min);      % x_i is in range (0,1)

    % 3D Constitutive matrix
    coeff = 1/((1 + nu)*(1 - 2*nu));
    C0_i = coeff*[(1 - nu), nu, nu, 0, 0, 0;
            nu, (1 - nu), nu, 0, 0, 0;
            nu, nu, (1 - nu), 0, 0, 0;
            0, 0, 0, (1 - 2*nu)/2, 0, 0;
            0, 0, 0, 0, (1 - 2*nu)/2, 0;
            0, 0, 0, 0, 0, (1 - 2*nu)/2];
    
    % Linear hexahedral element
    syms eta1 eta2 eta3
    
    N1 = 0.125*(1-eta1)*(1-eta2)*(1-eta3);
    N2 = 0.125*(1+eta1)*(1-eta2)*(1-eta3);
    N3 = 0.125*(1+eta1)*(1+eta2)*(1-eta3);
    N4 = 0.125*(1-eta1)*(1+eta2)*(1-eta3);
    N5 = 0.125*(1-eta1)*(1-eta2)*(1+eta3);
    N6 = 0.125*(1+eta1)*(1-eta2)*(1+eta3);
    N7 = 0.125*(1+eta1)*(1+eta2)*(1+eta3);
    N8 = 0.125*(1-eta1)*(1+eta2)*(1+eta3);
   
    grad_N1 = gradient(N1); grad_N2 = gradient(N2);
    grad_N3 = gradient(N3); grad_N4 = gradient(N4);
    grad_N5 = gradient(N5); grad_N6 = gradient(N6);
    grad_N7 = gradient(N7); grad_N8 = gradient(N8);

    B = 1/(dx*dy*dz)*...
        [grad_N1(1),0,0,grad_N2(1),0,0,grad_N3(1),0,0,grad_N4(1),0,0,...
        grad_N5(1),0,0,grad_N6(1),0,0,grad_N7(1),0,0,grad_N8(1),0,0;
        0,grad_N1(2),0,0,grad_N2(2),0,0,grad_N3(2),0,0,grad_N4(2),0,...
        0,grad_N5(2),0,0,grad_N6(2),0,0,grad_N7(2),0,0,grad_N8(2),0;
        0,0,grad_N1(2),0,0,grad_N2(2),0,0,grad_N3(2),0,0,grad_N4(2),...
        0,0,grad_N5(2),0,0,grad_N6(2),0,0,grad_N7(2),0,0,grad_N8(2);
        grad_N1(2),grad_N1(1),0,grad_N2(2),grad_N2(1),0,grad_N3(2),grad_N3(1),0,grad_N4(2),grad_N4(1),0,...
        grad_N5(2),grad_N5(1),0,grad_N6(2),grad_N6(1),0,grad_N7(2),grad_N7(1),0,grad_N8(2),grad_N8(1),0;
        0,grad_N1(3),grad_N1(2),0,grad_N2(3),grad_N2(2),0,grad_N3(3),grad_N3(2),0,grad_N4(3),grad_N4(2),...
        0,grad_N5(3),grad_N5(2),0,grad_N6(3),grad_N6(2),0,grad_N7(3),grad_N7(2),0,grad_N8(3),grad_N8(2);
        grad_N1(3),0,grad_N1(1),grad_N2(3),0,grad_N2(1),grad_N3(3),0,grad_N3(1),grad_N4(3),0,grad_N4(1),...
        grad_N5(3),0,grad_N5(1),grad_N6(3),0,grad_N6(1),grad_N7(3),0,grad_N7(1),grad_N8(3),0,grad_N8(1)];

    % Element stiffness matrix
    unsolved_k0_i = (B' * C0_i * B);
    k0_i = integrateStiffness(unsolved_k0_i,3);
    k0_i_final = E_i * k0_i;
    
end