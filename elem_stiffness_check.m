nu = 0.3;

k1 = -(6*nu - 4)/9;     k2 = 1/12;      k3 = -1/9;
k4 = -(4*nu - 1)/12;    k5 = (4*nu - 1)/12;     k6 = 1/18;
k7 = 1/24;      k8 = -1/12;     k9 = (6*nu - 5)/36;
k10 = -(4*nu - 1)/24;   k11 = -1/24;    k12 = (4*nu - 1)/24;
k13 = (3*nu - 1)/18;    k14 = (3*nu - 2)/18;

K1_mat = [k1,k2,k2,k3,k5,k5;
    k2,k1,k2,k4,k6,k7;
    k2,k2,k1,k4,k7,k6;
    k3,k4,k4,k1,k8,k8;
    k5,k6,k7,k8,k1,k2;
    k5,k7,k6,k8,k2,k1];

K2_mat = [k9,k8,k12,k6,k4,k7;
    k8,k9,k12,k5,k3,k5;
    k10,k10,k13,k7,k4,k6;
    k6,k5,k11,k9,k2,k10;
    k4,k3,k5,k2,k9,k12;
    k11,k4,k6,k12,k10,k13];

K3_mat = [k6,k7,k4,k9,k12,k8;
    k7,k6,k4,k10,k13,10;
    k5,k5,k3,k8,k12,k9;
    k9,k10,k2,k6,k11,k5;
    k12,k13,k10,k11,k6,k4;
    k2,k12,k9,k4,k5,k3];

K4_mat = [k14,k11,k11,k13,k10,k10;
    k11,k14,k11,k12,k9,k8;
    k11,k11,k14,k12,k8,k9;
    k13,k12,k12,k14,k7,k7;
    k10,k9,k8,k7,k14,k11;
    k10,k8,k9,k7,k11,k14];

K5_mat = [k1,k2,k8,k3,k5,k4;
    k2,k1,k8,k4,k6,k11;
    k8,k8,k1,k5,k11,k6;
    k3,k4,k5,k1,k8,k2;
    k5,k6,k11,k8,k1,k8;
    k4,k11,k6,k2,k8,k1];

K6_mat =[k14,k11,k7,k13,k10,k12;
    k11,k14,k7,k12,k9,k2;
    k7,k7,k14,k10,k2,k9;
    k13,k12,k10,k14,k7,k11;
    k10,k9,k2,k7,k14,k7;
    k12,k2,k9,k11,k7,k14];

coeff = 1/((1+nu)*(1-2*nu));
K = [K1_mat,K2_mat,K3_mat,K4_mat;
    K2_mat',K5_mat,K6_mat,K4_mat';
    K3_mat',K6_mat,K5_mat',K2_mat';
    K4_mat,K3_mat,K2_mat,K1_mat'];

