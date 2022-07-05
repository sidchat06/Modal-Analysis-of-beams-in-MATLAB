clc
clear all
close all
%Assembling Global Stiffness Matrix from Local Stiffness Matrix
L = input('Enter length of the beam: ');
n_ele = input('Enter the number of elements/members: ');
d = input('Enter density of the beam: ');
A = input('Enter area of cross-section of the beam: ');
EI = input('Enter the value of EI (constant for each element): ');
n_node = n_ele + 1;

% 'x' represents position of each node serially
Le = L/n_ele; % Length of each element is the same
K = zeros((n_node * 2));
M = zeros((n_node * 2));

for i = 1 : n_ele
    
    k_local = EI * [ 12/Le^3   6/Le^2  -12/Le^3   6/Le^2;
                     6/Le^2    4/Le    -6/Le^2    2/Le;
                    -12/Le^3  -6/Le^2   12/Le^3  -6/Le^2;
                     6/Le^2    2/Le    -6/Le^2    4/Le ];
    
    m_local = (d*A*Le)/420 * [ 156     -22*Le     54       13*Le;
                              -22*Le    4*Le^2   -13*Le   -3*Le^2;
                               54      -13*Le     156      22*Le;
                               13*Le   -3*Le^2    22*Le    4*Le^2 ];
                 
    K(2*i-1:2*i+2,2*i-1:2*i+2) = K(2*i-1:2*i+2,2*i-1:2*i+2) + k_local;
    M(2*i-1:2*i+2,2*i-1:2*i+2) = M(2*i-1:2*i+2,2*i-1:2*i+2) + m_local;
    
end


% FURTHER CALCULATIONS IS ONLY FOR A SIMPLY SUPPORTED BEAM


% ENTERING BOUNDARY CONDITIONS

K_22 = K([2:2*n_node-2,2*n_node:end],[2:2*n_node-2,2*n_node:end]);
M_22 = M([2:2*n_node-2,2*n_node:end],[2:2*n_node-2,2*n_node:end]);



disp(' ');
K_22
disp(' ');
M_22

% CALCULATING EIGENVALUES AND EIGENVECTORS

[evec, eval] = eig(K_22,M_22);
 
for i = 1 : n_ele
     
    eval(i,i) = sqrt(eval(i,i));
     
end

eval = diag(eval)
evec


% 'eval_transverse' IS THE EIGENVALUE OF EACH DISPLACEMENT (NOT ROTATION) BY EACH DOF
% 'evec_transverse' IS THE EIGENVECTOR OF EACH EIGENVALUE OF EACH DISPLACEMENT (NOT ROTATION) BY EACH DOF

eval_transverse = eval([2:2:(2*n_node - 2)]);
evec_transverse = evec(2:2:(2*n_node - 2),:);

eval_transverse
evec_transverse

% PLOTTING MODE SHAPES

plot(evec_transverse(:,1), 'r', 'linewidth',2)
hold on
plot(evec_transverse(:,2), 'b', 'linewidth',2)
plot(evec_transverse(:,3), 'g', 'linewidth',2)
plot(evec_transverse(:,4), 'c', 'linewidth',2)

legend({'Mode 1', 'Mode 2', 'Mode 3', 'Mode 4'}, 'Location','northeast')

hold off
