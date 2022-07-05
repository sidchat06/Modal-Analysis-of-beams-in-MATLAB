clc
clear all
close all
%Assembling Global Stiffness Matrix from Local Stiffness Matrix
L = input('Enter length of the beam: ');
n_ele = input('Enter the number of elements/members: ');
n_damage = input('Enter the range damaged members [1st and last position of the members]: ');
p = input('Enter the percentage of flexural rigidity damaged: '); 
d = input('Enter density of the beam: ');
A = input('Enter area of cross-section of the beam: ');
EI = input('Enter the value of EI (constant for each element): ');
a = input('Enter the mode shape to be plotted: ');

n_node = n_ele + 1;
EI_damaged  = EI*(1 - 0.01*p);
% 'x' represents position of each node serially
Le = L/n_ele; % Length of each element is the same
K = zeros((n_node * 2));
M = zeros((n_node * 2));
K_dam = zeros((n_node * 2));

for i = 1 : n_ele
    
    if i >= n_damage(1) && i <= n_damage(2)
        k_local_dam = EI_damaged * [ 12/Le^3   6/Le^2  -12/Le^3   6/Le^2;
                                    6/Le^2    4/Le    -6/Le^2    2/Le;
                                   -12/Le^3  -6/Le^2   12/Le^3  -6/Le^2;
                                    6/Le^2    2/Le    -6/Le^2    4/Le ]; 
    else
        k_local_dam = EI * [ 12/Le^3   6/Le^2  -12/Le^3   6/Le^2;
                             6/Le^2    4/Le    -6/Le^2    2/Le;
                            -12/Le^3  -6/Le^2   12/Le^3  -6/Le^2;
                             6/Le^2    2/Le    -6/Le^2    4/Le ];
    end

        k_local = EI * [ 12/Le^3   6/Le^2  -12/Le^3   6/Le^2;
                         6/Le^2    4/Le    -6/Le^2    2/Le;
                        -12/Le^3  -6/Le^2   12/Le^3  -6/Le^2;
                         6/Le^2    2/Le    -6/Le^2    4/Le ];    
                     
        m_local = (d*A*Le)/420 * [ 156     -22*Le     54       13*Le;
                                  -22*Le    4*Le^2   -13*Le   -3*Le^2;
                                   54      -13*Le     156      22*Le;
                                   13*Le   -3*Le^2    22*Le    4*Le^2 ];
    
                 
    K_dam(2*i-1:2*i+2,2*i-1:2*i+2) = K_dam(2*i-1:2*i+2,2*i-1:2*i+2) + k_local_dam;
    K(2*i-1:2*i+2,2*i-1:2*i+2) = K(2*i-1:2*i+2,2*i-1:2*i+2) + k_local;
    M(2*i-1:2*i+2,2*i-1:2*i+2) = M(2*i-1:2*i+2,2*i-1:2*i+2) + m_local;
    
end


% FURTHER CALCULATIONS IS ONLY FOR A CANTILEVER BEAM


% ENTERING BOUNDARY CONDITIONS

K_22_dam = K_dam(3:2*n_node,3:2*n_node);
K_22 = K(3:2*n_node,3:2*n_node);
M_22 = M(3:2*n_node,3:2*n_node);

% CALCULATING EIGENVALUES AND EIGENVECTORS

[evec, eval] = eig(K_22,M_22);
[evec_dam, eval_dam] = eig(K_22_dam,M_22);
 
for i = 1 : n_ele
     
    eval(i,i) = sqrt(eval(i,i));
    eval_dam(i,i) = sqrt(eval_dam(i,i));
     
end


eval = diag(eval);
eval_dam = diag(eval_dam);

% ALTERNATING TO NOT PLOT THE MODE SHAPES AGAINST THE ROTATIONAL DOF

evec_transverse = evec(1:2:(2*n_node - 2),:);
evec_transverse_dam = evec_dam(1:2:(2*n_node - 2),:);

% PLOTTING MODE SHAPES

% plot(evec_transverse(:,1), 'r', 'linewidth',2)
% hold on
% plot(evec_transverse(:,2), 'b', 'linewidth',2)
% plot(evec_transverse(:,3), 'g', 'linewidth',2)
% plot(evec_transverse(:,4), 'c', 'linewidth',2)
% 
% legend({'Mode 1', 'Mode 2', 'Mode 3', 'Mode 4'}, 'Location','northeast')
% 
% hold off

plot(evec_transverse(:,a), 'r', 'linewidth',2)
hold on
plot(evec_transverse_dam(:,a), 'b', 'linewidth',2)
plot(-evec_transverse_dam(:,a), 'b', 'linewidth',2)
legend({'Baseline', 'Damaged'}, 'Location','northeast')
hold off

