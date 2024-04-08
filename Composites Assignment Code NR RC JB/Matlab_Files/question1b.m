% Script to calculate the strain and stress distributions through the thickness 
% of a given laminate, then plots the strains and stresses in the principal
% coordinate system

close all
clear

%% Laminate Properties
E1 = 145.3e9;
E2 = 8.5e9;
G12 = 4.58e9;
v12 = 0.31;
Xt = 1932e6;
Yt = 108e6;
S = 132.8e6;

E1_std = 3.28e9;
E2_std = 1.28e9;
v12_std = 0.018;
G12_std = 0.83e9;
Xt_std = 128.3e6;
Yt_std = 8.2e6;
S_std = 6.21e6;

%% Layup Properties
theta = [0, 0, 90, 30, 90];

num_layers = length(theta);
t = 0.125e-3;

thickness = ones([1,num_layers]) * t;

%% Input Force
force = [0.2e2; 1.8e4; 0; 18e3; 0; 0];

%% ABD Generation
Q_lam = Qlam(E1, E2, v12, G12);
Q_array = Qarray(theta, Q_lam);
ABD = ABD_matrix(theta, thickness, Q_array);
abd = ABD^-1;

%% Midplane Strain Computation
midplane_strain = linsolve(ABD, force);

%% Ply Analysis
[z_list, angle_at_z] = point_comp(theta, thickness);

[strains_glob, strains_princ, stresses_glob,...
    stresses_princ] = ply_strains(midplane_strain, thickness, theta, Q_lam);

%% Plotting    
z = zk(thickness);
%% Principle Strains
strains = figure(2);
xtitles = {'Strain in 1'; 'Strain in 2'; 'Shear strain in 12'};
for i = 1:3
    subplot(1,3,i);
    plot(strains_princ(i,:), z_list, "-o", 'MarkerSize',3, 'Linewidth',0.8)
    yline(0)
    xline(0)
    yticks(z)
    grid on
    xlabel(xtitles(i))
end
ax=axes(strains,'visible','off'); 
ax.YLabel.Visible='on';
ylabel(ax,'Z (mm) - thickness');
yh1 = get(gca,'ylabel'); 
p = get(yh1,'position'); 
p(1) = 1.3*p(1) ;        
set(yh1,'position',p)  
   
%% Principle Stresses
stresses = figure(4);
xtitles = {'Sigma1'; 'Sigma2'; 'Sigma12'};
for i = 1:3
    subplot(1,3,i);
    plot(stresses_princ(i,:), z_list, "-o", 'MarkerSize',3, 'Linewidth',0.8)
    yline(0)
    xline(0)
    yticks(z)
    grid on
    xlabel(xtitles(i))
end
han=axes(stresses,'visible','off'); 
han.YLabel.Visible='on';
ylabel(han,'Z (mm) - thickness');
yh = get(gca,'ylabel'); 
p = get(yh,'position'); 
p(1) = 1.3*p(1) ;        
set(yh,'position',p) 
    


%% Functions

function [z] = zk(thickness)
    h = abs( sum(thickness)/2 );
    z = -h;
    for k = 1:numel(thickness)
        z = [z z(end)+thickness(k)];

    end
end

function [Q_lam] = Qlam(E1, E2, v12, G12)
    v21 = v12 * E2 / E1;

    Q = 1 - v12*v21;
    Q11 = E1 / Q;
    Q22 = E2 / Q;
    Q12 = v12 * E2 / Q;
    Q66 = G12;

    Q_lam = [Q11 Q12 0; Q12 Q22 0; 0 0 Q66];
end

function [Q] = Q_matrix(theta_angle, Q_lam)
    m = cos(deg2rad(theta_angle));
    n = sin(deg2rad(theta_angle));
    
    Q11 = Q_lam(1,1);
    Q12 = Q_lam(1,2);
    Q22 = Q_lam(2,2);
    Q66 = Q_lam(3,3);

    Qxx = Q11*m^4 + 2 * (Q12 + 2*Q66)*(m^2)*(n^2) + Q22*n^4;
    Qxy = (Q11 + Q22 - 4*Q66)*(m^2)*(n^2) + Q12*(m^4 + n^4);
    Qyy = Q11*n^4 + 2*(Q12 + 2*Q66)*(m^2)*(n^2) + Q22*m^4;
    Qxs = (Q11-Q12-2*Q66)*n*m^3 + (Q12-Q22+2*Q66)*n^3*m;
    Qys = (Q11-Q12-2*Q66)*m*n^3 + (Q12-Q22+2*Q66)*m^3*n;
    Qss = (Q11 + Q22 - 2*Q12 - 2*Q66)*(m^2)*(n^2) + Q66*(n^4 + m^4);
    
    Q = zeros(3);
    Q(1,1) = Qxx;
    Q(1,2) = Qxy;
    Q(2,1) = Qxy;
    Q(2,2) = Qyy;
    Q(3,3) = Qss;
    Q(3,1) = Qxs;
    Q(1,3) = Qxs;
    Q(3,2) = Qys;
    Q(2,3) = Qys;
    
end

function [Q_array] = Qarray(theta, Q_lam)
    num_lam = size(theta);
    Q_array = zeros(3,3,num_lam(2));
    
    for k = 1:numel(theta)
        angle = theta(k);
        Q_array(:,:,k) = Q_matrix(angle, Q_lam);
    end
end

function [ABD] = ABD_matrix(theta, thickness, Q_array)
    A = zeros(3);
    B = zeros(3);
    D = zeros(3);
    z = zk(thickness);

    for k = 1:numel(theta) 
        Q = Q_array(:,:,k);
    
        A = A + Q .* thickness(k);
        B = B + (1/2) * Q .* ( z(k+1)^2 - z(k)^2 );
        D = D + (1/3) * Q .* ( (z(k+1))^3 - (z(k))^3 );
    
    end

    ABD = [A B; B D]; 
end


function [prin_args] = direc_matrix(glob_args, theta)
    m = cos(deg2rad(theta));
    n = sin(deg2rad(theta));
    direc = [m^2, n^2, m*n; n^2, m^2, -m*n; -2*m*n, 2*m*n, (m^2 - n^2)];

    prin_args = direc * glob_args;
end

function [z_list, angle_at_z] = point_comp(theta, thickness)
    angle_at_z=[theta; theta];
    angle_at_z=angle_at_z(:)';
   
    z = zk(thickness);
    z_list=[z; z];
    z_list = z_list(:)';
    z_list = z_list(2:end-1);
end

function [strains_glob, strains_princ, stresses_glob, stresses_princ...
    ] = ply_strains(midplane_strain, thickness, theta, Q_lam)
    
    [z_list, angle_at_z] = point_comp(theta, thickness);

    %Strain Calc
    epsilonxx = midplane_strain(1) + z_list*midplane_strain(4);
    epsilonyy = midplane_strain(2) + z_list*midplane_strain(5);
    gammaxy = midplane_strain(3) + z_list*midplane_strain(6);

    strains_glob = [epsilonxx; epsilonyy; gammaxy];

    strains_princ = zeros(size(strains_glob));

    for i=1:numel(z_list)
        theta = angle_at_z(i);
        strains = strains_glob(:,i);
        strains_princ(:,i) = direc_matrix(strains, theta);
    end
    
    %Stress Calc
    stresses_glob = zeros(size(strains_glob));

    for i = 1:numel(z_list)     
        stresses_glob(:,i) = Q_matrix(angle_at_z(i), Q_lam) * strains_glob(:,i);
    end

    stresses_princ = zeros(size(stresses_glob));
    for k=1:numel(z_list)
        stresses_princ(:,k) = Q_lam * strains_princ(:,k);
    end
  
end

