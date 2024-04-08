%Script to estimate the reliability of a composite plate by calculating a
%probability of failure determined by a monte carlo simulation. 
close all
clear

%% Variables
%Units: N, m, Pa

E1_mu = 145.3e9;
E2_mu = 8.5e9;
G12_mu = 4.58e9;
v12_mu = 0.31;
Xt_mu = 1932e6;
Yt_mu = 108e6;
S_mu = 132.8e6;

Xc = 1480e6;
Yc = 220e6;

E1_std = 3.28e9;
E2_std = 1.28e9;
v12_std = 0.018;
G12_std = 0.83e9;
Xt_std = 128.3e6;
Yt_std = 8.2e6;
S_std = 6.21e6;

properties_mu = [E1_mu, E2_mu, G12_mu, v12_mu, Xt_mu, Yt_mu, S_mu];
properties_std = [E1_std, E2_std, G12_std, v12_std, Xt_std, Yt_std, S_std];


%% Initial Layup Definition
theta = [0, 90, 45, -45, -45, 45, 90, 0, 0, 90, 45, -45, -45, 45, 90, 0];
num_layers = length(theta);
t = 0.125e-3;
thickness = ones([1,num_layers]) * t;

%% Initial Load Definition
N = 1200*1000; % N/m
angle_applied = 30; % Deg

Nx = N * cos(deg2rad(angle_applied));
Ny = N * sin(deg2rad(angle_applied));

force = [Nx; Ny; 0; 0; 0; 0];

%% CDF Computation
%Define n: the number of points in distribution, and no_std: the number of
%standard deviations considered
n = 1000;
no_std = 5;

R = 1000;
f=0;
[cdfs, x] = gencdfs(properties_mu, properties_std, n, no_std);

props_store = [];
%% Run Sim

for i=1:R
    count = 0;
    fpf = false;

    %Loop through properties until failure is found
    while fpf == false 
        %Material Property Initialisation
        props = rand_properties(cdfs, x);
        E1 = props(1); E2  = props(2); v12 = props(4); G12 = props(3);
        Xt = props(5); Yt = props(6); S = props(7);

        %Property Recording
        props_store = [props_store; props];

        %ABD Matrix Generation and Midplane Strain
        Q_lam = Qlam(E1, E2, v12, G12);
        Q_array = Qarray(theta, Q_lam);
        ABD = ABD_matrix(theta, thickness, Q_array);
        midplane_strain = linsolve(ABD, force);
    
        %Ply strains and stresses
        [strains_glob, strains_princ, stresses_glob,...
                stresses_princ] = ply_strains(midplane_strain,...
                thickness, theta, Q_lam);
           
        %Failure Analysis
        for k=1:size(stresses_princ,2)
            sigma1 = stresses_princ(1,k); sigma2 = stresses_princ(2,k);
            sigma3 = 0; sigma12 = stresses_princ(3,k);
            mof = 1.1;
            fpfff = puck_ff(sigma1, sigma2, sigma3, Xt, Xc, v12, mof, E1);
            fpfiff = puck_iff(sigma2, sigma12, Yt, Yc, S);

            if fpfff == true
                fpf(k) = true;
            elseif fpfiff == true
                fpf(k) = true;
            else
                fpf(k) = false;
            end
        end
        f = f + 1;
        count = count + 1; 
 
    end
    
    %Runtime tracker
    if mod(i,1000) == 0
        disp(i)
    end
    
    %Probability of failure calcs  
    pf = 1/count;
    pfs(i) = pf;
    pf_r(i) = sum(pfs) / i;
    counts(i) = count;


end
    
    
    
%Sum all pf --> divide by R
pfR = sum(pfs) / R

total_sims = sum(counts) 
altpf = f / total_sims
alt = R/total_sims  


%% Plotting
figure(1);
plot(pfs)
ylim([0 0.04])
xlabel('R iterations')
ylabel('PfR')
figure(2);
plot(pf_r)
ylim([0 0.04])
xlabel('R iterations')
ylabel('PfR')

nbins=100;
figure(3);
xtitles = {'E1', 'E2', 'v12', 'G12', 'Xc', 'Yc', 'S'};
for i = 1:7
    propert = props_store(:,i);
    label = xtitles(i);
    if i == 7
        i = 8;
    end
    subplot(3,3,i);
    histogram(propert, nbins)
    xlabel(label)
end

%% Functions
function [props] = rand_properties(cdfs, x)
    %Initialise props array
    rows = size(cdfs);
    props = zeros(1, rows(1));

    %Loop through all properties in array
    for k = 1:numel(props)
        prob = rand;
        ind = find(min(abs(cdfs(k,:) - prob)) == abs(cdfs(k,:) - prob));
        property = x(k,ind);
        props(k) = property;
    end
end

function [cdfs, x] = gencdfs(props_mu, props_std, N, std_considered)
    %std_considered represents the range of values considered over N steps

    %Initialise array with all possible property values linearly distributed
    x = zeros(numel(props_mu), N);
    cdfs = zeros(numel(props_mu), N);
    
    for i = 1:numel(props_mu)
        mu = props_mu(i);
        sigma = props_std(i);
        x(i,:) = linspace(mu - std_considered*sigma,mu + std_considered*sigma,N);
        cdfs(i, :) = cdf('Normal',x(i,:),mu,sigma);
    end
end
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
    % k = numel(theta);
    Q_array = zeros(3,3,8);
    
    for k = 1:numel(theta)
        angle = theta(k);
        Q_array(:,:,k) = Q_matrix(angle, Q_lam);
    end
end

function [ABD] = ABD_matrix(theta, thickness, Q_array)
    tol = 1e-5;
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

    ABD(ABD<tol)=0;
end

function [global_prop] = glob_prop(A, thickness_array)
    Ex = ( A(1,1)*A(2,2) - A(1,2)^2) / (sum(thickness_array)*A(2,2));
    Ey = ( A(1,1)*A(2,2) - A(1,2)^2) / (sum(thickness_array)*A(2,2));
    vxy = A(1,2) / A(2,2);
    vyx = A(1,2) / A(1,1);
    Gxy = A(3,3) * sum(thickness_array);
    global_prop = [Ex, Ey, vxy, vyx, Gxy];
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

function [failure] = puck_ff(sigma1, sigma2, sigma3, Xt, Xc, v12, mof, E1)
    
    % Assume fibre poisson is 0.2, E_f = 225
    vf = 0.2;
    Ef = 225e9;

    %Initialise logical representing whether its failed. 
    failure = false;

    if sigma1 > 0
        R = Xt;
    else
        R = Xc;
    end

    a = (1/R);
    b = sigma1 - (v12 - vf*mof* E1/Ef)*(sigma2+sigma3);
    

    f = a * b;
    if f >= 1
        failure = true;
    end
end

function [failure] = puck_iff(sigma2, sigma12, Yt, Yc, S)
    failure = false;
    p12p = 0.3;
    p12n = 0.2;

    %Mode A
    if sigma2 > 0
        p12p = 0.3;
        a = (sigma12/S)^2;
        b = (1-p12p*Yt/S)^2;
        c = (sigma2/Yt)^2;
        d = p12p * sigma2 / S;

        f = sqrt(a + b * c) + d;

    %Mode B and C 
    elseif sigma2 < 0
        p23n = p12n * sigma23A / sigma12u;
        sigma23a = S / (2 * p12n) * (sqrt(1+2*p12n*Yc/S)-1);
        sigma12c = S* sqrt(1 + 2*p23n);

        if abs(sigma2/sigma12) >= 0 && abs(sigma2/sigma12) <= sigma23a / abs(sigma12c)
            %Mode B
            f = 1/S * (sqrt(sigma12^2 + (p12n * sigma2)^2) + p12n*sigma2);
        else
            %Mode C
            f = ((sigma12 / (2*(1 + p23n*S)))^2 + ...
                (sigma2 / Yc)^2) * Yc / (-sigma2);
        end
    end

    if f >= 1
        failure = true;         
    end
end

