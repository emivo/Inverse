%% Exercise 4 M2 
% Emil Airta

% Load the data
load data/deco03 m0 mn mn2
% Differentiate the ESF
psf_emp = diff(m0);
psf_emp = [0 psf_emp(2:end)]; % shift little bit to left
indices = find(psf_emp > 0.01);
dy = indices(end) - indices(1);

% Build A
A = zeros(256,256);
for k = 1:256
    for j = k:256
        A(j,k) = 1/256; 
    end
end

% Determine the SVD of matrix A
[U1,D1,V1] = svd(A);
svals1   = diag(D1);

% Compute truncated SVD reconstruction
r_alpha = 30;
Dp_alpha = zeros(size(A));
for iii = 1:r_alpha
    Dp_alpha(iii,iii) = 1/svals1(iii);
end
fn = V1*(Dp_alpha.')*(U1.')*mn(:);
fn = fn/sum(fn);

fn2 = V1*(Dp_alpha.')*(U1.')*mn2(:);
fn2 = fn2/sum(fn2);

% Normalization PSF
C = sum(psf_emp);
p = psf_emp/C;

% Theoretical PSF
xvec = linspace(-127, 127, 255);
psf = PSF(xvec, floor(dy/2));

% Symmetrice PSF
psf_sym1 = (p + fliplr(p))/2; % This one is better
psf_sym2 = (p +flipud(p))/2;

%% Plot curves in M2
figure(1)
clf
subplot(1,2,1);
plot(xvec, p, 'b-o', xvec, psf, 'r')
hold on
plot(xvec, fn(2:end), 'k-*')
xlim([-dy*2 dy*2])
legend('Low-noise PSF', 'Theoretical PSF', 'Medium noise tSVD PSF')
subplot(1,2,2);
plot(xvec, p, 'b-o', xvec, psf, 'r')
hold on
plot(xvec, fn2(2:end), 'g-p')
xlim([-dy*2 dy*2])
legend('Low-noise PSF', 'Theoretical PSF', 'High noise tSVD PSF')

figure(2)
clf
subplot(1,2,1);
plot(xvec, psf_sym1, 'b', xvec, psf_sym2, 'k', xvec, psf, 'r')
legend('Low-noise symmetric PSF', 'Low-noise symmetric PSF', 'Theoretical PSF')
xlim([-dy*2 dy*2])

subplot(1,2,2);
plot(xvec, psf, 'r', xvec, fn(2:end), 'b', xvec, fn2(2:end), 'g')
legend('Theoretical PSF', 'Medium noise tSVD PSF', 'High noise tSVD PSF')
xlim([-dy*2 dy*2])

%% M3
% Target f
f_zeros = zeros(128,1);
f_ones = ones(128,1);
f = [f_zeros; f_ones];
% compute new A
% Parameters
n = 256;
delta_x = 1/n;
v = dy/2;

% Calculate convolution matrix
conA2 = convmtx(psf_sym1(indices), n);
A2 = conA2(:, (v+1):(end-v));
[U, D, V] = svd(A2);
svals = diag(D);
% Init L
I = eye(n-1,n);
L = (1/delta_x)*(circshift(I,1,2) - I);

% Init alphas
na = 8;
alphas = zeros(2*na + 1, 1);
len = length(alphas);
for k = 1:2*na
    alphas(k) = 10^(na - k);
end
% Init result variables and m prime to L
f_alphas_0 = zeros(n, len);
f_alphas_n = zeros(n, len);
f_alphas_n2 = zeros(n, len);
rsnes_0 = zeros(len, 1); 
rsnes_n = zeros(len, 1); 
rsnes_n2 = zeros(len, 1); 
m_prime_0 = [m0'; zeros(n - 1, 1)];
m_prime_n = [mn'; zeros(n - 1, 1)];
m_prime_n2 = [mn2'; zeros(n - 1, 1)];
% Compute stacked form
for k = 1:len
    alpha = alphas(k);

    A_prime = [A2; sqrt(alpha)*L];
    
    f_alphas_0(:,k) = A_prime \ m_prime_0;
    f_alphas_n(:,k) = A_prime \ m_prime_n;
    f_alphas_n2(:,k) = A_prime \ m_prime_n2;
    % Relative square norm error
    rsnes_0(k)  = norm(f_alphas_0(:,k) - f)/norm(f) * 100;
    rsnes_n(k)  = norm(f_alphas_n(:,k) - f)/norm(f) * 100;
    rsnes_n2(k)  = norm(f_alphas_n2(:,k) - f)/norm(f) * 100;
end
% Tikhonov
% init result variables
halflen = ceil(len/2);
recs = zeros(n, halflen);
recs_n = zeros(n, halflen);
recs_n2 = zeros(n, halflen);
rsnes0 = zeros(halflen, 1); 
rsnesn = zeros(halflen, 1); 
rsnesn2 = zeros(halflen, 1); 
for k = halflen:len
    index = k - halflen + 1;
    % Regularization parameter
    delta = alphas(k); 

    % Build the matrix Dplus_alpha
    Dplus = zeros(size(D));
    Dplus(1:n,1:n) = diag(svals./(svals.^2+delta));

    % Compute reconstructions
    recs(:, index)    = V*Dplus*(U.')*m0';
    recs_n(:, index)   = V*Dplus*(U.')*mn';
    recs_n2(:, index)   = V*Dplus*(U.')*mn2';
    
    % Relative square norm error
    rsnes0(index)  = norm(recs(:,index) - f)/norm(f) * 100;
    rsnesn(index)  = norm(recs_n(:,index) - f)/norm(f) * 100;
    rsnesn2(index)  = norm(recs_n(:,index) - f)/norm(f) * 100;
end
% Minimums for stacked
[v0,iii] = min(rsnes_0);
[vn,iiin] = min(rsnes_n);
[vn2,iiin2] = min(rsnes_n2);
% Minimums for Tikhonov
[v_0,ind] = min(rsnes0);
[v_n,indn] = min(rsnesn);
[v_n2,indn2] = min(rsnesn2);

%% Plots for M3
figure(3)
clf
plot(xvec, f_alphas_0(2:end,iii), 'r')
hold on
plot(xvec, f_alphas_n(2:end,iiin), 'b-o', xvec, f_alphas_n2(2:end,iiin2), 'g-p')
plot(xvec, f(2:end), 'k', 'linewidth', 2);
xlabel(['RSNE low: ', num2str(v0), ' RSNE medium: ', num2str(vn), ' RSNE high: ', num2str(vn2) ])
legend('Low-noise','Medium-noise', 'High-noise', 'target f')
title('Reconstruction with stacked form')

figure(4)
clf
subplot(2,2,1);
plot(xvec, f(2:end), 'k', xvec, recs(2:end, ind), 'b', 'linewidth', 2)
title('Tikhonov regularization with low-noise')
subplot(2,2,2);
plot(xvec, f(2:end), 'k', xvec, recs_n(2:end, indn), 'r', 'linewidth', 2)
title('Tikhonov regularization with medium-noise')
subplot(2,2,3);
plot(xvec, f(2:end), 'k', xvec, recs_n2(2:end, indn2), 'g', 'linewidth', 2)
title('Tikhonov regularization with high-noise')
subplot(2,2,4);
plot(xvec, f(2:end), 'k', xvec, f_alphas_n2(2:end,iiin2), 'linewidth', 2)
title('stacked form Tikhonov regularization with high-noise')