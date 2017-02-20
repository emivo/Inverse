%% Exercise 4 M1
% Emil Airta

% Parameters
n = 128;
a = .025;
% Evalution points and point spread function
delta_x = 1/n;
xvec = delta_x*(0:(n-1)); 
v = floor(a*n);
vvec = linspace(-a+delta_x,a-delta_x, v*2 + 1);
p = (PSF(vvec,a));
f = targetf(xvec);
% Calculate convolution matrix
conA = convmtx(p, n);
A = conA(:, (v+1):(end-v))*delta_x;
m_tilde = A*f';
% Add some noise to the data
sigma = .01;
noise = sigma*randn(size(m_tilde));
mn = m_tilde + noise;


% Init Ls
I = eye(n-1,n);
L = (1/delta_x)*(circshift(I,1,2) - I);
I2 = eye(n+1,n);
L0 = (1/delta_x)*(I2 - 1*circshift(I2, [1 0]));
% Init alphas
na = 20;
alphas = zeros(2*na + 1, 1);
len = length(alphas);
for k = 1:2*na
    alphas(k) = 10^(na - k);
end

% Init result variables and m prime to L
f_alphas = zeros(n, len);
rsnes = zeros(len, 1); 
m_prime = [mn; zeros(n - 1, 1)];
% Init result variables and m prime to L0
f_alphas0 = zeros(n, len);
rsnes0 = zeros(len, 1); 
m_prime0 = [mn; zeros(n + 1, 1)];
%% Calculate k�n�keppi for each alpha
for k = 1:len
    alpha = alphas(k);

    A_prime = [A; sqrt(alpha)*L];
    A_prime0 = [A; sqrt(alpha)*L0];
    
    f_alphas(:,k) = A_prime \ m_prime;
    f_alphas0(:,k) = A_prime0 \ m_prime0;
    % Relative square norm error
    rsnes(k)  = norm(f_alphas(:,k) - f')/norm(f) * 100;
    rsnes0(k)  = norm(f_alphas0(:,k) - f')/norm(f) * 100;
end

%% Plots
figure(1) 
clf
semilogx(alphas,rsnes)
hold on
semilogx(alphas,rsnes0, 'r')
legend('L','L0')
title('Relative square norm errors')

figure(3)
clf
[v,iii] = min(rsnes);
[v0,iii0] = min(rsnes0);
plot(xvec, f_alphas(:,iii), 'r','linewidth', 2)
hold on
plot(xvec, f_alphas0(:,iii0), 'b-o','markersize', 5)
plot(xvec, f, 'k', 'linewidth', 2);
xlabel(['RSNE L: ', num2str(v), ' RSNE L0: ', num2str(v0)])
legend(['f_alpha with L and alpha ', num2str(alphas(iii))],['f_alpha with L0 and alpha ', num2str(alphas(iii0))], 'target f')
title('Reconstruction with minimal RSNE')

figure(4)
clf

plot(xvec, f_alphas(:,iii + 1), 'r','linewidth', 2)
hold on
plot(xvec, f_alphas0(:,iii0 - 1), 'b-o','markersize', 5)
plot(xvec, f, 'k', 'linewidth', 2);
xlabel(['RSNE L: ', num2str(rsnes(iii + 1)), ' RSNE L0: ', num2str(rsnes0(iii0 + 1))])
legend(['f_alpha with L and alpha ', num2str(alphas(iii + 1))],['f_alpha with L0 and alpha ', num2str(alphas(iii0 - 1))], 'target f')
title('Reconstruction by other choise of alpha')