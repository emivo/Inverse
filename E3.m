% Exercise 3
% Emil Airta
%% Init
n = 2^6;
deltax = 1/n;
xvec = deltax*(0:(n-1)); 

a = .025;
v = round((a*n)-1);
vvec = linspace(-a+deltax,a-deltax, v*2 + 1);
p = (PSF(vvec,a));
convA = convmtx(p, n);
convA = convA(:, (v+1):(end-v))*deltax;
m = convA*targetf(xvec)';
sigma = .01;
noise = sigma*randn(size(m));
mn = m + noise;
f = targetf(xvec);
%% Problem 1
% a
alpha = 1; 
A_prime = [convA; sqrt(alpha)*eye(n)];
m_prime = [mn; zeros(n, 1)];
f_alpha = A_prime\m_prime;
% size(A_prime)
% size(m_prime)
% b
nalphas = 12;
f_alphas = zeros(size(nalphas,n));
for iii = 1:nalphas
    alpha = 10^(iii - nalphas/2);
    A_prime = [convA; sqrt(alpha).*eye(n)];
    m_prime = [mn; zeros(n, 1)];
     size(A_prime)
     size(m_prime)
    f_alphas(iii) = A_prime\m_prime;  
end
rsne = zeros(nalphas,1);
for iii = 1:alphas
    rsne(iii) = norm(f_alphas(iii) - f)/norm(f) * 100;
end
%% Plot
figure(1) 
clf 
% plot a
plot(xvec, mn, 'r.', xvec, f, 'k', xvec, f_alpha, 'b')
legend('Discrete convolution data', 'f', 'f_alpha')
%  plot b
figure(2)
clf
plot(1:nalphas, rsne, 'r.', xvec, targetf(xvec), 'k', xvec, f_alpha, 'b')








