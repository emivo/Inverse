% Exercise 4 M2 
% Emil Airta

% Load the data
load data/deco03 m0 mn
% Differentiate the ESF
psf_emp = diff(m0);
psf_emp = [psf_emp(2:end) 0]; % shift little bit to left
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
[U,D,V] = svd(A);
svals   = diag(D);

% Compute reconstruction
r_alpha = 30;
Dp_alpha = zeros(size(A));
for iii = 1:r_alpha
    Dp_alpha(iii,iii) = 1/svals(iii);
end
fn = V*(Dp_alpha.')*(U.')*mn(:);
fn = fn/sum(fn);
fn = [fn(2:end); fn(1)]; %little shift

fn2 = V*(Dp_alpha.')*(U.')*mn(:);
fn2 = fn2/sum(fn2);
fn2 = [fn2(2:end); fn2(1)];%little shift

% Normalization PSF
C = sum(psf_emp);
p = psf_emp/C;


% Theoretical PSF
xvec = linspace(-127, 127, 255);
psf = PSF(xvec, floor(dy/2));

% Symmetrice PSF
psf_sym1 = (p + fliplr(p))/2;
psf_sym2 = (p +flipud(p))/2;

% Plot curves
figure(1)
clf
plot(xvec, p, 'b-o', xvec, psf, 'r')
hold on
plot(xvec, fn(2:end), 'k-*', xvec, fn2(2:end), 'g-p')
xlim([-dy*2 dy*2])
legend('Low-noise PSF', 'Theoretical PSF', 'Medium noise tSVD PSF', 'High noise tSVD PSF')
figure(2)
clf
plot(xvec, psf_sym1, 'b', xvec, psf_sym2, 'k', xvec, psf, 'r', xvec, fn(2:end), xvec, fn2(2:end))
legend('Low-noise symmetric PSF', 'Low-noise symmetric PSF', 'Theoretical PSF', 'Medium noise tSVD PSF', 'High noise tSVD PSF')
xlim([-dy*2 dy*2])
