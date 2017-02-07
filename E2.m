% Exercise 2
% Emil Airta

A = [1,1;0,1;(-1/3),1];
y = [1;-2;-2];
[U, D, V]  = svd(A);
Dplus = pinv(D);
x = V*Dplus*U'*y;
xgrid = linspace(-10,10);
l1 = 1 - xgrid;
l2 = -2*ones(size(xgrid));
l3 = -2 + 1/3*xgrid;
figure(1)
clf
plot(xgrid, l1, xgrid, l2, '--', xgrid, l3,'-.', 'LineWidth', 2)
hold on 
plot(x(1),x(2),'rx', 'MarkerSize', 15)

% M2
% a)
n = 2^14;
deltax = 1/n;
xvec = deltax*(0:(n-1)); 

a = .025;
Ny = 2048;
y = linspace(-a, a, Ny);
dy = y(2) - y(1);
v = round((a*n)-1);

vvec = linspace(-a+deltax,a-deltax, v*2 + 1);
p = (PSF(vvec,a));
conA = convmtx(p, n);
conA = conA(:, (v+1):(end-v))*deltax;
m = conA*targetf(xvec)';

mtilde = zeros(size(xvec));
for iii = 1:length(xvec) 
    x = xvec(iii);
    intg = PSF(y,a).*targetf(x-y);
    mtilde(iii) = dy*sum(intg);
end


figure(2)
clf
plot(xvec, mtilde, 'b.')
hold on
plot(xvec, m, 'r.')

rsne = norm(mtilde' - m)/norm(mtilde) * 100
% n should be between 2^10 and 2^11 to get 1% error
% n = 2^14 error was 0.09% so that would be enough