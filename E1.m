% exercise 1
% Emil Airta
% a)
p = [1/16, 3/16, 1/2, 3/16, 1/16]';
f = zeros(size(1:16))';
f(5) = 1;
f(6) = 1;
res = zeros(size(1:16))';
for iii = 1:16
    sum = 0;
    for jjj = 1:5
        if (iii-(jjj-3) > 0)&&(iii - (jjj-3) <= 16)
            sum = sum + p(jjj)*f(iii-(jjj-3));
        end    
    end
    res(iii) = sum;
end

% b)
m = 2;
A = convmtx(p',16);
A = A(:,(m+1):(end-m));
convo  = A*f;

% c)
c = conv2(f,p, 'same');

% Problem 2
res2 = zeros(size(1:16))';
for iii = 1:16
    sum = 0;
    for jjj = 1:5
        if (iii-(jjj-3) > 0)&&(iii - (jjj-3) <= 16)
            sum = sum + p(jjj)*f(iii-(jjj-3));
        else
            tmp = iii-(jjj-3) - 16;
            if (iii-(jjj-3) <= 0)
                tmp = 16 + iii-(jjj-3);
            end
            sum = sum + p(jjj)*f(tmp);
        end    
    end
    res2(iii) = sum;
end
% b)
temp = convmtx(p',16);
A = temp(:,(m+1):(end-m));
for iii = 1:m
    A(1:m, end-m+1:end) = temp(1:m,1:m);
    A((end-m+1):end,1:m) = temp((end-m+1):end,(end-m+1):end);
end
% c) 
shiftp = zeros(size(1:16))';
shiftp(7:11) = p;
shiftp = fftshift(shiftp);
fourier = ifft(fft(shiftp).*fft(f));
