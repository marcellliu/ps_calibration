clear
close all

load('calibrationSession.mat')
R = 80;
k = calibrationSession.CameraParameters.IntrinsicMatrix;
k(2,1) = 0;
k = k';

% ????
I = imread('Image__2019-02-28__13-37-58.png');
[I,newOrigin] = undistortImage(I,calibrationSession.CameraParameters);
[nrows,ncols] = size(I);

Ig = I;
th = 9;
Ig(Ig>th) = 255;
Ig(Ig<th) = 0;

Omega_padded = padarray(Ig,[1 1],0);

% Pixels who have bottom neighbor in mask
Ib = Omega_padded(3:end,2:end-1);
% Pixels who have top neighbor in mask
It = Omega_padded(1:end-2,2:end-1);
% Pixels who have right neighbor in mask
Ir = Omega_padded(2:end-1,3:end);
% Pixels who have left neighbor in mask
Il = Omega_padded(2:end-1,1:end-2);

%????
a = abs(Ib-Ig);
b = abs(Ir-Ig);
c = abs(It-Ig);
d = abs(Il-Ig);
e = a+b+c+d;
ia = find(a>0);
ib = find(b>0);
ic = find(c>0);
id = find(d>0);
ie = find(e>250);
% ind = [ia;ib;ic;id];
ind = ie;
ind = unique(ind);
[v,u] = ind2sub(size(I),ind);

cc = [1;1;1];
cc(1) = sum(v)/length(v);
cc(2) = sum(u)/length(u);

subplot(2,3,1);
imshow(I);
subplot(2,3,2);
imshow(I);
hold on;
plot(u,v,'dy','MarkerSize',1);
hold off;
subplot(2,3,3);
imshow(Ig);

% ????
p = [1 1 1 1 1 1];
p0 = p;
F=@(p,x)p(1)*x(:,1).^2+p(2)*x(:,1).*x(:,2)+p(3)*x(:,2).^2+p(4)*x(:,1)+p(5)*x(:,2)+p(6);
warning off

for i= 1:100
    uv = [u,v];
    p = nlinfit(uv,zeros(length(u),1),F,p0);
    p0 = p;
end

UpMinx=min(u);
UpMaxx=max(u);
UpMiny=min(v);
UpMaxy=max(v);

m = zeros(size(I));
for  i=1:nrows
    for j=1:ncols
        x = [j,i];
        if abs(F(p,x))< 5e-30
            m(i,j) = 1;
        end
    end
end

ind = find(m==1);
[v,u] = ind2sub(size(I),ind);

subplot(2,3,4);
imshow(I);
hold on;
plot(u,v,'dy','MarkerSize',1);

pt = zeros(length(u),3);
ft = zeros(length(u),3);
for j=1:length(u)
    uv = [u(j);v(j);1];
    r = sqrt(sum((inv(k)*uv).^2));
    pt(j,:) = inv(k)*uv/r;
    ft(j,:) = inv(k)*uv;
end

pc = sum(pt,1)/length(u);
uvc = k*pc'/pc(3);

theta = zeros(length(u),1);
for i=1:length(u)
    a = pt(1,:);
    b = pc;
    theta(i) = acos(dot(a,b)/norm(a,2)/norm(b,2));
end
theta = sum(theta)/length(u);

r = 0;
for i=1:length(u)
    r = r+sqrt(sum((pt(i,:)-pc).^2));
end
r = r/length(u);
lambda = R*cos(theta)/r;

pc = pc*lambda+R*sin(theta)*pc/norm(pc,2);
pt = pt*lambda;

t = zeros(length(u),1);
for i=1:length(u)
    t(i) = sqrt(sum((pt(i,:)-pc).^2));
end

p = zeros(length(I),3)*NaN;
n = zeros(length(I),3)*NaN;
thr = 20;
mask = zeros(size(I))*NaN;
for i=1:10:nrows
    for j=1:10:ncols
        if (I(i,j)>thr)
            mask(i,j) = 1;
            P = inv(k)*[j;i;1];
            A = sum(P.^2);
            B = -sum(2*P.*pc');
            C = sum(pc'.^2)-R^2;
            delta = B^2-4*A*C;
            if delta>=0
                z = (-B-sqrt(delta))/(2*A);
            else
                z = NaN;
            end
            p((j-1)*nrows+i,:) = P*z;
            n((j-1)*nrows+i,:) = (P*z-pc')/norm(P*z-pc',2);
        end
    end
end
% p=rmmissing(p);
% n=rmmissing(n);

ind = find(mask==1);
[mv,mu] = ind2sub(size(I),ind);

subplot(2,3,5);
imshow(I);
hold on
plot(u,v,'dy','MarkerSize',1);
plot(mu,mv,'dg','MarkerSize',1);

subplot(2,3,6);
scatter3(p(:,1),p(:,2),p(:,3));
hold on
scatter3(pt(:,1),pt(:,2),pt(:,3));
axis equal