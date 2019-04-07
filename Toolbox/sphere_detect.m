function [imageFileNames,pc,g,u,v,cos4a] = sphere_detect(images_folder,cameraParams,R,th)
images = dir(fullfile(images_folder, '*.tiff'));
nb_images = length(images);
for im = 1:nb_images
    if strcmp(images(im).name,'all.tiff')
        all = imread(fullfile(images_folder,images(im).name));
    else
        imageFileNames{im} = fullfile(images_folder,images(im).name);
    end
end
clear images

[all,~] = undistortImage(all,cameraParams);
[nrows,ncols] = size(all);

% all = double(all);

k = cameraParams.IntrinsicMatrix;
k = k';

ff = mean([k(1,1),k(2,2)]);
[xx,yy] = meshgrid(1:ncols,1:nrows);
xx = xx-k(1,3);
yy = yy-k(2,3);
cos4a = (ff./sqrt(xx.^2+yy.^2+ff^2)).^4;
% all = bsxfun(@rdivide,all,cos4a);

% all = int8(all);

Ig = double(all);
Ig = Ig/max(max(Ig));
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
ie = find(e>th);
% ind = [ia;ib;ic;id];
ind = ie;
ind = unique(ind);
[v,u] = ind2sub(size(all),ind);
clear a b c d e ia ib ic id ie;

cc = [1;1;1];
cc(1) = sum(v)/length(v);
cc(2) = sum(u)/length(u);

figure(221)
imshow(all);
figure(222)
imshow(all);
hold on;
plot(u,v,'dy','MarkerSize',1);
hold off;
figure(223)
imshow(Ig);
imwrite(Ig,'mask.png');

% ????
g = [1 1 1 1 1 1];
g0 = g;
F=@(g,x)g(1)*x(:,1).^2+g(2)*x(:,1).*x(:,2)+g(3)*x(:,2).^2+g(4)*x(:,1)+g(5)*x(:,2)+g(6);
warning off

u = u-k(1,3);
v = v-k(2,3);

uv = [u,v];
g = nlinfit(uv,zeros(length(u),1),F,g0);
g0 = g;

u = zeros(1000,1)*NaN;
v = zeros(1000,1)*NaN;

for i=1:500
    theta = pi/500*i;
    A = (g(1)*cos(theta)^2+g(2)*cos(theta)*sin(theta)+g(3)*sin(theta)^2);
    B = g(4)*cos(theta)+g(5)*sin(theta);
    C = g(6);
    delta = B^2-4*A*C;
    r = (-B+sqrt(delta))/(2*A);
    u(i) = r*cos(theta);
    v(i) = r*sin(theta);
    r = (-B-sqrt(delta))/(2*A);
    u(end-i+1) = r*cos(theta);
    v(end-i+1) = r*sin(theta);
end
clear A B C delta r;

u = u+k(1,3);
v = v+k(2,3);

pt = zeros(length(u),3);
for j=1:length(u)
    uv = [u(j);v(j);1];
    pt(j,:) = inv(k)*uv;
    pt(j,:) = pt(j,:)/norm(pt(j,:),2);
end

% ???????????????????? ?????????????? ???????? ???????????????
P=@(g,x)g(1)*x(:,1)+g(2)*x(:,2)+g(3)*x(:,3)+g(4);
h = [1 1 1 1];
h0 = h;
h = nlinfit(pt,zeros(length(u),1),P,h0);
n = [h(1) h(2) h(3)];
O = -h(4)/(h(1)^2+h(2)^2+h(3)^2);
O = O*n;
rf = sqrt(1-norm(O)^2);
i1 = [-h(4)/h(1) 0 0];
i2 = [0 -h(4)/n(2) 0];
a = i1-i2;
b = cross(a,n);
b = b/norm(b);
a = a/norm(a);
ptf = zeros(1000,3);
for i = 1:1000
    theta = i*2*pi/1000;
    ptf(i,:) = O+a*rf*cos(theta)+b*rf*sin(theta);
end

pt = ptf;
pc = mean(pt);
uvc = k*pc'/pc(3);

% f = pt;
% f = f(:,1).^2+f(:,2).^2+f(:,3).^2;
% f = sqrt(f);
% figure(67)
% plot(f)

r = 0;
for i=1:length(u)
    r = r+sqrt(sum((pt(i,:)-pc).^2));
end
r = r/length(u);

theta = zeros(length(u),1);
for i=1:length(u)
    a = pt(1,:);
    b = pc;
    theta(i) = acos(dot(a,b)/norm(a,2)/norm(b,2));
end
theta = sum(theta)/length(u);

lambda = R*cos(theta)/r;
l = lambda+R*sin(theta)/norm(pc,2);
pc = pc*l;
pt = pt*l;


u = round(u);
v = round(v);
uvc = round(uvc);

figure(224)
imshow(all);
hold on;
plot(u,v,'dy','MarkerSize',1);
plot(uvc(1),uvc(2),'*b','MarkerSize',3);
plot(k(1,3),k(2,3),'*r','MarkerSize',3);
end

