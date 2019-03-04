function [imageFileNames,pc,g,u,v,cos4a] = sphere_detect(images_folder,cameraParams,R,th,E,thr)
images = dir(fullfile(images_folder, '*.png'));
nb_images = length(images);
for im = 1:nb_images
    if images(im).name == "all.png"
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

Ig = all;
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
[v,u] = ind2sub(size(all),ind);

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

% ????
g = [1 1 1 1 1 1];
g0 = g;
F=@(g,x)g(1)*x(:,1).^2+g(2)*x(:,1).*x(:,2)+g(3)*x(:,2).^2+g(4)*x(:,1)+g(5)*x(:,2)+g(6);
warning off

for i= 1:10
    uv = [u,v];
    g = nlinfit(uv,zeros(length(u),1),F,g0);
    g0 = g;
end

m = zeros(size(all));
for  i=1:nrows
    for j=1:ncols
        x = [j,i];
        m(i,j) = F(g,x);
    end
end

m = m/max(max(abs(m)));

for  i=1:nrows
    for j=1:ncols
        if abs(m(i,j)) < E
            m(i,j) = 1;
        else
            m(i,j) = 0;
        end
    end
end

ind = find(m==1);
[v,u] = ind2sub(size(all),ind);

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

figure(224)
imshow(all);
hold on;
plot(u,v,'dy','MarkerSize',1);
plot(uvc(1),uvc(2),'*b','MarkerSize',3);
plot(k(1,3),k(2,3),'*r','MarkerSize',3);
end
