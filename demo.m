clear
close all

addpath 'Toolbox';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tunable param
% PARAMETERS TO SET
R = 80;
th = 0.05;
thr = 0.5;

ratio = 20;

do_display = 0;
% For near-sources, anisotropy parameter and rough position must be pre-calibrated
% S = [-50 -85 75; 55 -90 70; -75 -45 80; 85 -45 80; -75 45 85; 80 50 75; -45 85 70; 50 80 70]; % for near sources only : initial position of the source in mm, wrt camera (measured manually, or by triangulation using reflective spheres) - set to [0;0;0] to put the source at camera center - In this demo the LED is ahead of the camera (S(3)>0), on its left (S(1)<0) and above (S(2)<0)
theta_12 = pi/180*65; % for near anisotropic sources only : theta_12 is the angle such that the intensity of the emitted light is half that of the intensity in the principal direction. Usually this angle is provided by the LED manufacturer.

e = 2e-4;
camcalib_folder = 'cc6';
images_folder = 'data10';
squareSize = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CAMERA CALIBRATION
[imageFileNames,cameraParams] = camCalib(camcalib_folder,squareSize);
% View reprojection errors
h1=figure; showReprojectionErrors(cameraParams, 'BarGraph');

% Visualize pattern locations
h2=figure; showExtrinsics(cameraParams, 'CameraCentric');

[imageFileNames,pc,g,u,v,m] = sphere_detect(images_folder,cameraParams,R,th);

k = cameraParams.IntrinsicMatrix;
k = k';
F=@(g,x)g(1)*x(:,1).^2+g(2)*x(:,1).*x(:,2)+g(3)*x(:,2).^2+g(4)*x(:,1)+g(5)*x(:,2)+g(6);

Dir = zeros(length(imageFileNames),3);
S = zeros(length(imageFileNames),3);
Ss = zeros(length(imageFileNames),3);
Phi = ones(length(imageFileNames),1);
mu = zeros(length(imageFileNames),1);
%
% Im = imread(imageFileNames{1});
% w = fspecial('gaussian',[7,7],1);
% h = fspecial('average',7);
% Im = imfilter(Im,w,'replicate');
% Im = imfilter(Im,h,'corr','replicate');
% [nrows,ncols] = size(Im);
% Z = zeros(nrows,ncols);
% for i=1:nrows
%     for j=1:ncols
%         u = j;
%         v = i;
%         P = inv(k)*[u;v;1];
%         A = sum(P.^2);
%         B = -sum(2*P.*pc');
%         C = sum(pc'.^2)-R^2;
%         delta = B^2-4*A*C;
%         if delta>0
%             z = (-B-sqrt(delta))/(2*A);
%         else
%             z = NaN;
%         end
%         Z(i,j) = z;
%     end
% end

% N = zeros(nrows*ncols,3)*NaN;
% P = zeros(nrows*ncols,3)*NaN;
% for i = 1:ratio:nrows
%     for j = 1:ratio:ncols
%         uv = [j;i;1];
%         P((j-1)*nrows+i,:) = (inv(k)*Z(i,j)*uv)';
%         N((j-1)*nrows+i,:) = P((j-1)*nrows+i,:)-pc;
%         N((j-1)*nrows+i,:) = N((j-1)*nrows+i,:)/norm(N((j-1)*nrows+i,:));
%     end
% end
% t = linspace(0,pi,25);
% p = linspace(0,2*pi,25);
% [theta,phi] = meshgrid(t,p);
% x = R*sin(theta).*sin(phi)+pc(1);
% y = R*sin(theta).*cos(phi)+pc(3);
% z = R*cos(theta)+pc(2);
% figure(99)
% surf(x,y,z)
% scatter3(pc(1),pc(3),pc(2))
% hold on;
% scatter3(P(:,1),P(:,3),P(:,2));
% axis equal

% for im=1:length(imageFileNames)
maxti = 5;
% 
% for im=1:length(imageFileNames)
%     ti = 0;
%     e = 1;
%     while(ti<maxti&&e>1e-4)
%         ti = ti+1;
%         disp([num2str(ti) 'th iteration'])
%         Im = imread(imageFileNames{im});
%         maxIm = max(max(Im));
%         w = fspecial('gaussian',[7,7],1);
%         h = fspecial('average',7);
%         Im = imfilter(Im,w,'replicate');
%         Im = imfilter(Im,h,'corr','replicate');
%         figure(23)
%         imshow(Im);
%         hold on;
%         Im = double(Im);
%         [nrows,ncols] = size(Im);
%         ff = mean([k(1,1),k(2,2)]);
%         [xx,yy] = meshgrid(1:ncols,1:nrows);
%         xx = xx-k(1,3);
%         yy = yy-k(2,3);
%         cos4a = (ff./sqrt(xx.^2+yy.^2+ff^2)).^4;
%         clear xx yy;
%         Im = bsxfun(@rdivide,Im,cos4a);
%         
%         p = zeros(length(Im),3)*NaN;
%         n = zeros(length(Im),3)*NaN;
%         mask = zeros(size(Im))*NaN;
%         for i=1:ratio:nrows
%             for j=1:ratio:ncols
%                 uv = [j-k(1,3),i-k(2,3)];
%                 if (Im(i,j)>maxIm*thr)
%                     mask(i,j) = 1;
%                     P = inv(k)*[j;i;1];
%                     A = sum(P.^2);
%                     B = -sum(2*P.*pc');
%                     C = sum(pc'.^2)-R^2;
%                     delta = B^2-4*A*C;
%                     if delta>=0
%                         Z = (-B-sqrt(delta))/(2*A);
%                     else
%                         Z = NaN;
%                     end
%                     p((j-1)*nrows+i,:) = P*Z;
%                     n((j-1)*nrows+i,:) = (P*Z-pc')/norm(P*Z-pc',2);
%                 end
%             end
%         end
%         ind = find(mask==1);
%         [lv,lu] = ind2sub(size(Im),ind);
%         plot(lu,lv,'dy','MarkerSize',1);
%         hold off;
%         I = Im(ind);
%         p = p(ind,:);
%         n = n(ind,:);
%         clear i ind;
%         [e,x0]=cali(S(im,:),theta_12,I,p,n);
%         ps = x0(1:3);
%         th = x0(4);
%         ph = x0(5);
%         dir = [cos(th)*cos(ph);sin(th)*cos(ph);sin(ph)];
%         Ss(im,:) = ps;
%         Dir(im,:) = dir;
%         Phi(im) = x0(6);
%         mu(im) = x0(7);
%         
%         S = Ss;
%     end
%     figure(h2)
%     hold on
%     ps = Ss(im,:);
%     dir = Dir(im,:);
%     plot3(ps(1),ps(3),ps(2),'d','MarkerSize',15,'MarkerEdgeColor','k','LineWidth',1);
%     quiver3(ps(1),ps(3),ps(2),dir(1),dir(3),dir(2),60,'Linewidth',2,'Color','b','MaxHeadSize',5);
%     view(0,90)
%     drawnow
% end
% 
% figure(h2)
% hold on
% t = linspace(0,pi,25);
% p = linspace(0,2*pi,25);
% [theta,phi] = meshgrid(t,p);
% x = R*sin(theta).*sin(phi)+pc(1);
% y = R*sin(theta).*cos(phi)+pc(3);
% z = R*cos(theta)+pc(2);
% surf(x,y,z)
% axis equal
% %
% % load('xyz');
% % figure(42)
% % surf(x,y,z)
% % hold on
% % surf(XYZ(:,:,1),XYZ(:,:,3),XYZ(:,:,2))
% % axis equal
% %
% clearvars -EXCEPT S Dir Phi mu k