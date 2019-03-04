clear
close all

addpath 'Toolbox';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS TO SET
R = 80;

% Tunable param
th = 253;
thr = 50;

% For near-sources, anisotropy parameter and rough position must be pre-calibrated
S = [-40 -80 110; 40 -80 110; -80 -30 85; 70 -30 80]; % for near sources only : initial position of the source in mm, wrt camera (measured manually, or by triangulation using reflective spheres) - set to [0;0;0] to put the source at camera center - In this demo the LED is ahead of the camera (S(3)>0), on its left (S(1)<0) and above (S(2)<0)
theta_12 = pi*(1/3); % for near anisotropic sources only : theta_12 is the angle such that the intensity of the emitted light is half that of the intensity in the principal direction. Usually this angle is provided by the LED manufacturer.

e = 4e-4;
camcalib_folder = 'cc/';
images_folder = 'data0/';
squareSize = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CAMERA CALIBRATION
[imageFileNames,cameraParams] = camCalib(camcalib_folder,squareSize);
% View reprojection errors
h1=figure; showReprojectionErrors(cameraParams, 'BarGraph');

% Visualize pattern locations
h2=figure; showExtrinsics(cameraParams, 'CameraCentric');

[imageFileNames,pc,g,u,v,m] = sphere_detect(images_folder,cameraParams,R,th,e,thr);

k = cameraParams.IntrinsicMatrix;
k = k';
F=@(g,x)g(1)*x(:,1).^2+g(2)*x(:,1).*x(:,2)+g(3)*x(:,2).^2+g(4)*x(:,1)+g(5)*x(:,2)+g(6);

s = zeros(length(imageFileNames),5);
for im=1:4
    Im = imread(imageFileNames{im});
    figure('name',(imageFileNames{im}));
    imshow(Im)
    hold on;
    
    Im = double(Im);
    [nrows,ncols] = size(Im);
    ff = mean([k(1,1),k(2,2)]);
    [xx,yy] = meshgrid(1:ncols,1:nrows);
    xx = xx-k(1,3);
    yy = yy-k(2,3);
    cos4a = (ff./sqrt(xx.^2+yy.^2+ff^2)).^4;
    clear xx yy;
    Im = bsxfun(@rdivide,Im,cos4a);
    
    p = zeros(length(Im),3)*NaN;
    n = zeros(length(Im),3)*NaN;
    mask = zeros(size(Im))*NaN;
    for i=1:50:nrows
        for j=1:50:ncols
            x = [j,i];
            if (Im(i,j)>thr && F(g,x)<0)
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
    
    ind = find(mask==1);
    [mv,mu] = ind2sub(size(Im),ind);
    i = Im(ind);
    I = [mv,mu,i];
    p = p(ind,:);
    n = n(ind,:);
    clear i ind;
   
    plot(u,v,'dy','MarkerSize',1);
    plot(mu,mv,'dg','MarkerSize',1);
    [normalized_shading,x0]=cali(S(im,:),theta_12,I,p,n);
    s(im,:) = x0;
%     te = bsxfun(@rdivide,I(:,3),normalized_shading)
% X = p;
% N = n;
% X_minus_Xs = bsxfun(@minus,X,S'); % Vector from source to surface
% norm_X_minus_Xs = sqrt((sum(X_minus_Xs.^2,2))); % Distance from source to surface
% shading = sum(-N.*X_minus_Xs./norm_X_minus_Xs,2); % Shading: normal to surface times normalized lighting vector
% % normalized_shading = shading./(norm_X_minus_Xs.^2); % Shading divided by squared source-surface distance
% 
% quiver3(p(:,1),p(:,2),p(:,3),n(:,1),n(:,2),n(:,3),1,'Linewidth',1,'Color','b','MaxHeadSize',5)
% hold on
% quiver3(p(:,1),p(:,2),p(:,3),X_minus_Xs(:,1),X_minus_Xs(:,2),X_minus_Xs(:,3),1,'Linewidth',1,'Color','g','MaxHeadSize',5)

end

