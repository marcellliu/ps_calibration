function [normalized_shading,x0]=cali(S,theta_12,I,X,N)
X_minus_Xs = bsxfun(@minus,X,S); % Vector from source to surface
norm_X_minus_Xs = sqrt((sum(X_minus_Xs.^2,2))); % Distance from source to surface
shading = sum(-N.*X_minus_Xs./norm_X_minus_Xs,2); % Shading: normal to surface times normalized lighting vector
normalized_shading = shading./(norm_X_minus_Xs.^2); % Shading divided by squared source-surface distance

mu = -log(2)./log(cos(theta_12));
[Dir,int,Phi] = svds(pinv(X_minus_Xs./norm_X_minus_Xs)*((bsxfun(@rdivide,I(:,3),normalized_shading)).^(1/mu)),1);
Phi = int*Phi;
if(Phi <0)
    Phi = -Phi;
    Dir = -Dir;
end
Phi = Phi.^mu;
[th,ph] = cart2sph(Dir(1),Dir(2),Dir(3));

x0 = [S';th;ph]; % Initial position, intensities and orientations

n = 0;
fl = 0;

it_el = 50000;
e_el = 1e-6;
e = e_el*10;

et = 9999;
xt = [];

r = 0.3;
flage = 3;
tt = 0;

while(n<it_el&&e>e_el&&fl<9)
    [F,J,psi] = reprojection_error_anisotropic(I(:,3),N,X,mu,x0);
    e = sum(F.^2)/length(F);
    de = et-e;
    
    if (de<0)
        flage = rem(flage+1,5);
        x0 = xt;
        [F,J,psi] = reprojection_error_anisotropic(I(:,3),N,X,mu,x0);
        tt = tt+1;
        if tt == 4
            x0(1:3) = x0(1:3)+0.1*r*randn(size(x0(1:3)));
            tt = 0;
            fl = fl+1;
        end
    else
        tt = 0;
        et = e;
        xt = x0;
    end
    
    switch flage
        case 0
            update_x = sum(J(:,1));
            x0(1) = x0(1) - sign(update_x)*r;
        case 1
            update_y = sum(J(:,2));
            x0(2) = x0(2) - sign(update_y)*r;
        case 2
            update_z = sum(J(:,3));
            x0(3) = x0(3) - sign(update_z)*r;
        case 3
            update_th = sum(J(:,4));
            x0(4) = x0(4)+2*pi;
            x0(4) = rem(x0(4) - sign(sum(update_th))*pi/180*0.8*r,2*pi);
        case 4
            update_ph = sum(J(:,5));
            x0(5) = rem(x0(5) - sign(sum(update_ph))*pi/180*0.8*r,pi);
            if (x0(5)<0)
                x0(5) = -x0(5);
            end
    end
    
    disp(['flage: ' num2str(flage) ' de: ' num2str(de) ' e: ' num2str(e) ' x: ' num2str(x0(1)) ' y: ' num2str(x0(2)) ' z: ' num2str(x0(3)) ' th: ' num2str(x0(4)) ' ph: ' num2str(x0(5)) ' ]']);
    figure(47);
    plot(F);
    hold on;
    plot(sort(F));

    hold off;
    drawnow
    pause(0.1)
    n = n+1;
end

end