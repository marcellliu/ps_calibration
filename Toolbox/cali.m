function [normalized_shading,x0]=cali(S,theta_12,I,X,N)
X_minus_Xs = bsxfun(@minus,X,S); % Vector from source to surface
norm_X_minus_Xs = sqrt((sum(X_minus_Xs.^2,2))); % Distance from source to surface
shading = sum(-bsxfun(@times,N,X_minus_Xs./norm_X_minus_Xs),2); % Shading: normal to surface times normalized lighting vector
normalized_shading = shading./(norm_X_minus_Xs.^2); % Shading divided by squared source-surface distance

mu = -log(2)./log(cos(theta_12));
[Dir,~,Phi] = svds(pinv(X_minus_Xs./norm_X_minus_Xs)*((bsxfun(@rdivide,I,normalized_shading)).^(1/mu)),1);
if(Phi <0)
    Dir = -Dir;
end
[th,ph] = cart2sph(Dir(1),Dir(2),Dir(3));
x0 = [S';th;ph;0;mu]; % Initial position, intensities and orientations

n = 0;
it_el = 50000;
e_el = 1e-6;
e = e_el*10;

et = 9999;
fl = 1;
o = 0;
sumet = 0;
eb = [];
E = [];

T = 1;
Tl = 1e-03;
while(T>Tl)
    r = T;
    while(n<it_el&&e>e_el&&o<2)
        [F,J,phi,Ii] = reprojection_error_anisotropic(I,N,X,mu,x0);
        F(isnan(F)==1) = [];
        e = sum(F.^2)/length(F);
        x0(6) = phi;
        switch fl
            case 1
                update_x = sum(J(:,1));
                x0(1) = x0(1) - sign(update_x)*r;
            case 2
                update_y = sum(J(:,2));
                x0(2) = x0(2) - sign(update_y)*r;
            case 3
                update_z = sum(J(:,3));
                x0(3) = x0(3) - sign(update_z)*r;
            case 4
                update_th = sum(J(:,4));
                x0(4) = x0(4) - sign(sum(update_th))*pi/180*r;
            case 5
                update_ph = sum(J(:,5));
                x0(5) = x0(5) - sign(sum(update_ph))*pi/180*r;
        end
        if length(eb)<20
            eb = [eb e];
        else
            sume = sum(eb);
            if ((sume-sumet)/sume)<1e-03
                o = o+1;
            end
            eb = [];
            sumet = sume;
        end
        
%             disp([' e: ' num2str(e) ' x: ' num2str(x0(1)) ' y: ' num2str(x0(2)) ' z: ' num2str(x0(3)) ' th: ' num2str(x0(4)) ' ph: ' num2str(x0(5)) ' phi: ' num2str(x0(6)) ' ]']);
        figure(47);
            plot(I)
            hold on 
            plot(Ii)
            hold off
%         E = [E e];

        n = n+1;
        fl = rem(fl,5)+1;
    end
    T = T*0.9;
    n = 0;
    o = 0;
end
disp([' e: ' num2str(e) ' x: ' num2str(x0(1)) ' y: ' num2str(x0(2)) ' z: ' num2str(x0(3)) ' th: ' num2str(x0(4)) ' ph: ' num2str(x0(5)) ' phi: ' num2str(x0(6)) ' ]']);
end