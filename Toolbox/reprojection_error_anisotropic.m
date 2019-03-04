function [F, J, Psi] = reprojection_error_anisotropic(I,N,X,mu,var)

S = var(1:3); % Position
th = var(4);
ph = var(5);
Dir = [cos(th)*cos(ph) sin(th)*cos(ph) sin(ph)];

Xs_minus_X = bsxfun(@minus,S',X); % Vector from surface to source
norm_X_minus_Xs = sqrt((sum(Xs_minus_X.^2,2))); % Distance from source to surface

anis = (sum(bsxfun(@times,Dir,-Xs_minus_X),2)).^mu;
shading = sum(N.*Xs_minus_X,2);
shading_anis = shading.*anis;
norme = 1./(norm_X_minus_Xs.^(3+mu));
Ii = shading_anis.*norme;
Ii(Ii>4095)=4095;
mI = max(I);
mIi = max(Ii);

I = I/mI;
Ii = Ii/mIi;

sI = sum(I);
sIi = sum(Ii);

Ii = Ii/sIi*sI;

F = I-Ii;
Psi = mI/mIi*sIi/sI;

Ii = Ii*mI;

[sF,sId] = sort(F);
n = size(sF);
sF(1:ceil(0.1*n)) = 0;
sF(ceil(0.9*n):n) = 0;
F = sF(sId);

J = zeros(length(shading),5);

d_shading_x = N(:,1);
d_shading_y = N(:,2);
d_shading_z = N(:,3);
d_anis_x = -mu*(sum(bsxfun(@times,Dir,-Xs_minus_X),2)).^(mu-1).*Dir(1);
d_anis_y = -mu*(sum(bsxfun(@times,Dir,-Xs_minus_X),2)).^(mu-1).*Dir(2);
d_anis_z = -mu*(sum(bsxfun(@times,Dir,-Xs_minus_X),2)).^(mu-1).*Dir(3);
d_shading_anis_x = d_shading_x.*anis+shading.*d_anis_x;
d_shading_anis_y = d_shading_y.*anis+shading.*d_anis_y;
d_shading_anis_z = d_shading_z.*anis+shading.*d_anis_z;
d_norm_x = -(3+mu)*(S(1)-X(:,1))./(norm_X_minus_Xs.^(5+mu));
d_norm_y = -(3+mu)*(S(2)-X(:,2))./(norm_X_minus_Xs.^(5+mu));
d_norm_z = -(3+mu)*(S(3)-X(:,3))./(norm_X_minus_Xs.^(5+mu));

J(1:length(shading),1) = -(d_shading_anis_x.*norme+shading_anis.*d_norm_x);
J(1:length(shading),2) = -(d_shading_anis_y.*norme+shading_anis.*d_norm_y);
J(1:length(shading),3) = -2*(d_shading_anis_z.*norme+shading_anis.*d_norm_z).*F;


Dir_th =  [-sin(th)*cos(ph) cos(th)*cos(ph) 0];
Dir_ph =  [-cos(th)*sin(ph) -sin(th)*sin(ph) cos(ph)];
d_anis_th = shading.*norme.*((sum(bsxfun(@times,Dir,-Xs_minus_X),2)).^(mu-1)).*(sum(bsxfun(@times,Dir_th,-Xs_minus_X),2));
d_anis_ph = shading.*norme.*((sum(bsxfun(@times,Dir,-Xs_minus_X),2)).^(mu-1)).*(sum(bsxfun(@times,Dir_ph,-Xs_minus_X),2));

J(1:length(shading),4) = -2*mu*d_anis_th;
J(1:length(shading),5) = -2*mu*d_anis_ph;
end

