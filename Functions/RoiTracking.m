function roi = RoiTracking(rt,lpts, rpts,lbpt,rbpt,lept,rebt, camMatrix, w )
%rt, predefined affine radio
%lpts,rpts: the retrieved 3d pts respectively on left and right half faces
%except the pts on mouth and cheeks.
%lbpt,rbpt, the inner corner points of left and right eyebrows.
%lept,rebt, the inner corner points of left and right eyes.
% camMatrix, the parameters of two camerias acquired through the clibaration.
% w, scale factor  acquired through the clibaration.
roi = [];
midPts = (lpts+rpts)/2;

%the calculation of mid-plane
p = mean(midPts,1);
R = bsxfun(@minus,midPts,p);
[V,~] = eig(R'*R);
n = V(:,1);%the normal of the sysmetry plane
%V = V(:,2:end);

dis = norm(lept-rebt);
pt1 = (lept+rebt)/2 + dis*n/2;
pt2 = (lept+rebt)/2 - dis*n/2;

pt3 = (lbpt+rbpt)/2 + dis*n/2;
pt4 = (lbpt+rbpt)/2 - dis*n/2;

pt1_2d = w*[pt1;1]*camMatrix;
pt2_2d = w*[pt2;1]*camMatrix;
pt3_2d = w*[pt3;1]*camMatrix;
pt4_2d = w*[pt4;1]*camMatrix;

pt5_2d = pt3_2d+rt*(pt3_2d-pt1_2d);
pt6_2d = pt4_2d+rt*(pt2_2d-pt1_2d);
roi = [roi;pt3_2d;pt4_2d;pt5_2d;pt6_2d];