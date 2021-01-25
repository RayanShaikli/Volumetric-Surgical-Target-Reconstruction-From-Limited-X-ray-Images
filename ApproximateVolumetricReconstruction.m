function ApproximateVolumetricReconstruction()

SDD = 200;
SAD = 100;
ADD = 100;
Ock = [0; 0; 0];
S0 = Ock + [0 ; SAD ; 0];
C0 = Ock + [0 ; -ADD ; 0 ];

RzNeg = [cosd(45) sind(45) 0 ; -sind(45) cosd(45) 0 ; 0 0 1];
RzPos = [cosd(-45) sind(-45) 0 ; -sind(-45) cosd(-45) 0 ; 0 0 1];
Sa = RzPos*S0;
Sb = RzNeg*S0;
Ca = RzPos*C0;
Cb = RzNeg*C0;
Wa = (Ca-Sa)/norm(Ca-Sa);
Wb = (Cb-Sb)/norm(Cb-Sb);

image(Ca, Cb);

pts = [Ock.';Sa.';Sb.';Ca.';Cb.'];

hold on
scatter3(pts(:,1), pts(:,2), pts(:,3));
plotCircle3D(Ca.',Wa.',20);
plotCircle3D(Cb.',Wb.',20);
hold off

end

function plotCircle3D(center,normal,radius)
theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
plot3(points(1,:),points(2,:),points(3,:),'r-');
end

function image(Ov, Oy)
r = 15;
n = 100;

theta = (0:n-1)*(2*pi/n);
x1 = Ov(1) + r*cos(theta);
y1 = Ov(2) + r*sin(theta);
polyin1 = polyshape(x1,y1);

x2 = Oy(1) + r*cos(theta);
y2 = Oy(2) + r*sin(theta);
polyin2 = polyshape(x2,y2);

plot(polyin1)
plot(polyin2)
axis equal

[x11,y11] = centroid(polyin1);
[x22,y22] = centroid(polyin2);
plot(polyin1)
hold on
plot(polyin2)
hold on
plot(x11,y11,'r*')
hold on
plot(x22,y22,'r*')
hold off
end