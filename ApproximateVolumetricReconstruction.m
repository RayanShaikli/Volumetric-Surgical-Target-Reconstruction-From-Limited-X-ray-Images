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
