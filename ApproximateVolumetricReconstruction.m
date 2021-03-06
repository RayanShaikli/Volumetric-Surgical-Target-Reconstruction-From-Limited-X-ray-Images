function ApproximateVolumetricReconstruction()

% Environment setup
% SDD: source detector distance
% SAD: source axis distance
% ADD: axis detector distance
% S: source
% C: detector center
% W: normal vector from detector
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

polyin = CreateCircles(2,C0);
centroids = FindCentroids(polyin);
M = SymbolicIntersection(centroids, [Sb.';Sa.']);

env = [Ock.';Sa.';Sb.';Ca.';Cb.'];
vec = [Wa.';Wb.'];
plotSystem(env,vec,centroids,M)

end

function SI = SymbolicIntersection(centroids, sources)
v1 = (sources(1,:)-centroids(1,:))./abs(sources(1,:)-centroids(1,:));
v2 = (sources(2,:)-centroids(2,:))./abs(sources(2,:)-centroids(2,:));
v1(3) = 0;
v2(3) = 0;
v3 = cross(v1,v2);

P = centroids(1,:).' - centroids(2,:).';
V = [-(v1.'), v2.', v3.'];
t = V\P;

L1 = v1*t(1) + centroids(1,:);
L2 = v2*t(2) + centroids(2,:);

SI = (L1+L2)/2;
end

function plotSystem(env,vec,centroids,M)
hold on
scatter3(M(1),M(2),M(3),'filled','cyan');
scatter3(env(:,1), env(:,2), env(:,3),'black');
plotCircle3D(env(4,:),vec(1,:),15);
plotCircle3D(env(5,:),vec(2,:),15);

for i=1:size(centroids,1)
    plot(centroids(i,1),centroids(i,2),'r*')
end

Cone(env(2,:),env(4,:),[0 15],20,'none',0,1);
Cone(env(3,:),env(5,:),[0 15],20,'none',0,1);
hold off
end

function plotCircle3D(center,normal,radius)
theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
plot3(points(1,:),points(2,:),points(3,:),'r-');
end

function polyshapes = CreateCircles(NCircle,C0)
r = 15;
n = 100;
theta = (0:n-1)*(2*pi/n);
for i=1:NCircle
    x = C0(1) + r*cos(theta);
    y = C0(2) + r*sin(theta);
    polyin = polyshape(x,y);
    polyshapes(i,:) = polyin;
end
end

function centroids = FindCentroids(polyin)
RzNeg = [cosd(45) sind(45) 0 ; -sind(45) cosd(45) 0 ; 0 0 1];
RzPos = [cosd(-45) sind(-45) 0 ; -sind(-45) cosd(-45) 0 ; 0 0 1];

[x,y] = centroid(polyin(1,:));
centroids(1,:) = RzNeg*[x;y;0];

[x,y] = centroid(polyin(1,:));
centroids(2,:) = RzPos*[x;y;0];
end
