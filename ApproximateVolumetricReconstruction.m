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
% plotSystem(env,vec,centroids,M)

% Testing Base Case Generator
%  each row is sphere center
X = [0;0];
Y = [0;5];
Z = [0;5];
R = [10;5]; % radius
GenerateBaseCase(X,Y,Z,R,Sa,Sb,Ca,Cb,Wa,Wb)
end

function C = GenerateBaseCase(X,Y,Z,R,Sa,Sb,Ca,Cb,Wa,Wb)
RzNeg = [cosd(45) sind(45) 0 ; -sind(45) cosd(45) 0 ; 0 0 1];
RzPos = [cosd(-45) sind(-45) 0 ; -sind(-45) cosd(-45) 0 ; 0 0 1];

[x,y,z] = sphere;
% projection on detector A
projpts = FProjection(x*R(1)+X(1),y*R(1)+Y(1),z*R(1)+Z(1),Sa,[Ca.';Wa.']);
pgonA1 = CreatePolygon(projpts,RzNeg);
projpts = FProjection(x*R(2)+X(2),y*R(2)+Y(2),z*R(2)+Z(2),Sa,[Ca.';Wa.']);
pgonA2 = CreatePolygon(projpts,RzNeg);
pgonA = union(pgonA1,pgonA2);

% projection on detector B
projpts = FProjection(x*R(1)+X(1),y*R(1)+Y(1),z*R(1)+Z(1),Sb,[Cb.';Wb.']);
pgonB1 = CreatePolygon(projpts,RzPos);
projpts = FProjection(x*R(2)+X(2),y*R(2)+Y(2),z*R(2)+Z(2),Sb,[Cb.';Wb.']);
pgonB2 = CreatePolygon(projpts,RzPos);
pgonB = union(pgonB1,pgonB2);

figure('Name','Detector A: +45 deg about z')
plot(pgonA)

figure('Name','Detector B: -45 deg about z')
plot(pgonB)

figure('Name','Intersecting Spheres in 3D')
surf(x*R(1)+X(1),y*R(1)+Y(1),z*R(1)+Z(1));
hold on
surf(x*R(2)+X(2),y*R(2)+Y(2),z*R(2)+Z(2));
hold off
end

function pgon = CreatePolygon(P,R)
for i=1:size(P,1)
    P(i,:) = (R*P(i,:).').'  + [0 100 0];
end
P = [P(:,1),P(:,3)];
size(P)
k = convhull(P);
pgon = polyshape(P(k,1),P(k,2));
end

function projpts = FProjection(X,Y,Z,source,plane)
% pts in each row
% plane is the detector given by n and A
P = source;
A = plane(1,:).';
n = plane(2,:).';
ctr = 1;
for i=1:size(X,1)
    for j=1:size(X,1)
        pts = [X(i,j),Y(i,j),Z(i,j)];
        v = (P-(pts.'))/norm(P-(pts.'));
        t = dot((A-P),n)/dot(v,n);
        projpts(ctr,:) = P + v*t;
        ctr = ctr+1;
    end
end
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
axis equal
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
