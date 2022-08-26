function [T,N,B,k,t] = frenet(x,y,z)

if (nargin==2)
    z = zeros(size(x));
end

% CONVERT TO COLUMN VECTOR
x=x(:);
y=y(:);
z=z(:);

% SPEED OF CURVE
dx=gradient(x);
dy=gradient(y);
dz=gradient(z);
dr=[dx,dy,dz];
ddx=gradient(dx);
ddy=gradient(dy);
ddz=gradient(dz);
ddr=[ddx,ddy,ddz];

% TANGENT
T=dr./mag(dr,3);

% DERIVIATIVE OF TANGENT
dTx=gradient(T(:,1));
dTy=gradient(T(:,2));
dTz=gradient(T(:,3));
dT=[dTx,dTy,dTz];

% NORMAL
N=dT./mag(dT,3);

% BINORMAL
B=cross(T,N);

% CURVATURE
k=mag(cross(dr,ddr),1)./((mag(dr,1)).^3);

% TORSION
t=dot(-B,N,2);

function N = mag(T,n)
% MAGNATUDE OF A VECTOR (Nx3)
N=sqrt(sum(abs(T).^2,2));
d=find(N==0); 
N(d)=eps*ones(size(d));
N=N(:,ones(n,1));
