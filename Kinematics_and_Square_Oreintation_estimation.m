%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                                                                 %%%%%
%%%%%                ROBOTICS - Dynamics and Control                  %%%%%
%%%%%                           MID TERM                              %%%%%
%%%%%                                                                 %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MidTerm_2017
close all
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Initial Parameters.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Angles Measured. [q1 q2 q3; q1 q2 q3]
A  =[-1.143 -0.243 0.063; -1.140 0.728 -1.150]';
P_B=[-1.276 0.245 -0.504; -1.277 -0.171 0.062]';
P_C=[-1.130 -0.237 0.359; -1.131 0.245 -0.322]';
D  =[-1.010 -0.287 0.288; -1.009 -0.753 0.888]';
Pa=vpa(Solve(A),4)
b=vpa(Solve(P_B),4);
Pc=vpa(Solve(P_C),4);
Pd=vpa(Solve(D),4);
hold on
plot3(b(1,2),b(2,2),b(3,2),'*')
pause(0.1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms t
t=t;
[X_d,~]=eqn(D(:,2),t);
r=13*2^(1/2);
i=2;
Eqn=(X_d(1)-b(1,i))^2+(X_d(2)-b(2,i))^2+(X_d(3)-b(3,i))^2-2*(13)^2;
td=solve(Eqn==0,t);
d_1=vpa(subs(X_d,t,td(1)),4);
plot3([b(1,2) d_1(1)],[b(2,2) d_1(2)],[b(3,2) d_1(3)])
pause(0.1)
d_2=vpa(subs(X_d,t,td(2)),4);
plot3([b(1,2) d_2(1)],[b(2,2) d_2(2)],[b(3,2) d_2(3)])
pause(0.1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Soilving for second line
O1=(b(:,i)+d_1)/2;
O2=(b(:,i)+d_2)/2;
V1=(b(:,i)-d_1)/norm(b(:,i)-d_1);
V2=(b(:,i)-d_2)/norm(b(:,i)-d_2);
l=13/sqrt(2);
syms t1 t2
t1=t1;
t2=t2;
C1=circle(O1,V1,t1,l)';
C2=circle(O2,V2,t2,l)';
Theta=0:pi/35:2*pi;
Cir=subs(C1,t1,Theta)';
plot3(Cir(:,1),Cir(:,2),Cir(:,3))
pause(0.1)
Cir=subs(C2,t2,Theta)';
plot3(Cir(:,1),Cir(:,2),Cir(:,3))
pause(0.1)
syms t3
t3=t3;
L_C=eqn(P_C(:,2),t3);
[Opt1,Dist1]=fminsearch(@distance1,[0;1000]);
[Opt2,Dist2]=fminsearch(@distance2,[0;1000]);
if(Dist1>Dist2)
pc=vpa((subs(L_C,t3,Opt2(2))),5);
d=d_2;
du=d_1;
else
pc=vpa(subs(L_C,t3,Opt1(2)),3);
d=d_1; 
du=d_2;
end
plot3([d(1) pc(1)],[d(2) pc(2)],[d(3) pc(3)],'*')
pause(0.1)
fprintf('Point B \n')
disp(b(:,2))
fprintf('Point D \n')
disp(d)
fprintf('Point C \n')
disp(pc)
syms x y z
x=x;y=y;z=z;
Plane=plane(pc,d,b(:,2),x,y,z);
[Opt3,Delta]=fminsearch(@Finalpoint,double(b(:,2)+d+pc)/3);
pa=vpa(Opt3,5);view([1,3,0.25])
fprintf('Point A \n')
disp(pa)
plot3(pa(1),pa(2),pa(3),'*')
axis('equal')
pause(0.1)
Square=fill3([pa(1) b(1,2) pc(1) d(1) pa(1)],[pa(2) b(2,2) pc(2) d(2) pa(2)],[pa(3) b(3,2) pc(3) d(3) pa(3)],1)
Square.FaceAlpha=0.5
A=text(pa(1),pa(2),pa(3),'A')
B=text(b(1,2),b(2,2),b(3,2),'B')
C=text(pc(1),pc(2),pc(3),'C')
D=text(d(1),d(2),d(3),'D')
Du=text(du(1),du(2),du(3),'D fake')
A.FontSize=15;
B.FontSize=15;
C.FontSize=15;
D.FontSize=15;
Du.FontSize=8;
view([1,3,0.25])
axis('equal')
xlabel('X-axis');ylabel('Y-axis');zlabel('Z-axis');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% P0 and P for Point A.
function R=Solve(A)
syms T1 T2
[X_1,M1]=eqn(A(:,1),T1);
[X_2,M2]=eqn(A(:,2),T2);
X=X_1-X_2;
x1=dot(M1,X);
x2=dot(M2,X);
[C,B]=equationsToMatrix(x1,x2,[T1 T2]);
T=C\B;
T1=T(1);
T2=T(2);
x1=eval(X_1);
x2=eval(X_2);
R=[x1 x2];
end
function [X,M1]=eqn(A,t)
q1=A(1,1);q2=A(2,1);q3=A(3,1);
P_T=puma(q1,q2,q3);
P1=P_T(1:3,4);
M1=P_T(1:3,1);
X=P1+M1*t;
end
function [R]=puma(q1,q2,q3)
d1=666/25.4; 
d2=243.5/25.4;
a2=431.8/25.4;
d3=93.4/25.4;
b1=270/25.4;
k=52.1/25.4;
R=Dhcon(q1,d1,0,-pi/2)*Dhcon(q2,-d2,a2,0)*Dhcon(q3,d3+k,b1,0);
function [ R ] = Dhcon( thetaZ,dz,dx,thetaX )
%Returns the DH convection matrix
R=Rotz(thetaZ)*Transz(dz)*Transx(dx)*Rotx(thetaX);
function [ R ] = Rotx( theta )
%Return the rotation matrix by Rotation of X in a 3D frame where theta is
%in radians
R=[1 zeros(1,3);
    0 cos(theta) -sin(theta) 0;
    0 sin(theta) cos(theta) 0;
    zeros(1,3) 1;];
end
function [ R ] = Rotz( theta )
%Return the rotation matrix by Rotation of X in a 3D frame where theta is
%in radians
R=[cos(theta) -sin(theta) 0 0;
   sin(theta) cos(theta) 0 0;
   0 0 1 0; 
   zeros(1,3) 1;];
end
function [ R ] = Transx( d )
%Return the translation matrix by Translation of X in a 3D frame where theta is
R=[eye(3) [d;0;0;]
    zeros(1,3) 1];
end
function [ R ] = Transz( d )
%Return the translation matrix by Translation of Z in a 3D frame where theta is
R=[eye(3) [0;0;d;]
    zeros(1,3) 1];
end
end

end
function [R]=circle(O,V1,t,r)
A1=[V1(2);-V1(1);0];
A1=A1/norm(A1);
B1=cross(V1,A1);
R(1)=O(1)+r*cos(t)*A1(1)+r*sin(t)*B1(1);
R(2)=O(2)+r*cos(t)*A1(2)+r*sin(t)*B1(2);
R(3)=O(3)+r*cos(t)*A1(3)+r*sin(t)*B1(3);
end
function [dist]=distance1(T)
e=subs(C1,t1,T(1));
e2=subs(L_C,t3,T(2));
dist=vpa(norm((e-e2)),5);
end
function [dist]=distance2(T)
e=subs(C2,t2,T(1));
e2=subs(L_C,t3,T(2));
dist=vpa(norm(e-e2),5);
end
function R = plane(A,B,C,x,y,z)
    AB=A-B;
    BC=B-C;
    P=cross(AB,BC);
    R=P(1)*x+P(2)*y+P(3)*z-dot(P,C);
end
function [opt]=Finalpoint(A)
opt=eval(10*(13-norm(A-b(:,i)))^2+10*(13-norm(A-d))^2+5*(subs(Plane,[x;y;z],A))^2+10*(13*sqrt(2)-norm(A-pc))^2);    
end
end