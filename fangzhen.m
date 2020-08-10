clear all;
Ai1=1.56;Ai2=3; Ap1=1.56;Ap2=3;
h1=0.25;%  
s=h1/2;%%
k=0.0001;%time step

x=0:s:1;%L=1

r=k/s/s;%步长比
r1=1-2*r;
t=0:k:10;

m=length(x);

n=length(t);



L1= 29.3434;L2=27.7151;%%observer gain
K1= 31.4166; K2=31.0823;%6control gain 
U=zeros(m,n);%%初始化解矩阵y(x,t) without control

for i=1:m
      U(i,1)=0;
end


for j=1:n-1
        U(1,j)=0;
    U(m,j)=U(m-1,j);
     for i=2:m-1
         U(i,j+1)=r*U(i-1,j)+r1*U(i,j)+r*U(i+1,j)+k*(3*U(i,j)-U(i,j)^3+(3+i)*exp(-0.002*j));
     end
     
end

figure(1)
surf(t,x,U);
shading interp;
% title('open-loop of $ y(x,t)$','Interpreter','latex','Fontsize',12);
xlabel('$t$','Interpreter','latex','Fontsize',12);
ylabel('$x$','Interpreter','latex','Fontsize',12);
zlabel('$z$','Interpreter','latex','Rotation',0,'Fontsize',12);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U=zeros(m,n);%%初始化解矩阵y(x,t) under control
hh1=zeros(m,n);hh2=zeros(m,n);%%%hp(\phi(x,t))

u=zeros(m,n);%%初始化解矩阵\hat y(x,t) 
H1=zeros(m,n);H2=zeros(m,n);
for i=1:m
    U(i,1)=0;   
     u(i,1)=0;   
end
for j=1:n-1%shijian  t
   U(1,j)=0; 
 U(m,j)=U(m-1,j); 
 u(1,j)=0;
 u(m,j)= u(m-1,j);
 for i=2:(m-1)%%%weizhix, t=0, \hat{y}(\bar{x}j,tk)-y(\bar{x}j,tk)=0
      H1(i,j)=u(i,j)*u(i,j)/1.44;H2(i,j)=1-H1(i,j);%%h1(\hat{\phi}(x,t)),h2(\hat{\phi}(x,t))
  u(i,j+1)=r1*u(i,j)+r*(u(i-1,j)+u(i+1,j))+k*(H1(i,j)*(Ai1*u(i,j)-K1*u(2*fix(i/2),j)-L1*(u(2*fix(i/2),j)-U(2*fix(i/2),j)))+H2(i,j)*(Ai2*u(i,j)-K2*u(2*fix(i/2),j)-L2*(u(2*fix(i/2),j)-U(2*fix(i/2),j))));
              hh1(i,j)=U(i,j)*U(i,j)/1.44;hh2(i,j)=1-hh1(i,j);%%h1({\phi}(x,t)),h2({\phi}(x,t))     
  U(i,j+1)=r1*U(i,j)+r*(U(i-1,j)+U(i+1,j))+k*(hh1(i,j)*Ap1*U(i,j)+hh2(i,j)*Ap2*U(i,j)-H1(i,j)*K1*u(2*fix(i/2),j)-H2(i,j)*K2*u(2*fix(i/2),j)+(3+i)*exp(-0.002*j));

 end  
end




figure(3)
surf(t,x,U);
shading interp;
xlabel('$t$','Interpreter','latex','Fontsize',10);
ylabel('$x$','Interpreter','latex','Fontsize',10);
zlabel('$z$','Interpreter','latex','Rotation',0,'Fontsize',10);
figure(4)
surf(t,x,u);
shading interp;
xlabel('$t$','Interpreter','latex','Fontsize',10);
ylabel('$x$','Interpreter','latex','Fontsize',10);
zlabel('$\hat{z}$','Interpreter','latex','Rotation',0,'Fontsize',10);

%%%%%%J(T)%%%%%%%%%%%%%%%%%%%
 S=zeros(1,m);
 J=zeros(1,n);

    for j=1:n-1
    for i=1:m
        S(i+1)=S(i)+((u(i,j)-U(i,j))^2*2.5^2-0.19*0.19*((3+i)*exp(-0.002*j))^2)*s;
    end
    J(j+1)=J(j)+S(end)*k;
    end
figure(7)
plot(t,J);
xlabel('$T$','Interpreter','latex','Fontsize',12);
ylabel('$J(T)$','Interpreter','latex','Rotation',0,'Fontsize',10);
xlim([0,2])
