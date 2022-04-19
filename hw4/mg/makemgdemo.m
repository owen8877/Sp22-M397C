% Matlab code makemgdemo.m
% For "Applied Numerical Linear Algebra",  Question 6.16
% Written by James Demmel, Jul 10, 1993
%                Modified, Jun  2, 1997
%
% Initialize b, x, jac1, jac2, iter to demonstrate testfmg
%
figure(1)
hold off, clf
k=4;
kk=2^k+1;
[x,y]=meshgrid(1:kk,1:kk);
b=sin((x-1)*3*pi/kk).*sin((y-1)*pi/kk);
b=rand(size(b))-.5;
b=b*0;
p1=((kk-1)/8:3*(kk-1)/8)+1;
p2=(5*(kk-1)/8:7*(kk-1)/8)-1;
b(p1,p2)=b(p1,p2)  -2*ones(prod(size(p1)),prod(size(p2)));
b(p2,p1)=b(p2,p1) + ones(prod(size(p2)),prod(size(p1)));
[nb,mb]=size(b);
b(:,1)=zeros(nb,1);b(:,mb)=zeros(nb,1);
b(1,:)=zeros(1,mb);b(nb,:)=zeros(1,mb);
xtrue = b;
b=zeros(size(b));
b(2:kk-1,2:kk-1) = 4*xtrue(2:kk-1,2:kk-1) ...
                    -xtrue(1:kk-2,2:kk-1) ...
                    -xtrue(3:kk  ,2:kk-1) ...
                    -xtrue(2:kk-1,1:kk-2) ...
                    -xtrue(2:kk-1,3:kk  );
hold off; clf
% x=b;b=zeros(size(b));
subplot(1,2,1); mesh(xtrue), title('true solution')
subplot(1,2,2); mesh(b), title('right hand side')
x=zeros(size(b));
jac1=1;jac2=1;iter=20;
