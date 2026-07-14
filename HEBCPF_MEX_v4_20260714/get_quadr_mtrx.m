function [Mp,Mq]=get_quadr_mtrx(Y,n)

% get the quadratic matrix of each bus power balancing equation

% get G and B
Y=sparse(Y);
G=(Y+Y')/2;
B=(Y'-Y)/2*1i;

% find Mp and Mq
Mp1=zeros(2*n,2*n,n);
Mq1=zeros(2*n,2*n,n);
Mp=cell(n,1);
Mq=cell(n,1);
for j=1:1:n
    m1=zeros(2*n,2*n);
    m1(j,1:n)=G(j,:);
    m1(j,(n+1):2*n)=-B(j,:);
    m1(j+n,1:n)=B(j,:);
    m1(j+n,(n+1):2*n)=G(j,:);
    Mp1(:,:,j)=0.5*(m1+m1.');
    Mp{j}=sparse(Mp1(:,:,j));
    n1=zeros(2*n,2*n);
    n1(j,1:n)=-B(j,:);
    n1(j,(n+1):2*n)=-G(j,:);
    n1(j+n,1:n)=G(j,:);
    n1(j+n,(n+1):2*n)=-B(j,:);
    Mq1(:,:,j)=0.5*(n1+n1.');
    Mq{j}=sparse(Mq1(:,:,j));
end

end
