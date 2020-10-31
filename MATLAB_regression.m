clc;
clear all;
close all;

raw=input('Input the file name :','s');
fs=input(' F threshold value for identifying variables (e.g., 0.1, 0.05): ');
x=(xlsread(raw));
m=size(x,2); n=size(x,1);
xb=zeros(1,n); sg=zeros(1,n); ds=zeros(1,n);
a=zeros(n);
iss='';
for i=1:n
c=0;

for j=1:m
c=c+x(i,j);
end
xb(i)=c/m;
c=0;
for j=1:m
c=c+(x(i,j)-xb(i))^2;
end
sg(i)=sqrt(c);
end
h=sg(n);
for i=1:n-1
for j=i+1:n
c=0;
for k=1:m
c=c+(x(i,k)-xb(i))*(x(j,k)-xb(j));
end
a(i,j)=c/(sg(i)*sg(j)); a(j,i)=a(i,j);
end; end
for i=1:n
xb(i)=i; sg(i)=0; a(i,i)=1;
end
disp('Correlation matrix')
CorrMat=a
l=0; s=0;
while (n>=1)
if (l==n-1) break; end
ma=0;
for i=1:n
ds(i)=xb(i);
end
for i=1:n-1
if (ds(i)==0) continue; end
if (a(i,i)<1e-05) continue; end
v1=a(i,n)*a(n,i)/a(i,i);
if (v1>ma) ma=v1; k=i; end
end
f1=ma*(m-l-2)/(a(n,n)-ma);
if (f1<=fs) break; end
xb(k)=0; sg(k)=k;
l=l+1;
for i=1:n
for j=1:n
if ((i~=k) & (j~=k)) a(i,j)=a(i,j)-a(i,k)*a(k,j)/a(k,k); end
end; end

for j=1:n
if (j~=k) a(k,j)=a(k,j)/a(k,k); a(j,k)=-a(j,k)/a(k,k); end
end
a(k,k)=1/a(k,k);
r=sqrt(1-a(n,n));
yn=h*sqrt(a(n,n)/(m-l-1));
if (s==0) s=1; continue; end
lab=0;
while (n>=1)
ma=-1e+18;
for i=1:n
ds(i)=sg(i);
end
for i=1:n-1
if (ds(i)==0) continue; end
if (a(i,i)<1e-05) continue; end
v1=a(i,n)*a(n,i)/a(i,i);
if (v1>ma) ma=v1;k=i; end
end
f1=-ma*(m-l-1)/a(n,n);
if (f1>fs) lab=1; break; end
sg(k)=0; xb(k)=k;
l=l-1;
for i=1:n
for j=1:n
if ((i~=k) & (j~=k)) a(i,j)=a(i,j)-a(i,k)*a(k,j)/a(k,k); end
end; end
for j=1:n
if (j~=k) a(k,j)=a(k,j)/a(k,k); a(j,k)=-a(j,k)/a(k,k); end
end
a(k,k)=1/a(k,k);
r=sqrt(1-a(n,n));
yn=h*sqrt(a(n,n)/(m-l-1));
end;
if (lab==1) continue; end
end
for i=1:n-1
a(i,1)=sg(i);
end
for i=1:n
c=0;
for j=1:m
c=c+x(i,j);
end
xb(i)=c/m;

c=0;
for j=1:m
c=c+(x(i,j)-xb(i))^2;
end
sg(i)=sqrt(c);
end
h=sg(n);
c=0;
for i=1:n-1
if (a(i,1)==0) continue; end
ds(i)=a(i,n)*sg(n)/sg(i);
a(i,2)=ds(i);
c=c+ds(i)*xb(i);
end
s=xb(n)-c;
iss=strcat(iss,'Qualified variables: \n');
for i=1:n-1
if (a(i,1)==0) continue; end
if (ds(i)~=0) iss=strcat(iss,' X-',num2str(i)); end
if ((ds(i+1)~=0) & (i<n-1)) iss=strcat(iss,','); end
if (ds(i)~=0)
end; end
iss=strcat(iss,'\nStepwise regression equation:\n');
iss=strcat(iss,'Y = ',num2str(s));
for i=1:n-1
if (a(i,1)==0) continue; end
if (ds(i)>0) e1=num2str(ds(i)); end
if (ds(i)<0) e1=num2str(abs(ds(i))); end
if (ds(i)>0) iss=strcat(iss,' + ',e1,' X (',num2str(i),')'); end
if (ds(i)<0) iss=strcat(iss,' - ',e1,' X (',num2str(i),')'); end
end
iss=strcat(iss,'\nCorrelation coefficient R=',num2str(r),', ','F value=',num2str(fs),'\n');
fprintf(iss) 