function [z1] = smoothing_function(x,y,z,radius,pt,tol)

N=size(z,1);


z1(1:N,1:N)=z;
count=0;
f1=0.1;
for i=1+radius:N-radius
    for j=1+radius:N-radius
        zave=mean2(z(i-radius:i+radius,j-radius:j+radius));
        rms=std2(z(i-radius:i+radius,j-radius:j+radius));
        if( (z(i,j)>zave+f1*rms) || (z(i,j)<zave-f1*rms))
           z1(i,j)=zave;
           count=count+1;
        end
    end    
end
p1=count/(N*N);


z1(1:N,1:N)=z;
count=0;
f2=3;
for i=1+radius:N-radius
    for j=1+radius:N-radius
        zave=mean2(z(i-radius:i+radius,j-radius:j+radius));
        rms=std2(z(i-radius:i+radius,j-radius:j+radius));
        if( (z(i,j)>zave+f2*rms) || (z(i,j)<zave-f2*rms))
           z1(i,j)=zave;
           count=count+1;
        end
    end    
end
p2=count/(N*N);

k=1;
err(k)=(p2-pt)/pt;

while (abs(err(k))>tol && k<20)
    k=k+1;
    z1(1:N,1:N)=z;
    f(k)=f1+(p1-pt)*(f2-f1)/(p1-p2);
    count=0;
    for i=1+radius:N-radius
       for j=1+radius:N-radius
        zave=mean2(z(i-radius:i+radius,j-radius:j+radius));
        rms=std2(z(i-radius:i+radius,j-radius:j+radius));
        if( (z(i,j)>zave+f(k)*rms) || (z(i,j)<zave-f(k)*rms))
           z1(i,j)=zave;
           count=count+1;
        end
       end  
    end
    p(k)=count/(N*N);
    err(k)=(p(k)-pt)/pt;
    if(p(k)>pt)
       f1=f(k);
    else
       f2=f(k);
    end
end

z1(1:radius,:)=0;
z1(N-radius:N,:)=0;
z1(:,1:radius)=0;
z1(:,N-radius:N)=0;


% Activate for visualization purposes 
% 
% figure(1) % original surface
% mesh(x,y,z)
% axis equal
% 
% figure(2) % magnitude of the applied corrections
% mesh(x,y,z1-z)
% axis equal
% 
% figure(3) % modified surface
% mesh(x,y,z1)
% axis equal





