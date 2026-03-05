clear variables
format long

clc
clear 
clf


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   			                                                    %
%                 Input parameters & flags                    %
%                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nnodi=512;   % number of heights per profile
nprof=nnodi; % number of profiles
flag_smoothing=0;    % 0: no smoothing; 1: smoothing applied


Ns1=10;               % number of surfaces belonging to body 1
Ns2=10;               % number of surfaces belonging to body 2
flagi1=0;             % 1: invert surfaces belonging to body 1; 0: do not invert
flagi2=0;             % 1: invert surfaces belonging to body 2; 0: do not invert

xxx1='D:\PVS\smBZM10_E_'; % path for the surfaces of the mold/stone (body 1) 	
xxx2='D:\PVS\smBZM10_A_'; % path for the surfaces of the mold/stone (body 2)


if(flag_smoothing==1)
    radius=5;        % square window for smoothing of lateral size 2*radius
else
    radius=0;
end    
pt=0.1;              % max% of z to be corrected by the smoothing algorithm based on bisection
tolp=0.2;            % tolerance in pt (% of pt)


for g=1:Ns1+Ns2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%							                                                %
%       Read the rough surface from confocal data             %
%       Data are stored in micrometers and are converted      %
%       into meters                                           %
%                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   xyz=[];
 
   if(g<=Ns1) % body 1     
     lab=sprintf('%03d.dat', g);
     xxx=strcat(xxx1,lab);	
     fp=fopen(['',xxx],'r');
     xyz=fscanf(fp,'%f',[3,inf]);    % Units of Measure: micrometers
     fclose(fp);

     x=reshape(xyz(1,:),nnodi,nnodi);    
     y=reshape(xyz(2,:),nnodi,nnodi);
     z=reshape(xyz(3,:),nnodi,nnodi); 

     x=x*1e-6; % conversion to meters
     y=y*1e-6;
     z=z*1e-6;

     side=max(max(x));       % lateral side of the rough surface
     delta=side/(nnodi-1);   % x,y grid spacing

     % Remove form via polynomial (e.g. tilt, curvature) (polyfitn)
     order    = 1;                       % remove mean plane tilt with order = 1
     Z_form   = form_poly(x,y,z,order);
     Z_noform = z - Z_form;               
     z=Z_noform;  
   
     z=z-min(min(z));

     if(flagi1==1)
       z=max(max(z))-z;
     end

   else % body 2

     lab=sprintf('%03d.dat', g-Ns1);
     xxx=strcat(xxx2,lab);      
     fp=fopen(['',xxx],'r');
     xyz=fscanf(fp,'%f',[3,inf]); % Units of Measure: micrometers
     fclose(fp);

     x=reshape(xyz(1,:),nnodi,nnodi);    
     y=reshape(xyz(2,:),nnodi,nnodi);
     z=reshape(xyz(3,:),nnodi,nnodi); 

     x=x*1e-6; % conversion to meters
     y=y*1e-6;
     z=z*1e-6;

     side=max(max(x));       % lateral side of the rough surface
     delta=side/(nnodi-1);   % x,y grid spacing     
 
     z=z-min(min(z));

     if(flagi2==1)
       z=max(max(z))-z;
     end

   end %end if

   if(flag_smoothing==1)
      z1=smoothing_function(x,y,z,radius,pt,tolp);
      z=z1;
      z1=[];

      xyz(1,:)=reshape(x,nnodi*nnodi,1);    
      xyz(2,:)=reshape(y,nnodi*nnodi,1);
      xyz(3,:)=reshape(z,nnodi*nnodi,1); 

      yyy=strcat('sm',xxx);	
      fp=fopen(['',yyy],'w');
      fprintf(fp,'%e %e %e\r\n',xyz); 
      fclose(fp);
   end  
      
   zave=mean(mean(z));
   zmax=max(max(z));
   rms=std2(z);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%  Compute profile statistics: slopes, maxima 2D (peaks) curvatures %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m0(g)=(std2(z))^2;

% rms slopes in directions x and y

slopey=[];
ny=0;
for j=2+radius:nnodi-1-radius
	 for i=2+radius:nprof-1-radius
            ny=ny+1;
    	    slopey(i-1,j-1)=(z(i+1,j)-z(i-1,j))/(2*delta);
            slopeyv(ny)=slopey(i-1,j-1);
	 end
end

slopex=[];
nx=0;
for i=2+radius:nprof-1-radius
	for j=2+radius:nnodi-1-radius
            nx=nx+1;
  	        slopex(i-1,j-1)=(z(i,j+1)-z(i,j-1))/(2*delta);
            slopexv(nx)=slopex(i-1,j-1); 
	end
end

% peak curvatures of profiles

curvp=[];
zpeak=[];
n_peaks=0;
for j=2+radius:nnodi-1-radius
	for i=2+radius:nprof-1-radius
    	if (z(i,j)>z(i,j-1) && z(i,j)>z(i,j+1))
      	   n_peaks=n_peaks+1;
       	   curvp(n_peaks)=-(z(i,j+1)-2*z(i,j)+z(i,j-1))/(delta^2);
           zpeak(n_peaks)=z(i,j);
    	end
    end	
end

% statistics of profile slopes
rms_slopex(g)=std(slopexv);
rms_slopey(g)=std(slopeyv);

m2x(g)=rms_slopex(g)^2;
m2y(g)=rms_slopey(g)^2;

% statistics of the peak heights (2D maxima)
mean_z_peaks(g)=mean(zpeak);
rms_z_peaks(g)=std(zpeak);

% statistics of the peak curvatures (2D maxima)
mean_curv_peaks(g)=mean(curvp);
rms_curv_peaks(g)=std(curvp);

m4(g)=rms_curv_peaks(g)^2;

density_peaks(g)=n_peaks/((nnodi-2*radius)*(nprof-2*radius)); %density of peaks

alfa_x(g)=m0(g)*m4(g)/m2x(g)^2;
alfa_y(g)=m0(g)*m4(g)/m2y(g)^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%             Compute asperity (3D maxima) statistics               %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Curvx=[];
Curvy=[];
R=[];

for i=2+radius:nprof-1-radius
   for j=2+radius:nnodi-1-radius
       if (z(i,j)>z(i,j-1) && z(i,j)>z(i,j+1))
       Curvy(i,j)=-2*(-delta*z(i,j-1)+2*delta*z(i,j)-delta*z(i,j+1))/(-delta*y(i,j-1)^2+2*delta*y(i,j)^2-delta*y(i,j+1)^2);
       else
       Curvy(i,j)=0;
       end
   end
end

for j=2+radius:nnodi-1-radius
   for i=2+radius:nprof-1-radius
       if (z(i,j)>z(i-1,j) && z(i,j)>z(i+1,j))
       Curvx(i,j)=-2*(-delta*z(i-1,j)+2*delta*z(i,j)-delta*z(i+1,j))/(-delta*x(i-1,j)^2+2*delta*x(i,j)^2-delta*x(i+1,j)^2);
       else
       Curvx(i,j)=0;
       end
   end
end

stdcx=std2(Curvx);
for i=2+radius:nprof-1-radius
   for j=2+radius:nnodi-1-radius
      if(Curvx(i,j)/stdcx>5 || Curvx(i,j)/stdcx<0.1)
          Curvx(i,j)=0;
      end
   end
end

stdcy=std2(Curvy);
for j=2+radius:nnodi-1-radius
   for i=2+radius:nprof-1-radius
      if(Curvy(i,j)/stdcy>5 || Curvy(i,j)/stdcy<0.1)
          Curvy(i,j)=0;
      end
   end
end


ns=0;

H=[];
ro=[];
curv=[];

for i=2+radius:nprof-1-radius
   for j=2+radius:nnodi-1-radius
       if (Curvx(i,j)*Curvy(i,j)~=0)
       R(i,j)=1/sqrt(Curvx(i,j)*Curvy(i,j));
       ns=ns+1;
       H(ns)=z(i,j);
       ro(ns)=R(i,j);
       curv(ns)=1/R(i,j);
       else
       R(i,j)=0;
       end
   end
end

% statistics of asperity (3D maxima) heights

mean_z_asperities(g)=mean(H);
rms_z_asperities(g)=std(H);

% statistics of asperity (3D maxima) curvatures

mean_curv_asperities(g)=mean(curv(1:ns));
rms_curv_asperities(g)=std(curv(1:ns));

density_asperities(g)=ns/((nnodi-2*radius)*(nprof-2*radius));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     %
%               Histograms of slopes                  %
%                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure(5*(g-1)+2)
%hist(slopexv,100);

%figure(5*(g-1)+3)
%hist(slopeyv,100);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                      %
% Joint probability of asperity heights and curvatures %
%                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

norm_asp=(H-mean_z_asperities(g))/rms_z_asperities(g); 
norm_curv=curv/rms_curv_asperities(g);
norm_asp=norm_asp';
norm_curv=norm_curv'; 

P=[norm_asp norm_curv];
num=[15 15]; 
X=linspace(min(norm_asp),max(norm_asp),num(1));
Y=linspace(min(norm_curv),max(norm_curv),num(2));                                               
iso=hist3(P,{X,Y});
iso=iso';
tot=sum(sum(iso));
iso=iso/(tot*(X(2)-X(1))*(Y(2)-Y(1)));       % joint probability


figure(g)
subplot(2,2,2);
[mC hC] = contour(X,Y,iso);
hold on
clabel(mC,hC)
xlabel('Normalized asperity height','fontsize',20);
ylabel('Dimensionless asperity curvature','fontsize',20);
h1 = gca;
fontsize(h1,20,'points')
set(gcf,'color','w');

subplot(2,2,4);
[n1,ctr1] = hist(norm_asp,50);
bar(ctr1,n1/(sum(n1)*(ctr1(2)-ctr1(1))),1);
xlim([min(norm_asp) max(norm_asp)]);
h2 = gca;
fontsize(h2,20,'points')
set(gcf,'color','w');

subplot(2,2,1);
[n2,ctr2] = hist(norm_curv,50);
barh(ctr2,n2/(sum(n2)*(ctr2(2)-ctr2(1))),1);
ylim([min(norm_curv) max(norm_curv)]);
h3 = gca;
fontsize(h3,20,'points')
set(gcf,'color','w');

subplot(2,2,3);
mesh(x,y,z)
xlabel('x (m)','FontSize',15);
ylabel('y (m)','FontSize',15);
zlabel('z (m)','FontSize',15);
h4 = gca;
fontsize(h4,20,'points')
set(gcf,'color','w');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%					                                     %
%                    PSD                       %
%                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[PSD] = psd_2D(z, delta);

qx(:,g)=PSD.qx;
qy(:,g)=PSD.qy;
Cx(:,g)=PSD.Cx;
Cy(:,g)=PSD.Cy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%					                                     %
%         Save statistical parameters          %
%                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ris_stat=[m0',rms_slopex',rms_slopey',mean_z_peaks',rms_z_peaks',mean_curv_peaks',rms_curv_peaks',density_peaks',mean_z_asperities',rms_z_asperities',mean_curv_asperities',rms_curv_asperities',density_asperities',alfa_x',alfa_y'];
save ris_stat.dat ris_stat -ASCII -DOUBLE -TABS


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                              %
%          Plot results for comparisons        %
%                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(g+1)
plot(m0(1:Ns1),'k')
hold on
plot(m0(Ns1+1:Ns1+Ns2),'k--')
xlabel('Surface','FontSize',20);
ylabel('m_0 (m^2)','FontSize',20);
h = gca;
fontsize(h,20,'points')
set(gcf,'color','w');

figure(g+2)
plot((rms_slopex(1:Ns1)).^2,'k')
hold on
plot((rms_slopey(1:Ns1)).^2,'b')
plot((rms_slopex(Ns1+1:Ns1+Ns2)).^2,'k--')
plot((rms_slopey(Ns1+1:Ns1+Ns2)).^2,'b--')
xlabel('Surface','FontSize',20);
ylabel('\sigma_{m}','FontSize',20);
h = gca;
fontsize(h,20,'points')
set(gcf,'color','w');

figure(g+3)
plot((rms_curv_asperities(1:Ns1)).^2,'k')
hold on
plot((rms_curv_asperities(Ns1+1:Ns1+Ns2)).^2,'k--')
xlabel('Surface','FontSize',20);
ylabel('\sigma_{\kappa}','FontSize',20);
h = gca;
fontsize(h,20,'points')
set(gcf,'color','w');

% Power spectra in direction x and y

figure(g+4)
loglog(qx(:,1:Ns1),Cx(:,1:Ns1),'k') % surfaces body 1
hold on
loglog(qx(:,Ns1+1:Ns1+Ns2),Cx(:,Ns1+1:Ns1+Ns2),'r') % surfaces body 2
loglog(qy(:,1:Ns1),Cy(:,1:Ns1),'k--') % surfaces body 1
loglog(qy(:,Ns1+1:Ns1+Ns2),Cy(:,Ns1+1:Ns1+Ns2),'r--') % surfaces body 2
xlabel('q (m^{-1})','FontSize',20);
ylabel('C (m^4)','FontSize',20);
h = gca;
fontsize(h,20,'points')
set(gcf,'color','w');

end % end of cycle over g surfaces

