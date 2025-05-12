%function misc2D_modified

%Grid NxN
N   = 4096;
Pe= 1000;
h1  = Pe/(N-1);
A=2;
h2=Pe/(N-1);
x= -Pe/2:h1:Pe/2;x= x(1:N);
y= -Pe/(2*A):h1:Pe/(2*A);y=y(1:N/A);
[xx,yy]= meshgrid(x,y');

%Model parameters and initial conditions
Rs=-6.1;
Rf= 6.0;
d=10;
R=Rs+Rf;
t=0.1;
[cf,cfs,cfx,c,cs,csx]= ini_c(xx,N,t,Pe,d);
ce=c+cs;
ce2=cf+cfs;
NSTEP= 3200000;
Psi= 0*cf;
ux= 0*cf;
uy= 0*cf;

%Wave numbers...
k1= (2*pi/Pe)*[0:(N/2-1) -(N/2):(-1)];
k2=(2*pi*A/Pe)*[0:(N/(2*A)-1) -(N/(2*A)):(-1)];
[kx,ky]= meshgrid(k1,k2);
kx2= kx.^2;
ky2= ky.^2;
k2i= -(kx.^2 + ky.^2);k2i(1,1)= 1;k2i= k2i.^-1;

dt= 0.0005;


%dt= 1/(8*N);

H='DD_Fig4';
path1 = ['C:\Users\ARC\Desktop\MTP\FigNew\',H,'\',H,'_jet'];
mkdir(path1);
path2 = ['C:\Users\ARC\Desktop\MTP\FigNew\',H,'\',H,'_gray'];
mkdir(path2);

% %% Initialize video
% myVideo = VideoWriter('NN1_2Pe3000_R0DD'); %open video file
% myVideo.FrameRate = 5;  %can adjust this, 5 - 10 works well for me
% open(myVideo)

for istep= 1:NSTEP;
    
    %Advance saturations
    c0f= cf;
    ce2 = cf+cfs;
    c0s=c;
    ce = cs+c;
    %mu = exp(Rs*ce+Rf*ce2);
    %figure(2);
    %surf(xx,yy,mu,'facecolor','interp','edgecolor','none','facelighting','phong');
    for istage= 1:3;
        ff= f_c(cf,cfx,ux,uy,Pe,kx,ky,kx2,ky2,d);
        fs= f_c(c,csx,ux,uy,Pe,kx,ky,kx2,ky2,1);
        
        if istage==1;
            cf= cf + dt*( (8/15)*ff );
            c= c + dt*( (8/15)*fs );
            ff1= ff;
            fs1=fs;
            
        elseif istage==2;
            cf= cf + dt*( (5/12)*ff + (-17/60)*ff1 );
            c= c + dt*( (5/12)*fs + (-17/60)*fs1 );
            ff1= ff;
            fs1=fs;
        else
            cf= cf + dt*( (3/4)*ff + (-5/12)*ff1 );
            c= c + dt*( (3/4)*fs + (-5/12)*fs1 );
        end;
    end;
    
    %Update velocities
    [Psi,ux,uy]= streamf(N,ux,uy,c,cf,csx,cfx,Rs,Rf,kx,ky,k2i);
    
    t= t+dt;
    [temp,cfs,cfx,temp,cs,csx]= ini_c(xx,N,t,Pe,d);
    %tc=t*Pe;
    
   %dt= 1/(8*N*max(max(ux)));
    
    if (mod(istep,10)==0) | (istep==1);
        iptsetpref('ImshowBorder','tight');
        handle= figure(1);surf(xx,yy,ce,'facecolor','interp','edgecolor','none','facelighting','phong');
        view([0,0,1]);axis tight;axis equal; title(['t=  ' num2str(t)]);colormap(jet);drawnow
        pause(0.01) %Pause and grab frame
        frame = getframe(gca,[0 2  430 214]); %get frame
        L = frame2im(frame);        
        whereToStore=fullfile(path1,[H num2str(istep) '.png']);
        imwrite(L, whereToStore);
%        writeVideo(myVideo, frame);
        colormap(flipud(gray));
        %fpath = ['C:\Users\ARC\Desktop\MTP\FigNew\',H,'\',H,'_gray'];
        frame = getframe(gca,[0 2  430 214]); %get frame
        L = frame2im(frame);        
        whereToStore=fullfile(path2,[H num2str(istep) '.png']);
        imwrite(L, whereToStore);
%       writeVideo(myVideo, frame);
    end
end;
% close(myVideo)



function [cf,cfs,cfx,c,cs,csx]= ini_c(xx,N,t,Pe,d)

m= size(xx,1);
if N<300;delta= 2000000;else;delta= 30000;end


cs = 0.5*(1- erf((xx)/sqrt(4*t)));
ce=0*cs;
csx = -(1.0/sqrt(4*t*pi))*exp((-xx.*xx)/(4*t));
cfs = 0.5*(1- erf((-xx)/sqrt(4*t*d)));
cfx = (1.0/sqrt(4*t*pi*d))*exp((-xx.*xx)/(4*t*d));

for j=1:1:m
    ce(j,:) = 0.01*rand()*exp(-(xx(j,:)/Pe).^2/(.01*.01));
end

c=ce;
cf = ce;
%N12= N/2;
 %for i= 1:N12
  %    c(:,N-(i-1))=c(:,i);
 %end



end


function f= f_c(c,csx,ux,uy,Pe,kx,ky,kx2,ky2,d)

i= sqrt(-1);

chat= fft2(c);
cuxhat= fft2(c.*(ux));
cuyhat= fft2(c.*uy);
duxcx= real(ifft2(i*kx.*cuxhat));
duycy= real(ifft2(i*ky.*cuyhat));
J = duxcx + duycy + ux.*csx;
Jhat = fft2(J);

dcxx= real(ifft2(-kx2.*chat));
dcyy= real(ifft2(-ky2.*chat));

%f= -( duxcx + duycy + ux.*csx - d*(dcxx + dcyy) );
fhat= -(Jhat + d*(kx2.*chat+ky2.*chat) );
f=real(ifft2(fhat));


end

function [Psi,ux,uy]= streamf(N,ux,uy,cs,cf,csx,cfx,Rs,Rf,kx,ky,k2i)

i= sqrt(-1);

cshat= fft2(cs);
cfhat = fft2(cf);

dcsx= real(ifft2(i*kx.*cshat));
dcsy= real(ifft2(i*ky.*cshat));

dcfx= real(ifft2(i*kx.*cfhat));
dcfy= real(ifft2(i*ky.*cfhat));

%-Vorticity
omega= -Rf*( dcfx.*uy  + cfx.*uy - dcfy.*ux -dcfy ) + Rs*( dcsx.*uy +csx.*uy - dcsy.*ux -dcsy );

%Solve for the streamfunction
omegahat= fft2(omega);
Psi= real(ifft2(-omegahat.*k2i));
Psi= Psi-Psi(1,1);

%Derivatives of the stream function
Psihat= fft2(Psi);
dPsix= real(ifft2(i*kx.*Psihat));
dPsiy= real(ifft2(i*ky.*Psihat));

%Velocities
ux=  dPsiy;
uy= -dPsix;
end
