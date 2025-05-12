%function misc2D_modified

%Grid NxN
N   = 512;
h1  = 1/(N-1);
A=2;
x= -1/2:h1:1/2;x= x(1:N);
y= -1/(2*A):h1:1/(2*A);y=y(1:N/A);
[xx,yy]= meshgrid(x,y');

%Model parameters and initial conditions
Pe=1500;
n=0.2;
De=1.0;
R=2.77725;
t=0.00000125;
[c,cs,csx]= ini_c(xx,yy,t,Pe);
NSTEP= 100000;
Psi= 0*c;
ux= 0*c;
uy= 0*c;

%Wave numbers...
k1= (2*pi)*[0:(N/2-1) -(N/2):(-1)];
k2=(2*pi*A)*[0:(N/(2*A)-1) -(N/(2*A)):(-1)];
[kx,ky]= meshgrid(k1,k2);
kx2= kx.^2;
ky2= ky.^2;
k2i= -(kx.^2 + ky.^2);k2i(1,1)= 1;k2i= k2i.^-1;

dt= 1/(8*N);

H='NP02';
path1 = ['C:\Users\ARC\Desktop\MTP\FigNew\Reff_Shearcode\',H,'\',H,'_jet'];
mkdir(path1);
path2 = ['C:\Users\ARC\Desktop\MTP\FigNew\Reff_Shearcode\',H,'\',H,'_gray'];
mkdir(path2);



 %% Initialize video
% myVideo = VideoWriter('myVideoFile'); %open video file
% myVideo.FrameRate = 15;  %can adjust this, 5 - 10 works well for me
% open(myVideo)

for istep= 1:NSTEP;
    
    %Advance saturations
    c0= c;
    ce=cs+c;
    for istage= 1:3;
        f= f_c(c,cs,csx,ux,uy,Pe,kx,ky,kx2,ky2);
        
        if istage==1;
            c= c + dt*( (8/15)*f );
            f1= f;
        elseif istage==2;
            c= c + dt*( (5/12)*f + (-17/60)*f1 );
            f1= f;
        else
            c= c + dt*( (3/4)*f + (-5/12)*f1 );
        end;
    end;
    
    %Update velocities
    [Psi,ux,uy]= streamf(De,n,ux,uy,c,csx,ce,R,kx,ky,k2i);
    
     t= t+dt;
    [temp,cs,csx]= ini_c(xx,yy,t,Pe);
    %tc=t*Pe;
    
   %dt= 1/(8*N*max(max(ux)));
    
    if (mod(istep,50)==0) || (istep==1)
        iptsetpref('ImshowBorder','tight');
        handle= figure(3);surf(xx,yy,ce,'facecolor','interp','edgecolor','none','facelighting','phong');
        view([0,0,1]);axis tight equal;title(['t=  ' num2str(t)]);colormap(jet);drawnow;
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
        
%       if (istep==1024)||(istep==410)||(istep==615)||(istep==820)
%           figure(2);
%         contourf(ce,0.1:0.1:0.6,'fill','off')  
%         axis equal;
%        %frame = getframe(gca,[0 2  430 214]); %get frame
%        frame = getframe(gca);
%         L = frame2im(frame);        
%         fpath = 'C:\Users\ARC\Desktop\MTP\List of Figures\Shearthinning\';
%         whereToStore=fullfile(fpath,['braj_' num2str(istep) '.png']);
%         imwrite(L, whereToStore);
%     end
    end;
end;
% close(myVideo)



function [c,cs,csx]= ini_c(xx,yy,t,Pe)

N= size(xx,1);
if N<300;delta= 2000000;else;delta= 30000;end
% pos= 0*yy;
% w= [10 20 30 40 50 60];
% for iw= 1:length(w);
%    pos= pos + 0.004*(rand-0.5)*cos(w(iw)*pi*yy);
% end;
% 
% arg= (-0.71+xx-4/sqrt(delta)+pos)*sqrt(delta);
% c= 0.5*(1-erf(arg));

cs = 0.5*(1-erf((xx)/sqrt(4*t/Pe)));
csx = -(1.0/sqrt((4*t*pi/Pe)))*exp(-(xx.*xx)/(4*t/Pe));
%csxx = 0.5/(t*sqrt(4*pi*t)).*xx.*exp(-(xx.*xx)/(4*t));
ce=0*cs;
 for j=1:1:N
    ce(j,:) = 0.01*rand()*exp(-(xx(j,:)).^2/(0.01*0.01));
 end

c=ce;

end


function f= f_c(c,cs,csx,ux,uy,Pe,kx,ky,kx2,ky2)

i= sqrt(-1);

chat= fft2(c);
cuxhat= fft2(c.*ux);
cuyhat= fft2(c.*uy);
duxcx= real(ifft2(i*kx.*cuxhat));
duycy= real(ifft2(i*ky.*cuyhat));

dcxx= real(ifft2(-kx2.*chat));
dcyy= real(ifft2(-ky2.*chat));

f= -( duxcx + duycy + (ux).*csx - (1/Pe)*(dcxx + dcyy) );
end

function [Psi,ux,uy]= streamf(De,n,ux,uy,c,csx,ce,R,kx,ky,k2i)

i= sqrt(-1);
v2 = (1+ux).*(1+ux)+uy.*uy;
chat= fft2(c);
uxhat=fft2(ux);
uyhat=fft2(uy);

dcx= real(ifft2(i*kx.*chat));
dcy= real(ifft2(i*ky.*chat));
duxx = real(ifft2(i*kx.*uxhat));
duyx=real(ifft2(i*kx.*uyhat));
duxy=real(ifft2(i*ky.*uxhat));

%-Vorticity
%omega= R*( dcx.*uy - dcy.*ux -dcy );
N=( -( dcx.*uy + csx.*uy - dcy.*ux -dcy ));
M = ((ux+1).*(ux+1)).*duxy - uy.*uy.*duyx - 2*(ux+1).*uy.*duxx;

l=((n-1)/2)*(log(1+De^2.*v2));
m=(1-ce)*(n-1).*(De*De).*M./(1+De^2.*v2);

omega = -(R + ((n-1)/2)*(log(1+De^2.*v2))).*N + (1-ce)*(n-1).*(De*De).*M./(1+De^2.*v2);

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
