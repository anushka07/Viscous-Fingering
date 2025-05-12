%function misc2D_modified

%Grid NxN
N   = 256;
Pe= 500;
h1  = Pe/(N-1);
A=2;
h2=Pe/(N-1);
x= -Pe/2:h1:Pe/2;x= x(1:N);
y= -Pe/(2*A):h1:Pe/(2*A);y=y(1:N/A);
[xx,yy]= meshgrid(x,y');

%Model parameters and initial conditions
R=3.0;
t=0.00001;
[c,cs,csx]= ini_c(xx,yy,N,t,Pe);
NSTEP= 5000000;
Psi= 0*c;
ux= 0*c;
uy= 0*c;
J2=0*c;
J1=0*c;
%Wave numbers...
k1= (2*pi/Pe)*[0:(N/2-1) -(N/2):(-1)];
%k1= (2*pi/Pe)*(0:(N-1));
%k2=(2*pi*A/Pe)*(0:(N/A-1));
k2=(2*pi*A/Pe)*[0:(N/(2*A)-1) -(N/(2*A)):(-1)];
[kx,ky]= meshgrid(k1,k2);
kx2= kx.^2;
ky2= ky.^2;
k2i= -(kx.^2 + ky.^2);k2i(1,1)= 1;k2i= k2i.^-1;
%px=-mu/KP;
%px=-mu;
%py=0*c;
dt= 0.05;
J1 = J(ux,uy,c,csx,kx,ky);
%J2=J1;
H='Pe500_R3_';
path1 = ['C:\Users\ARC\Desktop\MTP\FigNew\',H,'\',H,'_jet'];
mkdir(path1);
path2 = ['C:\Users\ARC\Desktop\MTP\FigNew\',H,'\',H,'_gray'];
mkdir(path2);


% %dt= 1/(8*N);
% %t= 0;
% % %% Initialize video
% % myVideo = VideoWriter('NN3_Pe3000_R3'); %open video file
% % myVideo.FrameRate = 5;  %can adjust this, 5 - 10 works well for me
% % open(myVideo)

for istep= 1:NSTEP;
  if (istep==1)
    for istage= 1:3;
        f= f_c(c,csx,ux,uy,Pe,kx,ky,kx2,ky2);
        
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
   % J1 = J(ux,uy,c,csx,kx,ky);
  end 
    %Advance saturations
    c0= c;
    ce=cs+c;
 
    %compute N,J
     
     %N = -uy*(csx + dcx)+ux*(dcy)+dcy;
    
    %Update velocities
    [Psi,ux,uy]= streamf(ux,uy,c0,csx,R,kx,ky,k2i);
    
    %compute J(t)
    J2=J1;                             %put previous value in J2
    
    J1 = J(ux,uy,c0,csx,kx,ky);      %new value in J1
    
    %Predictor step c(t+dt)
    cp = AD(J1,J2,c,dt,kx,ky);
    
    [temp,cs,csx]= ini_c(xx,yy,N,t,Pe);
    %Predicted variables
   for j=1:3
    [Psi,ux,uy]= streamf(ux,uy,cp,csx,R,kx,ky,k2i);
    
    Jp = J(ux,uy,cp,csx,kx,ky);
    
    %Corrector step
    c=cor(Jp,J1,c0,kx,ky,dt);
    cp=c;
   end
    
    t= t+dt;
    %tc=t*Pe;
   %[temp,cs,csx]= ini_c(xx,yy,N,t,Pe);
   %dt= 1/(8*N*max(max(ux)));
    
    if (mod(istep,20)==0) | (istep==1);
        iptsetpref('ImshowBorder','tight');
        handle= figure(3);surf(xx,yy,ce,'facecolor','interp','edgecolor','none','facelighting','phong');
        view([0,0,1]);axis tight;axis equal; title(['t=  ' num2str(t)]);colormap(jet);drawnow
        pause(0.01) %Pause and grab frame
%      % if (istep==4000)
%        frame = getframe(gca,[0 2  430 214]); %get frame
%         L = frame2im(frame);        
%         fpath = 'C:\Users\ARC\Desktop\MTP\FigNew\NN2dt005\';
%         whereToStore=fullfile(fpath,['NN2' num2str(istep) '.png']);
%         imwrite(L, whereToStore);
% %         writeVideo(myVideo, frame);
    %  end
    
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
    end;
end;
%close(myVideo)

function[cp] =AD(J1,J2,c,dt,kx,ky)
chat=fft2(c);
J1hat=fft2(J1);J2hat=fft2(J2);
cp = (chat + (dt/2)*(-3*J1hat+J2hat.*exp(-(kx.^2+ky.^2)*dt))).*exp(-(kx.^2+ky.^2)*dt);

% c1 = -dt*(1.5*J1hat - 0.5*J2hat);
% cp = c1.^((kx.^2+ky.^2)*dt);
cp=real(ifft2(cp));

end

function[cc]=cor(Jp,J2,c,kx,ky,dt)
Jphat=fft2(Jp); J2hat=fft2(J2);
%cphat=fft2(cp);
chat=fft2(c);

cc = (chat.*exp(-(kx.^2+ky.^2)*dt) + (dt/2)*(-Jphat-J2hat.*exp(-(kx.^2+ky.^2)*dt)));

%cc = dt*(-(Jp+J2)/2 + (kx.^2+ky.^2).*(cp-c)/2 );

cc=real(ifft2(cc));
end

function[J1]=J(ux,uy,c,csx,kx,ky)
i= sqrt(-1);

chat= fft2(c);
dcx= real(ifft2(i*kx.*chat));
dcy= real(ifft2(i*ky.*chat));

J1=ux.*(csx + dcx)+uy.*(dcy);

end

function [c,cs,csx]= ini_c(xx,yy,N,t,Pe)

m= size(xx,1);

if N<300;delta= 10;else;delta= 30000;end
% w= [10 20 30 40 50 60];
% pos= 0*yy;
% for iw= 1:length(w);
%    pos= pos + 0.01*(rand-0.5)*cos(w(iw)*pi*yy);
% end;
% 
% arg= (xx-4*sqrt(4*t)+pos)/sqrt(4*t);



cs = 0.5*(1-erf((xx)/sqrt(4*t)));
%ce= 0.5*(1-erf(arg))-cs;
ce=0*cs;
xx2=xx.^2;
csx = -(1.0/sqrt(4*t*pi))*exp(-(xx2)/(4*t));
for j=1:1:m
     ce(j,:) = .01*rand()*exp(-(xx(j,:)/Pe).^2/(.01*.01));
end

%cs=ce+cs;
c=ce;
%N12= N/2;
 %for i= 1:N12
  %    c(:,N-(i-1))=c(:,i);
 %end



end


function f= f_c(c,csx,ux,uy,Pe,kx,ky,kx2,ky2)

i= sqrt(-1);

chat= fft2(c);
cuxhat= fft2(c.*ux);
cuyhat= fft2(c.*uy);
duxcx= real(ifft2(i*kx.*cuxhat));
duycy= real(ifft2(i*ky.*cuyhat));

dcxx= real(ifft2(-kx2.*chat));
dcyy= real(ifft2(-ky2.*chat));

f= -( duxcx + csx.*ux + duycy - (dcxx + dcyy) );

end

function [Psi,ux,uy]= streamf(ux,uy,c,csx,R,kx,ky,k2i)



i= sqrt(-1);

chat= fft2(c);

dcx= real(ifft2(i*kx.*chat));
dcy= real(ifft2(i*ky.*chat));

N = -uy.*(csx + dcx)+ux.*(dcy)+dcy;

Nhat = fft2(N);
%-Vorticity
%omega= R*( (dcx).*uy+ (csx).*uy - dcy.*ux -dcy );
omegahat = -R*Nhat;
%Solve for the streamfunction
%omegahat= fft2(omega);
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
