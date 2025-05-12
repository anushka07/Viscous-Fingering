%function misc2D_modified

%Grid NxN
N   = 512;
Pe= 1000;
h1  = Pe/(N-1);
A=2;
x= -Pe/2:h1:Pe/2;x= x(1:N);
y= -Pe/(2*A):h1:Pe/(2*A);y=y(1:N/A);
[xx,yy]= meshgrid(x,y');

%Model parameters and initial conditions
Rs = -5.7;
Rf = 4.8;
d = 10;
R=Rs+Rf;
De =0.034;
n = 0.6;
w1 = 0.375;
w2 = 0.05;
b1 = 0.1;
b2 = 0.01;
RPP = log(b2/b1);
RP = log(w2/w1);

t=0.5;
[cf,cfs,cfx,c,cs,csx]= ini_c(xx,N,t,Pe,d); % c -> cs
ce=c+cs;
ce2=cf+cfs;
Psi= 0*c;
ux= 0*c;
uy= 0*c;
uxi=0*c;
uyi=0*c;
uxt=0*c;
uyt=0*c;
w=0*c;
px=0*c;
py=0*c;
%Wave numbers...
k1= (2*pi/Pe)*[0:(N/2-1) -(N/2):(-1)];
k2=(2*pi*A/Pe)*[0:(N/(2*A)-1) -(N/(2*A)):(-1)];
[kx,ky]= meshgrid(k1,k2);
kx2= kx.^2;
ky2= ky.^2;
k2i= -(kx.^2 + ky.^2);k2i(1,1)= 1;k2i= k2i.^-1;


NSTEP= 3200000;

dt = 0.05;
%dt= 1/(8*N);

% 
% H='DD_Fig4';
% path1 = ['C:\Users\ARC\Desktop\MTP\FigNew\',H,'\',H,'_jet'];
% mkdir(path1);
% path2 = ['C:\Users\ARC\Desktop\MTP\FigNew\',H,'\',H,'_gray'];
% mkdir(path2);

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
   [mu,wi,b]=prop(ce2,ce,Rf,Rs,RP,RPP,w1,b1);
   
    for istage= 1:3;
        ff= f_c(cf,cfx,ux,uy,Pe,kx,ky,kx2,ky2,d);
        fs= f_c(c,csx,ux,uy,Pe,kx,ky,kx2,ky2,1);
        
        if istage==1;
            cf= cf + dt*( (8/15)*ff );
            c= c + dt*( (8/15)*fs );
            [pxt,pyt]=Press_p(mu,wi,b,ce,n,De,uxt,uyt,px,py,ux,uy);% calculating pxt and pyt with new py and py
            [g]= g_w(w,c0s,ce,csx,c0f,cfx,ux,uy,Rf,Rs,RP,RPP,b,wi,mu,pxt,pyt,kx,ky,uxt,uyt,De,n); %calculating dw/dt function with new pxt,pyt     
            
            px = px + dt*( (8/15)*pxt );       %1st step for RK method
            py = py + dt*( (8/15)*pyt );
            w = w + dt*( (8/15)*g );
        
            pxt1=pxt;                         %storing values before next step
            pyt1=pyt;          
     
            g1=g;
            ff1= ff;
            fs1=fs;
            
        elseif istage==2;
            cf= cf + dt*( (5/12)*ff + (-17/60)*ff1 );
            c= c + dt*( (5/12)*fs + (-17/60)*fs1 );
            [pxt,pyt]=Press_p(mu,wi,b,ce,n,De,uxt,uyt,px,py,ux,uy);           
            [g]= g_w(w,c0s,ce,csx,c0f,cfx,ux,uy,Rf,Rs,RP,RPP,b,wi,mu,pxt,pyt,kx,ky,uxt,uyt,De,n);
           
            px= px + dt*( (5/12)*pxt + (-17/60)*pxt1 );
            py= py + dt*( (5/12)*pyt + (-17/60)*pyt1 );
            w = w + dt*( (5/12)*g + (-17/60)*g1);    
            
          
            pxt1=pxt;
            pyt1=pyt;
           
            g1=g;
            ff1= ff;
            fs1=fs;
        else
            cf= cf + dt*( (3/4)*ff + (-5/12)*ff1 );
            c= c + dt*( (3/4)*fs + (-5/12)*fs1 );
            [pxt,pyt]=Press_p(mu,wi,b,ce,n,De,uxt,uyt,px,py,ux,uy);
            [g]= g_w(w,c0s,ce,csx,c0f,cfx,ux,uy,Rf,Rs,RP,RPP,b,wi,mu,pxt,pyt,kx,ky,uxt,uyt,De,n);
            
            px= px + dt*( (3/4)*pxt + (-5/12)*pxt1 );
            py= py + dt*( (3/4)*pyt + (-5/12)*pyt1 );
            w = w + dt*( (3/4)*g + (-5/12)*g1 );
                
        end;
    end;
    
    %Update velocities
    [Psi,ux,uy,uxt,uyt]= streamf(w,ux,uy,kx,ky,k2i,dt);
    
    t= t+dt;
    [temp,cfs,cfx,temp,cs,csx]= ini_c(xx,N,t,Pe,d);
    %tc=t*Pe;
    
   %dt= 1/(8*N*max(max(ux)));
    
    if (mod(istep,10)==0) | (istep==1);
        iptsetpref('ImshowBorder','tight');
        handle= figure(2);surf(xx,yy,ce,'facecolor','interp','edgecolor','none','facelighting','phong');
        view([0,0,1]);axis tight;axis equal; title(['t=  ' num2str(t)]);colormap(jet);drawnow
        pause(0.01) %Pause and grab frame
%         frame = getframe(gca,[0 2  430 214]); %get frame
%         L = frame2im(frame);        
%         whereToStore=fullfile(path1,[H num2str(istep) '.png']);
%         imwrite(L, whereToStore);
% %        writeVideo(myVideo, frame);
%         colormap(flipud(gray));
%         %fpath = ['C:\Users\ARC\Desktop\MTP\FigNew\',H,'\',H,'_gray'];
%         frame = getframe(gca,[0 2  430 214]); %get frame
%         L = frame2im(frame);        
%         whereToStore=fullfile(path2,[H num2str(istep) '.png']);
%         imwrite(L, whereToStore);
%     
%       writeVideo(myVideo, frame);
    end
end;
% close(myVideo)

function [Pxt,Pyt] = Press_p(mu,wi,b,ce,n,De,uxt,uyt,px,py,ux,uy)
v2  = (1+ux).*(1+ux)+ uy.*uy;
G   = (1+De^2.*v2).^(ce.*((n-1)/2));
Pxt =-(mu.*G.*((1+ux)+(1-b).*wi.*G.*uxt)+px)./(wi.*G);
Pyt =-(mu.*G.*(uy +(1-b).*wi.*G.*uyt)+py)./(wi.*G);
end

function [mu,wi,b]=prop(ce2,ce,Rf,Rs,RP,RPP,w1,b1)
mu = exp(Rs*(1-ce) + Rf*ce2);
wi = w1*exp(RP*(1-ce));
b = b1*exp(RPP*(1-ce));
end

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
fhat= -(Jhat + (d)*(kx2.*chat+ky2.*chat) );
f=real(ifft2(fhat));


end

function [g]= g_w(w,c,ce,csx,cf,cfsx,ux,uy,Rf,Rs,RP,RPP,b,wi,mu,pxt,pyt,kx,ky,uxt,uyt,De,n)

i= sqrt(-1);
uxhat=fft2(ux);
uyhat=fft2(uy);
uxx=real(ifft2(i*kx.*uxhat));
uxy = real(ifft2(i*ky.*uxhat));
uyy=real(ifft2(i*ky.*uyhat));
uyx = real(ifft2(i*kx.*uyhat));
chat= fft2(c);
cfhat=fft2(cf);

cx=real(ifft2(i*kx.*chat));
cy=real(ifft2(i*ky.*chat));
cfx=real(ifft2(i*kx.*cfhat));
cfy=real(ifft2(i*ky.*cfhat));
v2 = (1+ux).*(1+ux)+uy.*uy;

%G = (1+De^2.*v2).^(c.*((n-1)/2)); Had written this earlier, only
%perturbation was going into the G function
G = (1+De^2.*v2).^(ce.*((n-1)/2)); %After correction
Lx =((n-1)/2).*(cx + csx).*log(1+De^2.*v2)+((n-1).*ce.*De^2.*((ux+1).*uxx + uy.*uyx))./(1+De^2.*v2) ;
Ly =((n-1)/2).*cy.*log(1+De^2.*v2)+((n-1).*ce.*De^2.*((ux+1).*uxy + uy.*uyy))./(1+De^2.*v2); 



K=(wi.*(1-b)).^-1;

g = ((mu.*G.*(1-b)).^-1).*((pxt).*(Ly - RP.*cy)...
    + (-pyt).*(Lx - RP.*(cx+csx)))...
    + K.*( -w./G +...
    + ((ux+1)./G).*(Ly - Rs.*cy + Rf.*cfy)...
    - (uy./G).*(Lx - Rs.*(cx+csx) + Rf.*(cfx+cfsx) ))...
    + ((1-b).^-1).*(uxt.*(2*(1-b).*Ly + cy.*(-(Rs+RP).*(1-b)+RPP.*b) + cfy.*((Rf).*(1-b)+RPP.*b))...
    - uyt.*(2*(1-b).*Lx + (cx+csx).*(-(Rs+RP).*(1-b)+RPP.*b) + (cfx+cfsx).*((Rf).*(1-b)+RPP.*b)))...
    ;

end



function [Psi,ux,uy,uxt,uyt]= streamf(omega,ux,uy,kx,ky,k2i,dt)

i= sqrt(-1);

uxi=ux;
uyi=uy;


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

uxt=(ux-uxi)/dt;
uyt=(uy-uyi)/dt;


end

