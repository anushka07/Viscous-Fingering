%function misc2D_VE

%Grid NxN
N   = 512;
Pe= 5000;
%h   = 1/N;
%x= -0.5:h:0.5;x= x(1:N);
h1  = Pe/(N-1);
A=2;
x= -Pe/2:h1:Pe/2;x= x(1:N);
y= -Pe/(2*A):h1:Pe/(2*A);y=y(1:N/A);
[xx,yy]= meshgrid(x,y');

%Model parameters and initial conditions
R= 2.0;
De =1.0;
n=0.2;
w1 = 75.0;
w2=75.0;
b1 = 0.1;
b2 = 0.1;
NSTEP= 60000;
RPP = log(b2/b1);
RP = log(w2/w1);
KP=1;
t=0.1;
[c,cs,csx]= ini_c(xx,yy,t,Pe);
ce=c+cs;
[mu,wi,b]=prop(ce,R,RP,RPP,w1,b1);
Psi= 0*c;
ux= 0*c;
uy= 0*c;
uxi=0*c;
uyi=0*c;
uxt=0*c;
uyt=0*c;
w=0*c;
%T=0.000125;

%Wave numbers...
% k= (2*pi)*[0:(N/2-1) (-N/2):(-1)];
% [kx,ky]= meshgrid(k,k);
k1= (2*pi/Pe)*[0:(N/2-1) -(N/2):(-1)];
k2=(2*pi*A/Pe)*[0:(N/(2*A)-1) -(N/(2*A)):(-1)];
[kx,ky]= meshgrid(k1,k2);
kx2= kx.^2;
ky2= ky.^2;
k2i= -(kx.^2 + ky.^2);k2i(1,1)= 1;k2i= k2i.^-1;
px=0*c;
py=0*c;
%P1y=0*c;
dt=0.05;
%t= 0;
% H='nntrial2R_3';
% path1 = ['C:\Users\ARC\Desktop\MTP\FigNew\',H,'\',H,'_jet'];
% mkdir(path1);
% path2 = ['C:\Users\ARC\Desktop\MTP\FigNew\',H,'\',H,'_gray'];
% mkdir(path2);
% %% Initialize video
% myVideo = VideoWriter('PN3_Pe5000'); %open video file
% myVideo.FrameRate = 5;  %can adjust this, 5 - 10 works well for me
% open(myVideo)

[pxt,pyt]=Press_p(mu,wi,b,ce,n,De,uxt,uyt,px,py,ux,uy);
J1 = J(ux,uy,c,csx,kx,ky);
[F1x,F1y]=PJ1(ce,mu,wi,b,ux,uy,uxt,uyt,px,py,De,n);
[D1]=WJ1(c,cs,csx,w,mu,wi,b,pxt,pyt,ux,uy,uxt,uyt,De,n,kx,ky,R,RP,RPP);

for istep= 1:NSTEP
    if (istep==1)
        %Advance saturations
        c0= c;
        ce=c+cs;
        %ux0=ux;
        %uy0=uy;
        for istage= 1:3
            
            f= f_c(c,cs,csx,ux,uy,Pe,kx,ky,kx2,ky2);
            % c0=c;
            
            if istage==1
                c= c + dt*( (8/15)*f );
                [mu,wi,b]=prop(ce,R,RP,RPP,w1,b1);
                [pxt,pyt]=Press_p(mu,wi,b,ce,n,De,uxt,uyt,px,py,ux,uy);
                [g]= g_w(w,c0,ce,csx,ux,uy,R,RP,RPP,b,wi,mu,pxt,pyt,kx,ky,uxt,uyt,De,n);
                
                px = px + dt*( (8/15)*pxt );
                py = py + dt*( (8/15)*pyt );
                w = w + dt*( (8/15)*g );
                
                pxt1=pxt;
                pyt1=pyt;
                f1= f;
                g1=g;
            elseif istage==2
                c= c + dt*( (5/12)*f + (-17/60)*f1 );
                [mu,wi,b]=prop(ce,R,RP,RPP,w1,b1);
                [pxt,pyt]=Press_p(mu,wi,b,ce,n,De,uxt,uyt,px,py,ux,uy);
                [g]= g_w(w,c0,ce,csx,ux,uy,R,RP,RPP,b,wi,mu,pxt,pyt,kx,ky,uxt,uyt,De,n);
                
                px= px + dt*( (5/12)*pxt + (-17/60)*pxt1 );
                py= py + dt*( (5/12)*pyt + (-17/60)*pyt1 );
                w = w + dt*( (5/12)*g + (-17/60)*g1);
                
                
                pxt1=pxt;
                pyt1=pyt;
                f1= f;
                g1=g;
            else
                c= c + dt*( (3/4)*f + (-5/12)*f1 );
                [mu,wi,b]=prop(ce,R,RP,RPP,w1,b1);
                [pxt,pyt]=Press_p(mu,wi,b,ce,n,De,uxt,uyt,px,py,ux,uy);
                [g]= g_w(w,c0,ce,csx,ux,uy,R,RP,RPP,b,wi,mu,pxt,pyt,kx,ky,uxt,uyt,De,n);
                
                px= px + dt*( (3/4)*pxt + (-5/12)*pxt1 );
                py= py + dt*( (3/4)*pyt + (-5/12)*pyt1 );
                w = w + dt*( (3/4)*g + (-5/12)*g1 );
                
            end
            
        end
        t= t+dt;
       % J1 = J(ux,uy,c0,csx,kx,ky);
        [temp,cs,csx]= ini_c(xx,yy,t,Pe);
        W=w;
        Px=px;Py=py;
    end
    
    %Advance saturations
    c0= c;
    ce=cs+c;
    
    %Update new variables
    [mu,wi,b]=prop(ce,R,RP,RPP,w1,b1);  %compute new properties
    [Psi,ux,uy,uxt,uyt]= streamf(W,ux,uy,kx,ky,k2i,dt);     %n             %compute new ux and uy,uxt,uyt
    v2 = (1+ux).*(1+ux)+uy.*uy;
    G = (1+De^2.*v2).^(ce.*((n-1)/2));
    pxt =-(mu.*G.*((1+ux)+(1-b).*wi.*G.*uxt)+Px)./(wi.*G);
    pyt =-(mu.*G.*(uy +(1-b).*wi.*G.*uyt)+Py)./(wi.*G);                                      

    %Predictor step (t+dt)
    [Wp,D1] = WP(W,D1,c,cs,csx,mu,wi,b,pxt,pyt,ux,uy,uxt,uyt,De,n,dt,kx,ky,R,RP,RPP);%predict vorticity
    [Ppx,Ppy,F1x,F1y] = P(Px,Py,F1x,F1y,ce,mu,wi,b,ux,uy,uxt,uyt,De,n,dt);  %predict px,py
    
    J2=J1;
    J1 = J(ux,uy,c0,csx,kx,ky);                                             %new value in J1
    cp = AD(J1,J2,c0,dt,kx,ky);                                             %predict c
    
    
    [temp,cs,csx]= ini_c(xx,yy,t,Pe);
    %Predicted variables
    for j=1:3
        [mu,wi,b]=prop((cp+cs),R,RP,RPP,w1,b1);
        [Psi,upx,upy,uxt,uyt]= streamf(Wp,ux,uy,kx,ky,k2i,dt);
        Jp = J(upx,upy,cp,csx,kx,ky);
        [Fpx,Fpy]=PJ1((cp+cs),mu,wi,b,upx,upy,uxt,uyt,Ppx,Ppy,De,n);
        v2 = (1+upx).*(1+upx)+upy.*upy;
        G = (1+De^2.*v2).^((cp+cs).*((n-1)/2));
        ppxt =-(mu.*G.*((1+upx)+(1-b).*wi.*G.*uxt)+Ppx)./(wi.*G);
        ppyt =-(mu.*G.*(upy +(1-b).*wi.*G.*uyt)+Ppy)./(wi.*G);
        
        [Dp]=WJ1(cp,cs,csx,Wp,mu,wi,b,ppxt,ppyt,upx,upy,uxt,uyt,De,n,kx,ky,R,RP,RPP); %update pxt and pyt also with predicted
        
        %Corrector step
        [Pcx,Pcy] = Pc(Fpx,Fpy,F1x,F1y,Px,Py,dt);
        Wc=WC(W,Dp,D1,dt);
        cc=cor(Jp,J1,c0,kx,ky,dt);
        
        cp=cc;
        Ppx=Pcx;Ppy=Pcy;
        Wp=Wc;
        
    end
    c=cc;W=Wc;Px=Pcx;Py=Pcy;
    
    t= t+dt;
   
    %Update velocities
   %[Psi,ux,uy,uxt,uyt]= streamf(w,ux,uy,kx,ky,k2i,dt);    
   % [temp,cs,csx]= ini_c(xx,yy,t,Pe);
    %td=t*Pe;
    %dt= 1/(8*N*max(max(ux)));
    
    
    if (mod(istep,10)==0) || (istep==1)
        iptsetpref('ImshowBorder','tight');
        handle= figure(1);surf(xx,yy,ce,'facecolor','interp','edgecolor','none','facelighting','phong');
        view([0,0,1]);axis tight equal;title(['t=  ' num2str(t)]);colormap(jet);drawnow;
        %fpath = ['C:\Users\ARC\Desktop\MTP\FigNew\',H,'\',H,'_jet'];
          %savefig([fpath 'PN1' num2str(istep) '.fig']);
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
    end
end
% close(myVideo)
%% 
%% 
function [Ppx,Ppy,Fx,Fy] = P(Px,Py,F1x,F1y,ce,mu,wi,b,ux,uy,uxt,uyt,De,n,dt)
v2 = (1+ux).*(1+ux)+uy.*uy;
G = (1+De^2.*v2).^(ce.*((n-1)/2));
Fx = -(mu.*G.*((1+ux)+(1-b).*wi.*G.*uxt)+Px)./(wi.*G);
Fy = -(mu.*G.*(uy +(1-b).*wi.*G.*uyt)+Py)./(wi.*G);
%Jxhat=fft2(Jx);Jyhat=fft2(Jy);%J1xhat=fft2(J1x);J1yhat=fft2(J1y);

Ppx = Px + (dt/2)*(3*Fx - F1x);
Ppy = Py + (dt/2)*(3*Fy - F1y);

%Ppx=real(ifft2(Ppx));Ppy=real(ifft2(Ppy));
end

function [Fx,Fy]=PJ1(ce,mu,wi,b,ux,uy,uxt,uyt,px,py,De,n)
v2 = (1+ux).*(1+ux)+uy.*uy;
G = (1+De^2.*v2).^(ce.*((n-1)/2));
Fx = -(mu.*G.*((1+ux)+(1-b).*wi.*G.*uxt)+px)./(wi.*G);
Fy = -(mu.*G.*(uy +(1-b).*wi.*G.*uyt)+py)./(wi.*G);

%Jxhat=fft2(Jx);Jyhat=fft2(Jy);
end

function [Pcx,Pcy] = Pc(Fpx,Fpy,Fx,Fy,px,py,dt)
%v2 = (1+ux).*(1+ux)+uy.*uy;
%G = (1+De^2.*v2).^(ce.*((n-1)/2));

Pcx = px + (dt/2)*(Fpx + Fx);
Pcy = py + (dt/2)*(Fpy + Fy);

%Pcx = real(ifft2(Pcx));Pcy = real(ifft2(Pcy));
end

function [D1]=WJ1(c,cs,csx,w,mu,wi,b,pxt,pyt,ux,uy,uxt,uyt,De,n,kx,ky,R,RP,RPP)
i= sqrt(-1);
uxhat=fft2(ux);
uyhat=fft2(uy);
uxx=real(ifft2(i*kx.*uxhat));
uxy = real(ifft2(i*ky.*uxhat));
uyy=real(ifft2(i*ky.*uyhat));
uyx = real(ifft2(i*kx.*uyhat));
chat= fft2(c);

cx=real(ifft2(i*kx.*chat));
cy=real(ifft2(i*ky.*chat));

ce=c+cs;
v2 = (1+ux).*(1+ux)+uy.*uy;
G = (1+De^2.*v2).^((ce.*((n-1)/2)));
Lx =((n-1)/2).*(cx + csx).*log(1+De^2.*v2)+((n-1).*ce.*De^2.*((ux+1).*uxx + uy.*uyx))./(1+De^2.*v2) ;
Ly =((n-1)/2).*cy.*log(1+De^2.*v2)+((n-1).*ce.*De^2.*((ux+1).*uxy + uy.*uyy))./(1+De^2.*v2);
K=(wi.*(1-b)).^-1;

D1 = ((mu.*G.*(1-b)).^-1).*((pxt).*(Ly - RP.*cy)...
    + (-pyt).*(Lx - RP.*(cx+csx)))...
    + K.*( -w./G...
    + ((ux+1)./G).*(Ly - R.*cy)...
    - (uy./G).*(Lx - R.*(cx+csx)))...
    + ((1-b).^-1).*(uxt.*(2*(1-b).*Ly + cy.*(-(R+RP).*(1-b)+RPP.*b))...
    - uyt.*(2*(1-b).*Lx + (cx+csx).*(-(R+RP).*(1-b)+RPP.*b)))...
    ;
%J1hat=fft2(J1);
end


function [Wp,D]=WP(W,D1,c,cs,csx,mu,wi,b,pxt,pyt,ux,uy,uxt,uyt,De,n,dt,kx,ky,R,RP,RPP)

D=WJ1(c,cs,csx,W,mu,wi,b,pxt,pyt,ux,uy,uxt,uyt,De,n,kx,ky,R,RP,RPP);
%Jhat=fft2(J);
%W1=J;
Wp = W + (dt/2)*(3*D -D1);
end

function [Wc]=WC(W,Dp,D1,dt)
%v2 = (1+ux).*(1+ux)+uy.*uy;
%G = (1+De^2.*v2).^(ce.*((n-1)/2));
%K=(wi.*(1-b)).^-1;
%Lhat = fft2(K./G);

Wc = W + (dt/2)*(Dp + D1);

%Wc=real(ifft2(Wc));
end

function [Pxt,Pyt] = Press_p(mu,wi,b,c,n,De,uxt,uyt,px,py,ux,uy)
v2  = (1+ux).*(1+ux)+ uy.*uy;
G   = (1+De^2.*v2).^(c.*((n-1)/2));
Pxt =-(mu.*G.*((1+ux)+(1-b).*wi.*G.*uxt)+px)./(wi.*G);
Pyt =-(mu.*G.*(uy +(1-b).*wi.*G.*uyt)+py)./(wi.*G);
end
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


%cs = 0.5*(1-erf((xx)*sqrt(delta)));
cs = 0.5*(1-erf((xx)/sqrt(4*t)));
csx = -(1.0/sqrt((4*t*pi)))*exp(-(xx.*xx)/(4*t));
%csxx = 0.5/(t*sqrt(4*pi*t)).*xx.*exp(-(xx.*xx)/(4*t));
%cst = 0.25/(t*sqrt(pi*t))*xx.*exp(-(xx.*xx)/(4*t));

ce=0*cs;
 for j=1:1:N
    ce(j,:) = 0.01*rand()*exp(-(xx(j,:)/Pe).^2/(0.01*0.01));
 end

c=ce;
%c=cs+ce;
% N12= N/2;
% for i= 1:N12
%   c(:,i)= c(:,N-(i-1));
% end


end

function [mu,wi,b]=prop(c,R,RP,RPP,w1,b1)
mu = exp(R*(1-c));
wi = w1*exp(RP*(1-c));
b = b1*exp(RPP*(1-c));
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

%f= -( duxcx + duycy - (1/Pe)*(dcxx + dcyy ) );
f= -( duxcx + duycy + (ux).*csx - (dcxx + dcyy) );
end

function [g]= g_w(w,c,ce,csx,ux,uy,R,RP,RPP,b,wi,mu,pxt,pyt,kx,ky,uxt,uyt,De,n)

i= sqrt(-1);
uxhat=fft2(ux);
uyhat=fft2(uy);
uxx=real(ifft2(i*kx.*uxhat));
uxy = real(ifft2(i*ky.*uxhat));
uyy=real(ifft2(i*ky.*uyhat));
uyx = real(ifft2(i*kx.*uyhat));
chat= fft2(c);

cx=real(ifft2(i*kx.*chat));
cy=real(ifft2(i*ky.*chat));
v2 = (1+ux).*(1+ux)+uy.*uy;

G = (1+De^2.*v2).^(ce.*((n-1)/2));
Lx =((n-1)/2).*(cx + csx).*log(1+De^2.*v2)+((n-1).*ce.*De^2.*((ux+1).*uxx + uy.*uyx))./(1+De^2.*v2) ;
Ly =((n-1)/2).*cy.*log(1+De^2.*v2)+((n-1).*ce.*De^2.*((ux+1).*uxy + uy.*uyy))./(1+De^2.*v2); 



K=(wi.*(1-b)).^-1;

g = ((mu.*G.*(1-b)).^-1).*((pxt).*(Ly - RP.*cy)...
    + (-pyt).*(Lx - RP.*(cx+csx)))...
    + K.*( -w./G +...
    + ((ux+1)./G).*(Ly - R.*cy)...
    - (uy./G).*(Lx - R.*(cx+csx)))...
    + ((1-b).^-1).*(uxt.*(2*(1-b).*Ly + cy.*(-(R+RP).*(1-b)+RPP.*b))...
    - uyt.*(2*(1-b).*Lx + (cx+csx).*(-(R+RP).*(1-b)+RPP.*b)))...
    ;
% g = K.*( -w +...
%     (wi./mu).*((-pxt).*(Ly - RP.*cy)...
%     + (pyt).*(Lx - RP.*(cx+csx)))...
%     + (ux+1).*(Ly - R.*cy)...
%     - uy.*(Lx - R.*(cx+csx))...
%     + uxt.*(wi.*(1-b).*Ly + wi.*cy.*(-(R+RP).*(1-b)+RPP.*b))...
%     - uyt.*(wi.*(1-b).*Lx + wi.*(cx+csx).*(-(R+RP).*(1-b)+RPP.*b)))...
%     ;
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


