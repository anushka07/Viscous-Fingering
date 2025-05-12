%function misc2D_VE

%% Grid NxN
N   = 256;
Pe= 1500;
h  = 1/(N-1);
A=2;
x= -1/2:h:1/2;x= x(1:N);% Discretizing in x direction
y= -1/(2*A):h:1/(2*A);y=y(1:N/A); % Discretizing in x direction
[xx,yy]= meshgrid(x,y'); % Creating a rectangular grid with aspect ratio A 

%% Model parameters and initial conditions
R = 3.1437;                     %Log Mobility ratio
n = 0.2;                        %Power law index
De = 2.0;                       %Deborah number
w1 = 0.05;                      %Weissenberg number of displacing fluid
w2 = 0.05;                      %Weissenberd number of displaced fluid
b1 = 0.1;                       %Beta of displacing fluid
b2 = 0.1;                       %Beta of displaced fluid
RPP = log(b2/b1);
RP = log(w2/w1);
%KP=1;
t=0.000125;                     %t0 for base state
[c,cs,csx]= ini_c(xx,yy,t,Pe);  %Base state profiles
ce=c+cs;                        %Total base state concentration profile
%[mu,wi,b]=prop(ce,R,RP,RPP,w1,b1);% Base state parameter profiles

%% Initalizing with zero
Psi= 0*c;   %Streamfunction
ux= 0*c;    % x direction velocity
uy= 0*c;    % y direction velocity
uxi=0*c;    % previous value of ux
uyi=0*c;    % previous value of uy
uxt=0*c;    % time derivative of ux
uyt=0*c;    % time derivative of uy
w=0*c;      % vorticity
px=0*c;     % gradient of pressure in x direction
py=0*c;     % gradient of pressure in y direction

%% Wave numbers for Fourier space...
k1= (2*pi)*[0:(N/2-1) -(N/2):(-1)];             
k2=(2*pi*A)*[0:(N/(2*A)-1) -(N/(2*A)):(-1)];
[kx,ky]= meshgrid(k1,k2);
kx2= kx.^2;                                     % kx^2
ky2= ky.^2;                                     % ky^2
k2i= -(kx.^2 + ky.^2);k2i(1,1)= 1;k2i= k2i.^-1; % 1/(kx^2+ky^2)

%% Creating Directories
H='AB_NP3_Pe3k';
path1 = ['C:\Users\ARC\Desktop\MTP\FigNew\',H,'\',H,'_jet'];    % To store colored images
mkdir(path1);
path2 = ['C:\Users\ARC\Desktop\MTP\FigNew\',H,'\',H,'_gray'];   % To store grayscale images   
mkdir(path2);

%% Initialize video
myVideo = VideoWriter('PN3_Pe5000');                %open video file
myVideo.FrameRate = 5;                              %can adjust this, 5 - 10 works well
open(myVideo)
%% Time step & Iterations 
dt=1/(8*N);
NSTEP= 60000;

%% Advancing in Time 

% Computing functions at time t0
J1 = J(ux,uy,c,csx,kx,ky);
[F1x,F1y]=PJ1(ce,mu,wi,b,ux,uy,uxt,uyt,px,py,De,n);
[D1]=WJ1(c,cs,csx,w,mu,wi,b,pxt,pyt,ux,uy,uxt,uyt,De,n,kx,ky,R,RP,RPP);

%Advancing Saturations
for istep= 1:NSTEP
    if (istep==1)        % First Iteration done using explicit RK3 method
        
        c0= c;           % Initial Perturbation stored in c0
        ce=c+cs;         % Initial total concentration stored in ce
        
        [mu,wi,b]=prop(ce,R,RP,RPP,w1,b1); % Initial parameters
        
        for istage= 1:3   
            
            f= f_c(c,cs,csx,ux,uy,Pe,kx,ky,kx2,ky2);    % Updating f
            
            if istage==1
                c= c + dt*( (8/15)*f );                 % Updating c for first stage               
                     
                [pxt,pyt]=Press_p(mu,wi,b,ce,n,De,uxt,uyt,px,py,ux,uy); % Calculating pxt,pyt
                [g]= g_w(w,c0,ce,csx,ux,uy,R,RP,RPP,b,wi,mu,pxt,pyt,kx,ky,uxt,uyt,De,n); % Updating g
                
                px = px + dt*( (8/15)*pxt );    % Stepping forward in time 
                py = py + dt*( (8/15)*pyt );    % for px, py and w for first stage
                w = w + dt*( (8/15)*g );
                
                pxt1=pxt;   % Storing current values
                pyt1=pyt;
                f1= f;
                g1=g;
                
            elseif istage==2
                c= c + dt*( (5/12)*f + (-17/60)*f1 );   % Updating c for second stage 
                
                [pxt,pyt]=Press_p(mu,wi,b,ce,n,De,uxt,uyt,px,py,ux,uy); % Updating pxt,pyt
                [g]= g_w(w,c0,ce,csx,ux,uy,R,RP,RPP,b,wi,mu,pxt,pyt,kx,ky,uxt,uyt,De,n); % Updating g
                
                px= px + dt*( (5/12)*pxt + (-17/60)*pxt1 ); % Stepping forward
                py= py + dt*( (5/12)*pyt + (-17/60)*pyt1 ); % second stage
                w = w + dt*( (5/12)*g + (-17/60)*g1);
                
                
                pxt1=pxt;   % Storing current values
                pyt1=pyt;
                f1= f;
                g1=g;
            else
                c= c + dt*( (3/4)*f + (-5/12)*f1 );     % Updating c for third stage
               
                [pxt,pyt]=Press_p(mu,wi,b,ce,n,De,uxt,uyt,px,py,ux,uy); % Updating pxt,pyt
                [g]= g_w(w,c0,ce,csx,ux,uy,R,RP,RPP,b,wi,mu,pxt,pyt,kx,ky,uxt,uyt,De,n); % Updating g
                
                px= px + dt*( (3/4)*pxt + (-5/12)*pxt1 );   % Stepping forward
                py= py + dt*( (3/4)*pyt + (-5/12)*pyt1 );   % third stage
                w = w + dt*( (3/4)*g + (-5/12)*g1 );
                
            end
            
        end
        
        t= t+dt;                         % First time step completed
        [~,cs,csx]= ini_c(xx,yy,t,Pe);   % Update base state profile with new time t
        
    end
    
% Advance saturations with Adam Bashforth corrector-predictor method
    c0= c;
    ce=cs+c;
    
    %Update new variables
    [mu,wi,b]=prop(ce,R,RP,RPP,w1,b1);                      % Update parameters
    [Psi,ux,uy,uxt,uyt]= streamf(w,ux,uy,kx,ky,k2i,dt);     % Compute new ux,uy,uxt,uyt
    v2 = (1+ux).*(1+ux)+uy.*uy;                             % Calculating Shear-thinning 
    G = (1+De^2.*v2).^((1-ce).*((n-1)/2));                  % expressions
    pxt =-(mu.*G.*((1+ux)+(1-b).*wi.*G.*uxt)+px)./(wi.*G);  % Updating pxt,pyt
    pyt =-(mu.*G.*(uy +(1-b).*wi.*G.*uyt)+py)./(wi.*G);                                      

    %Predictor step (t+dt)
    [wp,D1] = WP(w,D1,c,cs,csx,mu,wi,b,pxt,pyt,ux,uy,uxt,uyt,De,n,dt,kx,ky,R,RP,RPP);% Predict vorticity
    [ppx,ppy,F1x,F1y] = P(px,py,F1x,F1y,ce,mu,wi,b,ux,uy,uxt,uyt,De,n,dt);           % Predict px,py
    
    J2=J1;                                                      % Storing old value
    J1 = J(ux,uy,c,csx,kx,ky,Pe);                               % New value in J1
    cp = AD(J1,J2,c,dt,kx,ky,Pe);                               % Predict c
    
    t= t+dt;                                                    % Update base cs, csx with
    [temp,cs,csx]= ini_c(xx,yy,t,Pe);                           % next time step
   
    % Iterations with predictor-corrector method
    for j=1:3
        %Predicted variables 
        [Psi,upx,upy,uxt,uyt]= streamf(wp,ux,uy,kx,ky,k2i,dt);  % Predict ux,uy with wp
        Jp = J(upx,upy,cp,csx,kx,ky);                           % Predict J with upx,upy,cp
        [Fpx,Fpy]=PJ1(cp,mu,wi,b,upx,upy,uxt,uyt,ppx,ppy,De,n); % Predict Fp 
        
        v2 = (1+upx).*(1+upx)+upy.*upy;                         % Predict Shear-thinning parameters
        G = (1+De^2.*v2).^((1-(cp+cs)).*((n-1)/2));
        
        ppxt =-(mu.*G.*((1+ux)+(1-b).*wi.*G.*uxt)+ppx)./(wi.*G);    % Predict pxt,pyt
        ppyt =-(mu.*G.*(uy +(1-b).*wi.*G.*uyt)+ppy)./(wi.*G);
        
        [Dp]= WJ1(cp,cs,csx,wp,mu,wi,b,ppxt,ppyt,upx,upy,uxt,uyt,De,n,kx,ky,R,RP,RPP); % Predict WJ1 
        
        %Corrector step
        [pcx,pcy] = Pc(Fpx,Fpy,F1x,F1y,px,py,dt);               % Corrected value of px, py
        wc=WC(w,Dp,D1,dt);                                      % Corrected value of w
        cc=cor(Jp,J1,c,kx,ky,dt,Pe);                            % Corrected calue of c
        
        %Update predicted values with corrected for next iteration
        cp=cc;
        ppx=pcx;ppy=pcy;
        wp=wc;
        
    end
    % Update actual values with corrected values for advanced time t
    c=cc;w=wc;px=pcx;py=pcy;
    ce=cs+c;
    
   % Printing and Storing contour plot
    if (mod(istep,50)==0) || (istep==1)
        iptsetpref('ImshowBorder','tight');
        handle= figure(1);surf(xx,yy,ce,'facecolor','interp','edgecolor','none','facelighting','phong');
        view([0,0,1]);axis tight equal;title(['t=  ' num2str(t)]);colormap(jet);drawnow;
        
        %path1 = ['C:\Users\ARC\Desktop\MTP\FigNew\',H,'\',H,'_jet'];        
        frame = getframe(gca,[0 2  430 214]); %get frame
        L = frame2im(frame);        
        whereToStore=fullfile(path1,[H num2str(istep) '.png']);
        imwrite(L, whereToStore);

%       writeVideo(myVideo, frame);
        colormap(flipud(gray));     %flipud used to get black displacing white
        
        %path2 = ['C:\Users\ARC\Desktop\MTP\FigNew\',H,'\',H,'_gray'];
        frame = getframe(gca,[0 2  430 214]); %get frame
        L = frame2im(frame);        
        whereToStore=fullfile(path2,[H num2str(istep) '.png']);
        imwrite(L, whereToStore);
    end
end
% close(myVideo)

%% Functions
%% For AB method

function [ppx,ppy,Fx,Fy] = P(px,py,F1x,F1y,ce,mu,wi,b,ux,uy,uxt,uyt,De,n,dt)
v2 = (1+ux).*(1+ux)+uy.*uy;
G = (1+De^2.*v2).^((1-ce).*((n-1)/2));
Fx = -(mu.*G.*((1+ux)+(1-b).*wi.*G.*uxt)+px)./(wi.*G);
Fy = -(mu.*G.*(uy +(1-b).*wi.*G.*uyt)+py)./(wi.*G);

ppx = px + (dt/2)*(3*Fx - F1x);
ppy = py + (dt/2)*(3*Fy - F1y);
end

function [Fx,Fy]=PJ1(ce,mu,wi,b,ux,uy,uxt,uyt,px,py,De,n)
v2 = (1+ux).*(1+ux)+uy.*uy;
G = (1+De^2.*v2).^((1-ce).*((n-1)/2));

Fx = -(mu.*G.*((1+ux)+(1-b).*wi.*G.*uxt)+px)./(wi.*G);
Fy = -(mu.*G.*(uy +(1-b).*wi.*G.*uyt)+py)./(wi.*G);
end

function [pcx,pcy] = Pc(Fpx,Fpy,Fx,Fy,px,py,dt)

pcx = px + (dt/2)*(Fpx + Fx);
pcy = py + (dt/2)*(Fpy + Fy);

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
G = (1+De^2.*v2).^(((1-ce).*((n-1)/2)));
Lx =(-(n-1)/2).*(cx + csx).*log(1+De^2.*v2)+((n-1).*(1-ce).*De^2.*((ux+1).*uxx + uy.*uyx))./(1+De^2.*v2) ;
Ly =(-(n-1)/2).*cy.*log(1+De^2.*v2)+((n-1).*(1-ce).*De^2.*((ux+1).*uxy + uy.*uyy))./(1+De^2.*v2);
K=(wi.*(1-b)).^-1;

D1 = ((mu.*G.*(1-b)).^-1).*((pxt).*(Ly - RP.*cy)...
    + (-pyt).*(Lx - RP.*(cx+csx)))...
    + K.*( -w./G...
    + ((ux+1)./G).*(Ly - R.*cy)...
    - (uy./G).*(Lx - R.*(cx+csx)))...
    + ((1-b).^-1).*(uxt.*(2*(1-b).*Ly + cy.*(-(R+RP).*(1-b)+RPP.*b))...
    - uyt.*(2*(1-b).*Lx + (cx+csx).*(-(R+RP).*(1-b)+RPP.*b)))...
    ;

end


function [Wp,D]=WP(W,D1,c,cs,csx,mu,wi,b,pxt,pyt,ux,uy,uxt,uyt,De,n,dt,kx,ky,R,RP,RPP)
D=WJ1(c,cs,csx,W,mu,wi,b,pxt,pyt,ux,uy,uxt,uyt,De,n,kx,ky,R,RP,RPP);

Wp = W + (dt/2)*(3*D -D1);
end

function [Wc]=WC(W,Dp,D1,dt)

Wc = W + (dt/2)*(Dp + D1);

end

function [Pxt,Pyt] = Press_p(mu,wi,b,ce,n,De,uxt,uyt,px,py,ux,uy)
v2  = (1+ux).*(1+ux)+ uy.*uy;
G   = (1+De^2.*v2).^((1-ce).*((n-1)/2));
Pxt =-(mu.*G.*((1+ux)+(1-b).*wi.*G.*uxt)+px)./(wi.*G);
Pyt =-(mu.*G.*(uy +(1-b).*wi.*G.*uyt)+py)./(wi.*G);
end

function[cp] =AD(J1,J2,c,dt,kx,ky,Pe)
% Solved in Fourier space
chat=fft2(c);
J1hat=fft2(J1);J2hat=fft2(J2);
cp = (chat + (dt/2)*(-3*J1hat+J2hat.*exp(-(kx.^2+ky.^2)*dt)/Pe)).*exp(-(kx.^2+ky.^2)*dt/Pe);

cp=real(ifft2(cp));

end

function[cc]=cor(Jp,J2,c,kx,ky,dt,Pe)
% Solved in Fourier space
Jphat=fft2(Jp); J2hat=fft2(J2);
chat=fft2(c);

cc = (chat.*exp(-(kx.^2+ky.^2)*dt/Pe) + (dt/2)*(-Jphat-J2hat.*exp(-(kx.^2+ky.^2)*dt/Pe)));

cc=real(ifft2(cc));
end

function[J1]=J(ux,uy,c,csx,kx,ky,Pe)
i= sqrt(-1);

chat= fft2(c);
dcx= real(ifft2(i*kx.*chat));
dcy= real(ifft2(i*ky.*chat));

J1=ux.*(csx + dcx)+uy.*(dcy);

end

function [c,cs,csx]= ini_c(xx,yy,t,Pe)

N= size(xx,1);
%% For perturbation using specified wavelengths
%if N<300;delta= 2000000;else;delta= 30000;end
% pos= 0*yy;
% w= [10 20 30 40 50 60];
% for iw= 1:length(w);
%    pos= pos + 0.004*(rand-0.5)*cos(w(iw)*pi*yy);
% end;
% 
% arg= (-0.71+xx-4/sqrt(delta)+pos)*sqrt(delta);
% c= 0.5*(1-erf(arg));

%% For perturbation with all wavelengths
cs = 0.5*(1-erf((xx)/sqrt(4*t/Pe)));
csx = -(1.0/sqrt((4*t*pi/Pe)))*exp(-(xx.*xx)/(4*t/Pe));

%csxx = 0.5/(t*sqrt(4*pi*t)).*xx.*exp(-(xx.*xx)/(4*t));
%cst = 0.25/(t*sqrt(pi*t))*xx.*exp(-(xx.*xx)/(4*t));

% Perturbation at interface
c=0*cs;
 for j=1:1:N
    c(j,:) = 0.01*rand()*exp(-(xx(j,:)).^2/(0.01*0.01));
 end

end

function [mu,wi,b]=prop(c,R,RP,RPP,w1,b1)
mu = exp(R*(1-c));
wi = w1*exp(RP*(1-c));
b = b1*exp(RPP*(1-c));
end
%% For RK3 method
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

G = (1+De^2.*v2).^((1-ce).*((n-1)/2));
Lx =(-(n-1)/2).*(cx + csx).*log(1+De^2.*v2)+((n-1).*(1-ce).*De^2.*((ux+1).*uxx + uy.*uyx))./(1+De^2.*v2) ;
Ly =(-(n-1)/2).*cy.*log(1+De^2.*v2)+((n-1).*(1-ce).*De^2.*((ux+1).*uxy + uy.*uyy))./(1+De^2.*v2); 



K=(wi.*(1-b)).^-1;

g = ((mu.*G.*(1-b)).^-1).*((pxt).*(Ly - RP.*cy)...
    + (-pyt).*(Lx - RP.*(cx+csx)))...
    + K.*( -w./G +...
    + ((ux+1)./G).*(Ly - R.*cy)...
    - (uy./G).*(Lx - R.*(cx+csx)))...
    + ((1-b).^-1).*(uxt.*(2*(1-b).*Ly + cy.*(-(R+RP).*(1-b)+RPP.*b))...
    - uyt.*(2*(1-b).*Lx + (cx+csx).*(-(R+RP).*(1-b)+RPP.*b)))...
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


