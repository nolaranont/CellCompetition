%implicit nutrient equation
function [r,fe_r,fe_t,P,V,VT,R,radial,hoop,C,Gam,residual,LamM,LamA] = solve_cells_competition2D_inv(param)
%1-R0 %2-gamma_c %3-D %4-eta %5-k
%6-Nr %7-tspan from 0 %8-mu %9-c_H %10-alpha
%11-dt %12-record_every

R0 = param(1);
gamma_c = param(2);
D = param(3);
scale_v = param(10);
scale_v_ext = 0;
dt = param(11);
record_every = param(12);
% p, v iteration
use_nonlinear = 0;
max_vp_ite = 100; % max iterations of v,p
tol = 1e-6;

disp_progress = 0;
Nr=param(6);
nrsize = 4;
Nrext = (nrsize-1)*(Nr-1)+1;
Nrall = Nr+Nrext;
r = [linspace(0,1,Nr)';linspace(1,nrsize,Nrext)'];
%r = linspace(0,1,Nr+Nrext)';
Nt = round((param(7))/dt);
dx = r(2)-r(1);
dx2 = r(Nr+2)-r(Nr+1);
RT_index = round((length(r)-1)/2);
Pbar = param(16);
bvalue =param(17);
fsolve_ops = optimoptions(@fsolve,'Display','off','algorithm','trust-region');
tol2 = max([tol,dx,dx2]);
%% nonlinear solver
function F = fvp(x)
    % x = [v; p]
    F = ones(size(x));

    % v tumor
    gamT= get_Gamma_T(fer,fet,x(Nr+1:2*Nr));
    gamExt= get_Gamma_Ext(fer,fet,x(Nrext+2*Nr+1:2*Nrall));

    F(1) = x(1); % v(1) = 0; boundary condition for v
    F(2) = -0.5*r(2).*x(1)+dx.*x(2)+0.5*r(2).*x(3)-dx*R(n-1).*r(2).*gamT(2);
    F(3:Nr) = 0.5*r(3:Nr).*x(1:Nr-2)-2*r(3:Nr).*x(2:Nr-1)+(1.5*r(3:Nr)+dx).*x(3:Nr)...
        -dx*R(n-1).*r(3:Nr).*gamT(3:Nr);
    % p tumor
    F(Nr+1:2*Nr-2) = -1.5*x(Nr+1:2*Nr-2)+2*x(Nr+2:2*Nr-1)-0.5*x(Nr+3:2*Nr) ...
        + dx*scale_v*R(n-1)*x(1:Nr-2) - dx*fr(1:Nr-2);
    F(2*Nr-1) = -0.5*x(2*Nr-2)+0.5*x(2*Nr) ...
        + dx*scale_v*R(n-1)*x(Nr-1) - dx*fr(Nr-1);
    F(2*Nr) = x(2*Nr)-x(2*Nr+Nrext+1)+sr(Nr+1)-sr(Nr); %stress jump condition
    %F(2*Nr) = x(2*Nr) - pEndT; %p = sr-Fext; boundary condition for p
    
    %v exterior
    F(2*Nr+1) = x(2*Nr+1)-x(Nr); %v_{Ext}(0) = v_T(R) at the interface
    F(2*Nr+2) = -0.5.*x(2*Nr+1)+dx2.*nrsize./(r(Nr+2).*(nrsize-r(Nr+2))).*x(2*Nr+2)+0.5.*x(2*Nr+3)-dx2*R(n-1).*nrsize.*(nrsize-1).*gamExt(2)./(nrsize-r(Nr+2)).^2;
    F(2*Nr+3:2*Nr+Nrext) = 0.5*x(2*Nr+1:2*Nr+Nrext-2)-2*x(2*Nr+2:2*Nr+Nrext-1)+(1.5+dx2*nrsize./(r(Nr+3:Nrall).*(nrsize-r(Nr+3:Nrall)))).*x(2*Nr+3:2*Nr+Nrext)...
        -dx2*R(n-1).*nrsize.*(nrsize-1).*gamExt(3:Nrext)./(nrsize-r(Nr+3:Nrall)).^2;
    F(2*Nr+Nrext) = 0; %v=0 at r infinity
    %p exterior
    %F(3*Nr+1) = x(3*Nr+1)-x(2*Nr)-sr(Nr+1)+sr(Nr);
    F(2*Nr+Nrext+1:2*Nrall-2) = -1.5*x(2*Nr+Nrext+1:2*Nrall-2)+2*x(2*Nr+Nrext+2:2*Nrall-1)-0.5*x(2*Nr+Nrext+3:2*Nrall) ...
        + dx2*scale_v_ext*R(n-1)*x(2*Nr+1:2*Nr+Nrext-2) - dx2*fr(Nr+1:Nrall-2);
    F(2*Nrall-1) = -0.5*x(2*Nrall-2)+0.5*x(2*Nrall) ...
        + dx2*scale_v_ext*R(n-1)*x(2*Nrall-1) - dx2*fr(Nrall-1);
    F(2*Nrall) = x(2*Nrall)-pEndExt;
end

%% Growth Function
function [Gamma,LamM,LamA] = get_Gamma_T(fer,fet,p)
     W = 1/2*param(8)*(fer(1:Nr).^2+fet(1:Nr).^2-2);
     Ptil = -param(8)*(fer(1:Nr).^2+fet(1:Nr).^2)/2+p+Pbar;
     Sig_inv = 1/2*param(8)*(fer(1:Nr).^2+fet(1:Nr).^2)-p;
     LamM = 0.5*param(5)*ones(Nr,1);
     %LamM = param(5)*(1+param(15).*abs(Sig_inv).*(Sig_inv>-tol2)./(1+abs(Sig_inv).*(Sig_inv>-tol2)));
     %LamA = param(4).*param(19)*abs(Sig_inv).*(Sig_inv<tol2)./(1+0.1*abs(Sig_inv).*(Sig_inv<tol2));
     %LamA = param(4).*param(19)*(Sig_inv).^2.*(Sig_inv<tol2)./(1+0.1*(Sig_inv).^2.*(Sig_inv<tol2));
     LamA = (Ptil+W);
     %Gamma = param(4)*(0.5*param(5)*c(1:Nr).^2-p-W);
     
     %Gamma = param(4)*(0.5*param(5).*lambdaT.*c(1:Nr).^2-(Ptil+W)).*c(1:Nr);
     %Gamma = param(4)*(0.5*param(5).*c(1:Nr).^2-(Ptil+W)).*c(1:Nr);
     %Gamma = param(4)*(0.5*param(5)-(Ptil+W));
     Gamma = LamM-LamA;
     %Gamma = (LamM-(Ptil+W));
end

function [Gamma,LamM,LamA] = get_Gamma_Ext(fer,fet,pext)
    
    %Gamma = zeros(size(pext));
    Wext = 1/2*param(9)*(fer(Nr+1:Nrall).^2+fet(Nr+1:Nrall).^2-2);
    Ptil = -param(9)*(fer(Nr+1:Nrall).^2+fet(Nr+1:Nrall).^2)/2+pext+Pbar;
 %   epsilon = 0.01;
    %alpha = param(15);
     Sig_inv = 1/2*param(9)*(fer(Nr+1:Nrall).^2+fet(Nr+1:Nrall).^2)-pext;
     LamM = param(14)*(param(20)*abs(Sig_inv).*(Sig_inv>-tol2)./(1+abs(Sig_inv).*(Sig_inv>-tol2)));
     %LamA = param(13).*param(21)*abs(Sig_inv).*(Sig_inv<tol2)./(1+abs(Sig_inv).*(Sig_inv<tol2));
     LamA = param(13).*param(21)*(Sig_inv).^2.*(Sig_inv<tol2)./(1+0.1*(Sig_inv).^2.*(Sig_inv<tol2));

    %lambdaE = 1-0.5*(tanh((r(Nr+1:Nrall)-1)/epsilon)+1);
    %lambdaE = exp(-alpha*(r(Nr+1:Nrall)-1));
    %Gamma = param(13)*(0.5*param(14)-(pext+Wext));
    %Gamma = param(13)*(0.5*param(14).*lambdaE-(Ptil+Wext));
    %Gamma(1)
    Gamma = LamM-LamA;
    %Gamma = (LamM-(Ptil+Wext));
end


%% initialize or restart


R = zeros(Nt+1,1); R(1) = R0;
fer = ones(length(r),1);
fet = ones(length(r),1);
c = ones(length(r),1);
%c = (R(1)*sinh(r))./(r*sinh(R(1)));
fe_r(:,1) = fer;
fe_t(:,1) = fet;
C(:,1) = c;
Gam(:,1) = [get_Gamma_T(ones(length(r),1),ones(length(r),1),0);get_Gamma_Ext(ones(length(r),1),ones(length(r),1),0)];
p= zeros(Nr,1);
pext = zeros(Nrext,1);
nStart = 2;
residual = zeros(Nt,1); residual(1) = NaN;

%% finite difference loop
time = cputime;
for n=nStart:Nt+1
%currTime = dt*(n-1);
%r_ext = R(n-1)*(nrsize*R(n-1)-1).*r(Nr+1:Nrall)./(nrsize*R(n-1)-r(Nr+1:Nrall));
% if disp_progress && mod(currTime-tspan(1),disp_progress)==0
%   disp(['T=' num2str(currTime)]);
% end

%% calculate stress from lagged p. Note: fet is current
%fe_tr = gradient(fet, dx);

%fe_tr = [(-3*fet(1)+4*fet(2)-fet(3))/2/dx;(fet(3)-fet(1))/2/dx;(fet(1:Nr-2)-4*fet(2:Nr-1)+3*fet(3:Nr))/2/dx;...
%    (-3*fet(Nr)+4*fet(Nr+1)-fet(Nr+2))/2/dx2;(fet(Nr+3)-fet(Nr+1))/2/dx2;(fet(Nr+1:Nrall-2)-4*fet(Nr+2:Nrall-1)+3*fet(Nr+3:Nrall))/2/dx2];
%fe_tr = [(-3*fet(1)+4*fet(2)-fet(3))/2/dx;(fet(3)-fet(1))/2/dx;(fet(1:Nr-2)-4*fet(2:Nr-1)+3*fet(3:Nr))/2/dx;...
%    (-3*fet(Nr+1:Nrall-2)+4*fet(Nr+2:Nrall-1)-fet(Nr+3:Nrall))/2/dx2;(fet(Nrall)-fet(Nrall-2))/2/dx2;(fet(Nrall-2)-4*fet(Nrall-1)+3*fet(Nrall))/2/dx2];

fe_tr = [(-3*fet(1)+4*fet(2)-fet(3))/2/dx;(fet(3:Nr)-fet(1:Nr-2))/2/dx;(fet(Nr-2)-4*fet(Nr-1)+3*fet(Nr))/2/dx;...
   (-3*fet(Nr+1)+4*fet(Nr+2)-fet(Nr+3))/2/dx2;(fet(Nr+3:Nrall)-fet(Nr+1:Nrall-2))/2/dx2;(fet(Nrall-2)-4*fet(Nrall-1)+3*fet(Nrall))/2/dx2];


sr = [param(8)*(fer(1:Nr).^2);param(9)*(fer(Nr+1:Nrall).^2)];
st = [param(8)*(fet(1:Nr).^2);param(9)*(fet(Nr+1:Nrall).^2)];
%sr = [param(8)*0.5*(fer(1:Nr).^2-fet(1:Nr).^2);param(9)*0.5*(fer(Nr+1:2*Nr).^2-fet(Nr+1:2*Nr).^2)];
%st = [param(8)*0.5*(fet(1:Nr).^2-fer(1:Nr).^2);param(9)*0.5*(fet(Nr+1:2*Nr).^(2)-fer(Nr+1:2*Nr).^2)];

[gamT,lamMT,lamAT] = get_Gamma_T(fer,fet,p);
[gamExt,lamMExt,lamAExt] = get_Gamma_Ext(fer,fet,pext);
%% update v and p together
gel_stress = param(9).*(log(R0./R(n-1))+0.5.*(1-(R0./R(n-1)).^2));
sigma_r3 = [param(8).*(-2.*(fet(1:Nr)).^(-3)).*fe_tr(1:Nr);param(9).*(-2.*(fet(Nr+1:Nrall)).^(-3)).*fe_tr(Nr+1:Nrall)];
%sigma_r3 = [param(8)*(-(fet(1:Nr)).^(-3)-fet(1:Nr)).*fe_tr(1:Nr);param(9)*(-(fet(Nr+1:2*Nr)).^(-3)-fet(Nr+1:2*Nr)).*fe_tr(Nr+1:2*Nr)];

%pEndT = sr(Nr)- gel_stress;
pEndExt = sr(Nrall);

fr = 1.*(sigma_r3 + 1./r.*(sr - st));
fr(1) = sigma_r3(1);
fr(Nr+1:Nrall) = sigma_r3(Nr+1:Nrall)+nrsize./(r(Nr+1:Nrall).*(nrsize-r(Nr+1:Nrall))).*(sr(Nr+1:Nrall)-st(Nr+1:Nrall));
N = 2*Nrall; % 2x large matrix


diagV = [1; dx.*1; (1.5*r(3:Nr)+dx).*1];
subDiagV = [-0.5*r(2)*1;-2*r(3:Nr)*1];
subsubDiagV = 0.5*r(3:Nr)*1;
supDiagV = [0;0.5*r(2)*1;zeros(Nr-2,1)];

diagP = [-1.5*ones(Nr-2,1);0;1];
subDiagP = [zeros(Nr-2,1);-0.5;0];
supDiagP = [2*ones(Nr-2,1);0.5];
supsupDiagP = -0.5*ones(Nr-2,1);

VPcoupled = [ones(Nr-1,1)*scale_v*R(n-1)*dx; 0]; % lower left submatrix. the p_r+R*v bit

% diagVExt = [1; dx2.*1; (1.5*r(Nr+3:Nrall)+dx2).*1];
% subDiagVExt = [0;-0.5*r(Nr+2)*1;-2*r(Nr+3:Nrall)*1];
% subsubDiagVExt = 0.5*r(Nr+3:Nrall)*1;
% supDiagVExt = [0;0;0.5*r(Nr+2)*1;zeros(Nrext-2,1)];
diagVExt = [1; dx2.*nrsize./(r(Nr+2).*(nrsize-r(Nr+2))); (1.5+dx2*nrsize./(r(Nr+3:Nrall-1).*(nrsize-r(Nr+3:Nrall-1)))).*1;1];
subDiagVExt = [0;-0.5;-2*ones(length(Nr+3:Nrall-1),1)*1;0];
subsubDiagVExt = [0.5*ones(length(Nr+3:Nrall-1),1)*1;0];
supDiagVExt = [0;0;0.5*1;zeros(Nrext-2,1)];

diagPExt = [-1.5*ones(Nrext-2,1);0;1];
subDiagPExt = [zeros(Nrext-2,1);-0.5;0];
supDiagPExt = [2*ones(Nrext-2,1);0.5];
supsupDiagPExt = [-0.5*ones(Nrext-2,1)];

VPcoupledExt = [ones(Nrext-1,1)*scale_v_ext*R(n-1)*dx2; 0]; % lower left submatrix. the p_r+R*v bit

diag = [diagV; diagP; diagVExt; diagPExt];
supDiag = [supDiagV; supDiagP;supDiagVExt;supDiagPExt];
subDiag = [subDiagV; subDiagP;subDiagVExt;subDiagPExt];
A = sparse(1:N,1:N,diag);
A = A + sparse(1:N-1,2:N,supDiag,N,N);
A = A + sparse(2:N,1:N-1,subDiag,N,N);
A = A + sparse(Nr+1:2*Nr,1:Nr,VPcoupled,N,N);
A = A + sparse(Nr+1:2*Nr-2,Nr+3:2*Nr,supsupDiagP,N,N);
A = A + sparse(3:Nr,1:Nr-2,subsubDiagV,N,N);
A = A + sparse(2*Nr+Nrext+1:2*Nrall,2*Nr+1:2*Nr+Nrext,VPcoupledExt,N,N);
A = A + sparse(2*Nr+Nrext+1:N-2,2*Nr+Nrext+3:N,supsupDiagPExt,N,N);
A = A + sparse(2*Nr+3:2*Nr+Nrext,2*Nr+1:2*Nr+Nrext-2,subsubDiagVExt,N,N);
A = A + sparse(2*Nr+1,Nr,-1,N,N); %jump cond for v
A = A + sparse(2*Nr,2*Nr+Nrext+1,-1,N,N);%jump cond for p

if ~use_nonlinear
    for i_vp=1:max_vp_ite % iterate a few times to eliminate oscillation
        % RHS of v and p

        rhsV = [0; (dx*R(n-1)).*(r(2:Nr)*1).*gamT(2:Nr)];
        rhsP = [dx.*fr(1:Nr-1); sr(Nr)-sr(Nr+1)];
        %rhsVExt = [0; (dx2*R(n-1)).*(r(Nr+2:Nrall)*1).*gamExt(2:Nrext)];
        rhsVExt = [0; dx2*R(n-1).*nrsize.*(nrsize-1).*gamExt(2:Nrext-1)./(nrsize-r(Nr+2:Nrall-1)).^2;0];
        rhsPExt = [dx2.*fr(Nr+1:Nrall-1); pEndExt];

        rhs = [rhsV; rhsP;rhsVExt;rhsPExt];
        
        
        % solution is [v; p;v_{Ext};p_{Ext}]
        sol = A\rhs;
        v = sol(1:Nr);
        p = sol(Nr+1:2*Nr);
        vext = sol(2*Nr+1:2*Nr+Nrext);
        pext = sol(2*Nr+Nrext+1:N);
        [gamT,lamMT,lamAT] = get_Gamma_T(fer,fet,p);
        [gamExt,lamMExt,lamAExt] = get_Gamma_Ext(fer,fet,pext);
        
        res = norm(fvp([v; p; vext; pext]));

        if res<tol
            break
        else
            if i_vp==max_vp_ite, use_nonlinear = 1; end
        end
    end
else
    x0 = [v; p; vext; pext];
    sol = fsolve(@fvp, x0, fsolve_ops);
    v = sol(1:Nr);
    p = sol((Nr+1):(2*Nr));
    vext = sol(2*Nr+1:2*Nr+Nrext);
    pext = sol(2*Nr+Nrext+1:N);
    [gamT,lamMT,lamAT] = get_Gamma_T(fer,fet,p);
    [gamExt,lamMExt,lamAExt] = get_Gamma_Ext(fer,fet,pext);
end

residual(n) = norm(fvp([v; p;vext;pext]));

%% update R

R(n) = R(n-1) + dt*v(Nr);
if R(n) < tol
    R(n)=0;
    break;
end
%% solve fer
drdt = v(Nr);
vt = (v-r(1:Nr)*drdt)/R(n);
vtext = (vext.*(nrsize-r(Nr+1:Nrall)).^2./R(n)./nrsize./(nrsize-1)...
         -(nrsize-r(Nr+1:Nrall)).*r(Nr+1:Nrall)*drdt./R(n)./nrsize);

v_r = zeros(Nr,1);
v_r_ext = zeros(Nr,1);
v_r(1) = gamT(1)/2;
v_r(2:Nr) = gamT(2:Nr)-v(2:Nr).*(1./(r(2:Nr).*R(n-1)));
v_r_ext(2:Nrext) = gamExt(2:Nrext)-vext(2:Nrext).*(1./(r(Nr+2:Nrall).*R(n-1)));
v_r_ext(1) = v_r(Nr);
%setup for 2nd order upwind solver
fetnew = zeros(Nrall,1);
fetnew(1)  = fet(1);
upwc = zeros(Nr,1);
upwc(1) = 0;
upwc(end) = 0;

Fet_r_p = zeros(Nr,1);
Fet_r_p(1:Nr-2) = (-3.*fet(1:Nr-2)+4.*fet(2:Nr-1)-fet(3:Nr))./2./dx;
Fet_r_p(Nr-1) = (fet(Nr)-fet(Nr-2))./2./dx;
Fet_r_m = zeros(Nr,1);
Fet_r_m(3:Nr) = (3.*fet(3:Nr)-4.*fet(2:Nr-1)+fet(1:Nr-2))./2./dx;
Fet_r_m(2) = (fet(3)-fet(1))./2./dx;

Fet_r_pExt = zeros(Nrext,1);
Fet_r_pExt(1:Nrext-2) = (-3.*fet(Nr+1:Nrall-2)+4.*fet(Nr+2:Nrall-1)-fet(Nr+3:Nrall))./2./dx2;
Fet_r_pExt(Nrext-1) = (fet(Nrall)-fet(Nrall-2))./2./dx2;
Fet_r_mExt = zeros(Nrext,1);
Fet_r_mExt(3:Nrext) = (3.*fet(Nr+3:Nrall)-4.*fet(Nr+2:Nrall-1)+fet(Nr+1:Nrall-2))./2./dx2;
Fet_r_mExt(2) = (fet(Nr+3)-fet(Nr+1))./2./dx2;

%b = ones(Nrext,1).*bvalue;
b = (linspace(0,1,Nrext).^param(18).*bvalue).'; %constant when param(18)=0;
%b = [zeros(Nrext-200,1);ones(200,1)].*bvalue;
D2 = gamT/2;%+1/3*b.*(fet.^2-fet.^(-4));
D2ext = gamExt/2 + b/2.*(fet(Nr+1:Nrall).^2-fet(Nr+1:Nrall).^(-2));

fetnew(1) = fet(1);
fetnew(2:Nr) = fet(2:Nr)-dt*(max(vt(2:Nr),0).*Fet_r_m(2:Nr)+min(vt(2:Nr),0).*Fet_r_p(2:Nr))...
    +dt*(v(2:Nr)./r(2:Nr)./R(n-1)-D2(2:Nr)).*fet(2:Nr);
fetnew(Nr+1:Nrall) = fet(Nr+1:Nrall)-dt*(max(vtext(1:Nrext),0).*Fet_r_mExt(1:Nrext)+min(vtext(1:Nrext),0).*Fet_r_pExt(1:Nrext))...
    +dt*(vext(1:Nrext).*(nrsize-r(Nr+1:Nrall))./r(Nr+1:Nrall)./R(n-1)./(nrsize-1)-D2ext(1:Nrext)).*fet(Nr+1:Nrall);
%fetnew(Nr+1) = fetnew(Nr);
%fetnew(Nr) = fetnew(Nr+1);
fernew = fetnew.^(-1);

% c_r_p = zeros(Nr,1);
% c_r_p(1:Nr-2) = (-3.*c(1:Nr-2)+4.*c(2:Nr-1)-c(3:Nr))./2./dx;
% c_r_p(Nr-1) = (c(Nr)-c(Nr-1))./dx;
% c_r_m = zeros(Nr,1);
% c_r_m(3:Nr) = (3.*c(3:Nr)-4.*c(2:Nr-1)+c(1:Nr-2))./2./dx;
% c_r_m(2) = (c(2)-c(1))./dx;
% upwc(2:Nr-1) = (max(vt(2:Nr-1),0).*c_r_m(2:Nr-1)+min(vt(2:Nr-1),0).*c_r_p(2:Nr-1));
% 
%iterate the nutrient solver n times at the first step
if n == 2
    max_inter = 10;
else
    max_inter = 1;
end

% for i = 1:max_inter
%     diagC    = [-1.5;(1+dt*(gamma_c+2*D2(2:Nr-1))+dt*2*D/(dx)^2/R(n-1)^2);1];
%     supDiagC = [2;-dt*D/((dx)*R(n-1)^2).*(1/dx+1./r(2:Nr-1))];
%     subDiagC = [-dt*D/(dx*R(n-1)^2).*(1/dx-1./r(2:Nr-1));0];
%     Amatrix_c = sparse(1:Nr,1:Nr,diagC,Nr,Nr)...
%         + sparse(1:Nr-1,2:Nr,supDiagC,Nr,Nr)...
%         + sparse(2:Nr,1:Nr-1,subDiagC,Nr,Nr)...
%         +sparse(1,3,-0.5,Nr,Nr);
%     Amatrix_c = full(Amatrix_c);
%     Bmatrix_c = [0;c(2:Nr-1);1]-dt*upwc;
%     cnew = Amatrix_c\Bmatrix_c;
%     %c = [cnew;ones(Nrext,1)];
%     
% end

 c = [ones(Nr,1);ones(Nrext,1)];
%Update fer fet
fet = fetnew;
fer = fernew;

if mod(n-1,record_every)==0
    nn = (n-1)/record_every+1;
    fe_r(:,nn) = fer;
    fe_t(:,nn) = fet;
    P(:,nn) = [p;pext];
    V(:,nn) = [v;vext];
    sr2(:,nn) = sr;
    st2(:,nn) = st;
    VT(:,nn) = [vt;vtext];
    C(:,nn) = c;
    Gam(:,nn) = [gamT;gamExt];
    LamM(:,nn) = [lamMT;lamMExt];
    LamA(:,nn) = [lamAT;lamAExt];

end
currenttime = cputime-time;
% if currenttime > 60*15/dt/dx/200/10^3
%     break;
% end

end
%radial and hoop stress
radial = sr2 - P;
hoop = st2 - P;
end

