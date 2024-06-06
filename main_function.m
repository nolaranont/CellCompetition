function [] = main_function(cH,K,KE,alpha,Pbar,etae,beta,power,R0,sen1,sen2,sen3)
%==List of parameters==%  
%1-R0
%2-gamma_c
%3-D
%4-eta lamAbar
%5-k
%6-Nr
%7-from 0 to tspan 
%8-mu
%9-c_H
%10-alpha
%11-dt
%12-record_every
%13-eta_{ext} lamAbar
%14-k_{ext}
%15-alpha or delta
%16-Pbar
%17-Beta_ext
%======================%
param = zeros(12,1);
param(1) = R0; param(2)=1; param(3) = 1.^2; param(4)=1; param(5)=K;
param(6) = 101; param(7)=100; param(8) = 1; param(9)=cH; param(10)= 1;
param(11) = (1e-3);param(12) = (1e3);param(13)=etae;param(14) = KE;
param(15) = alpha; param(16) = Pbar; param(17)=beta; param(18)= power;
param(19) = sen1; param(20) = sen2; param(21) = sen3;
resume = 0;

%save('4RCH1case1v2.mat','r','fe_r','fe_t','P','V','R','radial','hoop','C','Gam','param')
%save('4RCH1case2etae0.1v2.mat','r','fe_r','fe_t','P','V','R','radial','hoop','C','Gam','param')
%filename = sprintf('4RMCH%.2fETAE%.2fKE%.2fALPHA%d',param(9),param(13),param(14),param(15));
%filename = sprintf('C1K%.2fPbar%.2fCH%.2fETAE%.2fKE%.2fALPHA%dBeta%.2fPower%dR0%.2f',param(5),param(16),param(9),param(13),param(14),param(15),param(17),param(18),param(1));
%filename = sprintf('4RMPbar%.2fCH%.2fETAE%.2fKE%.2fALPHA%dBeta%.2fTail200',param(16),param(9),param(13),param(14),param(15),param(17));
%filename = sprintf('PWhighC1K%.2fPbar%.2fCH%.2fETAE%.2fKE%.2fALPHA%.2fBeta%.2fPower%dR0%.2f',param(5),param(16),param(9),param(13),param(14),param(15),param(17),param(18),param(1));
%filename = sprintf('200LAM%dC1K%.2fPbar%.2fCH%.2fETAE%.2fKE%.2fALPHA%.2fBeta%.2fPower%dR0%.2f',param(4),param(5),param(16),param(9),param(13),param(14),param(15),param(17),param(18),param(1));
%filename = sprintf('PWalpha1%2dLAM%dC1K%.2fPbar%.2fCH%.2fETAE%.2fKE%.2fALPHA%.2fBeta%.2fPower%dR0%.2fsen1%.2fsen2%.2fsen3%.2f',param(6)-1,param(4),param(5),param(16),param(9),param(13),param(14),param(15),param(17),param(18),param(1),param(19),param(20),param(21));
%filename = sprintf('2D%2d-4LAM%dC1K%.2fPbar%.2fCH%.2fETAE%.2fKE%.2fALPHA%.2fBeta%.2fPower%dR0%.2fsen1%.2fsen2%.2fsen3%.2f',param(6)-1,param(4),param(5),param(16),param(9),param(13),param(14),param(15),param(17),param(18),param(1),param(19),param(20),param(21));
filename = sprintf('G%ddt10M%2dLAM%dC1K%.2fPbar%.2fCH%.2fETAE%.2fKE%.2fALPHA%.2fBeta%.2fPower%dR0%.2fsen1%.2fsen2%.2fsen3%.2f',param(11),param(6)-1,param(4),param(5),param(16),param(9),param(13),param(14),param(15),param(17),param(18),param(1),param(19),param(20),param(21));

if resume
data = load(filename);
[r,fe_r,fe_t,P,V,VT,R,radial,hoop,C,Gam,residual,LamM,LamA] = solve_cells_competition2D_FixExtB_resume(data);    
else
%if ~isfile([filename,'.mat'])
%[r,fe_r,fe_t,P,V,VT,R,radial,hoop,C,Gam,residual,LamM,LamA] = solve_cells_competition2D(param);
[r,fe_r,fe_t,P,V,VT,R,radial,hoop,C,Gam,residual,LamM,LamA] = solve_cells_competition2D_FixExtB(param);
%[r,fe_r,fe_t,P,V,VT,R,radial,hoop,C,Gam,residual] = solve_cells_competition2D_v2(param);
%[r,fe_r,fe_t,P,V,VT,R,radial,hoop,C,Gam,residual,LamM,LamA] = solve_cells_competition2D_inv(param);
%[r,fe_r,fe_t,P,V,VT,R,radial,hoop,C,Gam,residual] = solve_cells_2D(param);
%[r,fe_r,fe_t,P,V,VT,R,radial,hoop,C,Gam,residual] = solve_cells_3D(param);
save([filename,'.mat'],'-v7.3','r','fe_r','fe_t','P','V','R','radial','hoop','C','Gam','param','LamM','LamA')
%save([filename,'.mat'],'-v7.3','r','fe_r','fe_t','P','V','R','radial','hoop','C','Gam','param')

%end
end
end

%main_function(1,1,1,0,0,1,0,0,1,0,0,0)
%main_function(2,1,1,0,0,1,0,0,1,0,0,0)
%main_function(0.5,1,1,0,0,1,0,0,1,0,0,0)
%main_function(0.1,1,1,0,0,1,0,0,1,0,0,0)
%main_function(1,1,2,0,0,1,0,0,1,0,0,0)
%main_function(1,1,0.5,0,0,1,0,0,1,0,0,0)
%main_function(1,1,4,0,0,1,0,0,1,0,0,0)
%main_function(1,1,0.1,0,0,1,0,0,1,0,0,0)
%main_function(1,3,0,0,0,0,0,0,1,0,0,0)