close all;
clear all;
clc
load('c.mat')
load('D.mat')
load('tau.mat')

c(1)=1;
D(1)=5000;
tau(1)=20*86400;
length_t = 30;
line_w = 1.5;
figure,
dt = (0:0.2:length_t )*86400;
dx = ones(1,length(0:0.2:length_t));
subplot(2,2,1)

 [ R1 ] = corrModel_AD_eq( [c(1),D(1),tau(1)] , dt , dx*0 );
 plot(dt/86400,R1,'Linewidth',line_w)
  hold on
 [ R2 ] = corrModel_AD_eq( [c(1),D(1),tau(1)] ,dt , dx*50000 );
 plot(dt/86400,R2,'Linewidth',line_w)
  hold on
 [ R3 ] = corrModel_AD_eq( [c(1),D(1),tau(1)] , dt  , dx*200000 );
  plot(dt/86400,R3,'Linewidth',line_w)
  hold on
  [ R4 ] = corrModel_AD_eq( [c(1),D(1),tau(1)] , dt , dx*500000 );
   plot(dt/86400,R4,'Linewidth',line_w)
  hold on
  xlim([0,25])
  ylabel("\rho_{Q}")
  hold on
    subtitle('(a)')
  xlabel("\Deltat (days)")
  set(gca,'FontSize',25)
  legend({"\Deltax = 0 km","\Deltax = 50 km","\Deltax = 200 km","\Deltax = 500 km"})
  
  subplot(2,2,2)

 [ R1 ] = corrModel_AD_eq( [0.1,D(1),tau(1)] , dt , dx*100000 );
 plot(dt/86400,R1,'Linewidth',line_w)
  hold on
 [ R2 ] = corrModel_AD_eq( [c(1),D(1),tau(1)] ,dt , dx*100000 );
 plot(dt/86400,R2,'Linewidth',line_w)
  hold on
 [ R3 ] = corrModel_AD_eq( [2,D(1),tau(1)] , dt  , dx*100000 );
  plot(dt/86400,R3,'Linewidth',line_w)
  hold on
  [ R4 ] = corrModel_AD_eq( [3,D(1),tau(1)] , dt , dx*100000 );
   plot(dt/86400,R4,'Linewidth',line_w)
  hold on
    xlim([0,25])
  ylabel("\rho_{Q}")
  hold on
  xlabel("\Deltat (days)")
  subtitle('(b)')
  set(gca,'FontSize',25)
  legend({"c_0 = 0.1 m/s","c_0 = 1 m/s","c_0 = 2 m/s","c_0 = 3 m/s"})
  
    subplot(2,2,3)

 [ R1 ] = corrModel_AD_eq( [c(1),1000,tau(1)] , dt , dx*100000 );
 plot(dt/86400,R1,'Linewidth',line_w)
  hold on
 [ R2 ] = corrModel_AD_eq( [c(1),D(1),tau(1)] ,dt , dx*100000 );
 plot(dt/86400,R2,'Linewidth',line_w)
  hold on
 [ R3 ] = corrModel_AD_eq( [c(1),50000,tau(1)] , dt  , dx*100000);
  plot(dt/86400,R3,'Linewidth',line_w)
  hold on
  [ R4 ] = corrModel_AD_eq( [c(1),500000,tau(1)] , dt , dx*100000 );
   plot(dt/86400,R4,'Linewidth',line_w)
  hold on
    xlim([0,25])
  ylabel("\rho_{Q}")
  hold on
  xlabel("\Deltat (days)")
  subtitle('(c)')
  set(gca,'FontSize',25)
  legend({"D_0 = 1000 m^2/s","D_0 = 5000 m^2/s","D_0 = 50000 m^2/s","D_0 = 500000 m^2/s"})
    subplot(2,2,4)
  [ R1 ] = corrModel_AD_eq( [c(1),D(1),10*86400] , dt , dx*100000 );
 plot(dt/86400,R1,'Linewidth',line_w)
  hold on
 [ R2 ] = corrModel_AD_eq( [c(1),D(1),20*86400] ,dt , dx*100000 );
 plot(dt/86400,R2,'Linewidth',line_w)
  hold on
 [ R3 ] = corrModel_AD_eq( [c(1),D(1),30*86400] , dt  , dx*100000);
  plot(dt/86400,R3,'Linewidth',line_w)
  hold on
  [ R4 ] = corrModel_AD_eq( [c(1),D(1),40*86400] , dt , dx*100000 );
   plot(dt/86400,R4,'Linewidth',line_w)
  hold on
    xlim([0,25])
  ylabel("\rho_{Q}")
  hold on
  xlabel("\Deltat (days)")
  subtitle('(d)')
  set(gca,'FontSize',25)
  legend({"\tau = 10 days","\tau = 20 days","\tau = 30 days","\tau = 40 days"})