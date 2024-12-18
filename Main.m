function Main(n)
clc
clearvars -except n
close all
% Wrote by zhe zhang EPFL LWE Nov. 29 2024, with Adrien, Romane, Pierre.
% For CC network calculation, only square, its experimental realization is QSHE coupled ring resonator networks
% For Checking MF in CC network at the critical point

% Functionalities in this code: Cavity spectra, finite transport, supercell bandstructure/ pumping
% One unitcell contains 4 nodes, and 8 ports.
% Unitcell is numbered by (ix,iy), They have the rank of cell, Rc=Mx*(ix-1)+iy
% As one unitcell contains four nodes, in one unitcell they are indexed by Index_node_incell=1~4
% As one unitcell contains eight ports, in one unitcell they are indexed by Index_port_incell=1~8
% Node index: Index_node: (Rc-1)*4+i
% Port index: (Rc-1)*8+i
% The key unitary matrix and all the informations are based on the ports.
% But it is also important to have Node and unitcell for the debug and checking, also plotting for eigenstates.


%% Parameter Scan
%Scaling
L_line=[50 75 125 175 250 375 550];%For Mx My, L = 2*Mx
N_L=length(L_line);%50 100 200 250 300 350 400

N_realization_line=[500 200 60 30 12 5 2]*10^4;%4
Num_eigenstate=1;

q_line=linspace(-1,3,81);
% q_line(q_line==0)=[];q_line(q_line==1)=[];% no 0 and 1, if there is delete
N_q=length(q_line);%center at 1/2


d_edge=[2,2,1,0,2,1];
Data_N=[N_q,Num_eigenstate];
%Setting
Theta_clean=45*pi/180;%in angle unit, uniform, 0-90
Eps_center=0;%take the most Num_eigenstate eigenstates whose quasienergy value is nearly at Eps_center.
% Disorder
Disorder_phase_level=2*pi;%in radius range, uniform 0-2pi
Disorder_Theta_level=0;%in radius unit, uniform, 0-pi/2
Graph_off=2;

%Statistic setting
LDOS_Range=[-20 5];
N_bar=200;


%% Results&Data check
% % Probaility analysis for fixed L. Do not need variations, no relation with moment q
% % Dimension N_L,N_bar 2, storing, avg and std
% 
% Probaility_LDOS_all_avg_std_L=zeros(N_bar,2);
% Probaility_LDOS_Bulk_avg_std_L=zeros(N_bar,2);
% Probaility_LDOS_Edge_avg_std_L=zeros(N_bar,2);
% Probaility_LDOS_Corner_avg_std_L=zeros(N_bar,2);
% Probaility_LDOS_Bulk_cylinder_avg_std_L=zeros(N_bar,2);
% Probaility_LDOS_Edge_cylinder_avg_std_L=zeros(N_bar,2);
% 
% %MF exponents Tau_qL (Delta_qL,Gama_q, in post data processing)
% %And Legendre transform MF singularity spectra, Alpha and f(alpha)
% %The output is three dimensional data: L_line*q_line*N_realization
% %Each basic calculation contain a q line.
% %Therefore the loop is only about the L_line*N_realization, cluster is for later one.
% % Dimension, N_L by 6 by N_q, by 2
% % 6： %_all,_Bulk,_Edge,_Corner for each column, each element in this cell data is of N_q*N_realization dimension
% % 2： %avg & std.
% 
% Tau_Lq_disorder_avg_std_L=zeros(6,N_q,2);
% Delta_Lq_disorder_avg_std_L=zeros(6,N_q,2);
% Gama_Lq_disorder_avg_std_L=zeros(6,N_q,2);
% Alpha_Lq_disorder_avg_std_L=zeros(6,N_q,2);
% f_alpha_Lq_disorder_avg_std_L=zeros(6,N_q,2);

%% Calculation parts Random part, data storage

delete(gcp('nocreate'));
if isempty(gcp('nocreate'))
    c=parcluster;
    Corenum=c.NumWorkers;
    threads=c.NumThreads;
    c.NumThreads = 1;
    parpool(Corenum);%36
else
    delete(gcp('nocreate'));
    c=parcluster;
    Corenum=c.NumWorkers;
    threads=c.NumThreads;
    c.NumThreads=1;
    %     c.NumWorkers=Corenum;%
    parpool(44);%36
end

% Random stream
stream=RandStream('mrg32k3a','Seed','shuffle');%'mlfg6331_64','mrg32k3a','philox4x32_10','threefry4x64_20',mldfg6331
File_name=['L',num2str(L_line(1)),'_',num2str(L_line(end)),'_Theta',num2str(floor(100*Theta_clean/pi)),...
    '_N_disorder_10_',num2str(floor(log10(N_realization_line(end)))),'_Index',num2str(n)];


%Calculations
tic
[Probaility_LDOS_all_avg_std_L,Probaility_LDOS_Bulk_avg_std_L,Probaility_LDOS_Edge_avg_std_L,...
    Probaility_LDOS_Corner_avg_std_L,Probaility_LDOS_Bulk_cylinder_avg_std_L,Probaility_LDOS_Edge_cylinder_avg_std_L,...
    Tau_Lq_disorder_avg_std_L,Delta_Lq_disorder_avg_std_L,Gama_Lq_disorder_avg_std_L,...
    Alpha_Lq_disorder_avg_std_L,f_alpha_Lq_disorder_avg_std_L]=Eigen_cavity_several(L_line(n),N_realization_line(n),Eps_center,Theta_clean,Disorder_phase_level,Disorder_Theta_level,...
    q_line,Data_N,LDOS_Range,N_bar,d_edge,Graph_off,n,stream);
toc
save(File_name,'Probaility_LDOS_all_avg_std_L','Probaility_LDOS_Bulk_avg_std_L','Probaility_LDOS_Edge_avg_std_L','Probaility_LDOS_Corner_avg_std_L','Probaility_LDOS_Bulk_cylinder_avg_std_L','Probaility_LDOS_Edge_cylinder_avg_std_L',...
    'Tau_Lq_disorder_avg_std_L','Delta_Lq_disorder_avg_std_L','Gama_Lq_disorder_avg_std_L','Alpha_Lq_disorder_avg_std_L','f_alpha_Lq_disorder_avg_std_L')
end

