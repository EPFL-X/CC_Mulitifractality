
function [Probaility_LDOS_all_avg_std,Probaility_LDOS_Bulk_avg_std,Probaility_LDOS_Edge_avg_std,...
     Probaility_LDOS_Corner_avg_std,Probaility_LDOS_Bulk_cylinder_avg_std,Probaility_LDOS_Edge_cylinder_avg_std,...
     Tau_Lq_disorder_avg_std,Delta_Lq_disorder_avg_std,Gama_Lq_disorder_avg_std,...
     Alpha_Lq_disorder_avg_std,f_alpha_Lq_disorder_avg_std]=Eigen_cavity_several(L,N_realization,...
    Eps_center,Theta_clean,Disorder_phase_level,Disorder_Theta_level,...
    q_line,Data_N,LDOS_Range,N_bar,d_edge,Graph_off,stream_index,stream)

N_q=Data_N(1);Num_eigenstate=Data_N(2);

%% Size setting
Mx=L;My=L;%Mx My are the number of unitcell for x direction and y direction
dL=1/max([Mx,My]);% Size of unitcell, keep in the range of [0 1 0 1], prefactor for all the position data in the ploting
N_port=Mx*My*8;
N_node=Mx*My*4;
%% Positions of unitcells, nodes, and ports
% Posi_node=[(ix-1),(iy-1)];%Left bottom corner coordination
Posi_node_incell=[1/4,1/4;1/4,3/4;3/4,3/4;3/4,1/4];
Posi_port_incell_inwave=[0,1/4;1/2,1/4;1/4,1/2;1/4,1;...
    1/2,3/4;1,3/4;3/4,0;3/4,1/2];
Posi_port_incell_outwave=[1/4,1/2;1/4,0;0,3/4;1/2,3/4;...
    3/4,1;3/4,1/2;1/2,1/4;1,1/4];
Rc_line=(1:Mx*My)';iy_line=kron(ones(Mx,1),(1:My)');ix_line=kron((1:Mx)',ones(My,1));
Array_index=reshape(Rc_line,My,Mx)';%Rc=Array_index(ix,iy);
Posi_unitcell_total=[(ix_line-1),(iy_line-1)];%Left bottom corner coordination, shift 1/2 and 1/2 is the center of unitcell
% Posi_node_total=kron(Posi_unitcell_total,ones(4,1))+kron(ones(Mx*My,1),Posi_node_incell);

% Position of a and b can be different
Posi_port_total_in=kron(Posi_unitcell_total,ones(8,1))+kron(ones(Mx*My,1),Posi_port_incell_inwave);
% Posi_port_total_out=kron(Posi_unitcell_total,ones(8,1))+kron(ones(Mx*My,1),Posi_port_incell_outwave);

% For check
% figure(1)
% scatter(Posi_node_total(:,1),Posi_node_total(:,2))
% hold on
% scatter(Posi_unitcell_total(:,1)+1/2,Posi_unitcell_total(:,2)+1/2)
% hold on
% scatter(Posi_port_total_in(:,1),Posi_port_total_in(:,2))

%% Connection matrix: |b>=P|a>e^(-i\phi),|a> incoming wave, |b> outgoing wave
% For fully randomed \phi, P's phase should be randomlized.
% For acceleration, using spharse matrix calculations.

% Connections intra Unitcell
P_incell=[1,3;4,5;6,8;7,2];
P_incell_all=kron((Rc_line-1)*8,ones(4,2))+kron(ones(Mx*My,1),P_incell);%connection data, 4*Mx*My,from the second column one (OUT) to the first column one (IN)

% Connections inter Unitcell: this should be careful
% Two cases: x and y connections
% 1 ix=1~Mx-1& iy=1~My, take x connections
% x connection: b_(Rc+My,3)=a_(Rc,6) e^(-1i*\phi) &  b_(Rc,8)=a_(Rc+My,1) e^(-1i*\phi)

% 2 ix=1~Mx & iy=1~My-1, take y connections
% y connection: b_(Rc+1,2)=a_(Rc,4) e^(-1i*\phi) &  b_(Rc,5)=a_(Rc+1,7) e^(-1i*\phi)

Rc_line_inrange_1=sort(reshape(Array_index(1:(Mx-1),1:My),[],1),'ascend');%Array_index(2:(Mx-1),2:(My-1));
P_1_x=[[(Rc_line_inrange_1+My-1)*8+3,(Rc_line_inrange_1-1)*8+6];[(Rc_line_inrange_1-1)*8+8,(Rc_line_inrange_1+My-1)*8+1]];

% 2 ix=1 and Mx & iy=2~My-1
% only y connections
Rc_line_inrange_2=sort(reshape(Array_index(1:Mx,1:(My-1)),[],1),'ascend');%Array_index(2:(Mx-1),2:(My-1));
P_2_y=[[(Rc_line_inrange_2+1-1)*8+2,(Rc_line_inrange_2-1)*8+4];[(Rc_line_inrange_2-1)*8+5,(Rc_line_inrange_2+1-1)*8+7]];

% Total inter unitcell connection
P_intercell_all=[P_1_x;P_2_y];%4*Mx*My-2*Mx-2*My
% Total networks' bulk connections
P_all=[P_incell_all;P_intercell_all];
P_all_C=parallel.pool.Constant(P_all);
%% Boundary ports
%  which is going to be applied with boundary conditions
% For cavity case, it is full reflections;
% For finite transport: most of them fully reflected but some are external (open ports);
% For band structure or supercell: they can be TBC/Flux in them.
% Therefore, we need to know, the rank for in and out, also, separate the top and bottom, left and right, put in cell format

%Boundary detection
%Left right input and output port index from bottom to top
Port_boundary_L=[reshape(Array_index(1,1:My)-1,[],1)*8+3,reshape(Array_index(1,1:My)-1,[],1)*8+1];%first column out and second column in
Port_boundary_R=[reshape(Array_index(Mx,1:My)-1,[],1)*8+8,reshape(Array_index(Mx,1:My)-1,[],1)*8+6];

%Top bottom input and output port index from left to right
Port_boundary_T=[reshape(Array_index(1:Mx,My)-1,[],1)*8+5,reshape(Array_index(1:Mx,My)-1,[],1)*8+4];%first column out and second column in
Port_boundary_B=[reshape(Array_index(1:Mx,1)-1,[],1)*8+2,reshape(Array_index(1:Mx,1)-1,[],1)*8+7];

%% Scattering part for only the scattering itselft formation
% Be careful, if there is no disorder on scattering node, it is fine.
Theta_all=Theta_clean*ones(N_node,1)+(rand(N_node,1)-1/2)*Disorder_Theta_level;
SU2_cell=cell(N_node,1);
parfor i=1:N_node
    SU2_cell{i}=sparse(SU2(Theta_all(i)));
end
SN=blkdiag(SU2_cell{:});
SN_C=parallel.pool.Constant(SN);
%% Key calculations this should be calculate in another code.Boundary condition is set there
% Data prepare
Probaility_LDOS_all=zeros(N_bar,N_realization);Probaility_LDOS_Bulk=zeros(N_bar,N_realization);
Probaility_LDOS_Edge=zeros(N_bar,N_realization);Probaility_LDOS_Corner=zeros(N_bar,N_realization);
Probaility_LDOS_Bulk_cylinder=zeros(N_bar,N_realization);Probaility_LDOS_Edge_cylinder=zeros(N_bar,N_realization);

%There is no L in the name due to this for fixed L with different realizations.
% Tau_q_disorder,Alpha_q_disorder,f_alpha_q_disorder
Mean_Psi_2q_all=zeros(N_q,N_realization);Mean_Psi_2q_Bulk=zeros(N_q,N_realization);
Mean_Psi_2q_Edge=zeros(N_q,N_realization);Mean_Psi_2q_Corner=zeros(N_q,N_realization);
Mean_Psi_2q_Bulk_cylinder=zeros(N_q,N_realization);Mean_Psi_2q_Edge_cylinder=zeros(N_q,N_realization);

Mean_Psi_alpha_q_top_all=zeros(N_q,N_realization);Mean_Psi_alpha_q_top_Bulk=zeros(N_q,N_realization);
Mean_Psi_alpha_q_top_Edge=zeros(N_q,N_realization);Mean_Psi_alpha_q_top_Corner=zeros(N_q,N_realization);
Mean_Psi_alpha_q_top_Bulk_cylinder=zeros(N_q,N_realization);Mean_Psi_alpha_q_top_Edge_cylinder=zeros(N_q,N_realization);

L_2DPBC=zeros(1,N_realization);%_all,_Bulk,_Edge,_Corner,Bulk_cylinder,Edge_Cylinder 
L_2DOBC=zeros(3,N_realization);
L_1DPBC=zeros(2,N_realization);

parfor j=1:N_realization
%     tic
    [Probaility_LDOS_all(:,j),...
        Mean_Psi_2q_all(:,j),...
        Mean_Psi_alpha_q_top_all(:,j),...
        L_2DPBC(j)]=Single_setting_cal(0,0,SN_C.Value,P_all_C.Value,Eps_center,N_q,Disorder_phase_level,q_line,...
        Port_boundary_L,Port_boundary_R,Port_boundary_T,Port_boundary_B,...
        Mx,My,Array_index,N_port,Posi_port_total_in,...
        Num_eigenstate,3*(N_realization*(stream_index-1)+j)+1,stream,Graph_off,LDOS_Range,N_bar);

%     [Probaility_LDOS_Bulk(:,j),Probaility_LDOS_Edge(:,j),Probaility_LDOS_Corner(:,j),...
%         Mean_Psi_2q_Bulk(:,j),Mean_Psi_2q_Edge(:,j),Mean_Psi_2q_Corner(:,j),...
%         Mean_Psi_alpha_q_top_Bulk(:,j),Mean_Psi_alpha_q_top_Edge(:,j),Mean_Psi_alpha_q_top_Corner(:,j),...
%         L_2DOBC(:,j)]=Single_setting_cal(1,1,SN_C.Value,P_all_C.Value,Eps_center,N_q,Disorder_phase_level,q_line,...
%         Port_boundary_L,Port_boundary_R,Port_boundary_T,Port_boundary_B,...
%         Mx,My,Array_index,N_port,Posi_port_total_in,...
%         Num_eigenstate,3*(N_realization*(stream_index-1)+j)+2,stream,Graph_off,LDOS_Range,N_bar);

    [Probaility_LDOS_Bulk_cylinder(:,j),Probaility_LDOS_Edge_cylinder(:,j),...
        Mean_Psi_2q_Bulk_cylinder(:,j),Mean_Psi_2q_Edge_cylinder(:,j),...
        Mean_Psi_alpha_q_top_Bulk_cylinder(:,j),Mean_Psi_alpha_q_top_Edge_cylinder(:,j),...
        L_1DPBC(:,j)]=Single_setting_cal(0,1,SN_C.Value,P_all_C.Value,Eps_center,N_q,Disorder_phase_level,q_line,...
        Port_boundary_L,Port_boundary_R,Port_boundary_T,Port_boundary_B,...
        Mx,My,Array_index,N_port,Posi_port_total_in,...
        Num_eigenstate,3*(N_realization*(stream_index-1)+j)+3,stream,Graph_off,LDOS_Range,N_bar);
%     toc
end

%Probability statistics
%As the statistics in theory is over all the realizations.
Probaility_LDOS_all_avg_std=zeros(N_bar,2);
Probaility_LDOS_all_avg_std(:,1)=mean(Probaility_LDOS_all,2);
Probaility_LDOS_all_avg_std(:,2)=std(Probaility_LDOS_all,1,2);

Probaility_LDOS_Bulk_avg_std=zeros(N_bar,2);
Probaility_LDOS_Bulk_avg_std(:,1)=mean(Probaility_LDOS_Bulk,2);
Probaility_LDOS_Bulk_avg_std(:,2)=std(Probaility_LDOS_Bulk,1,2);

Probaility_LDOS_Edge_avg_std=zeros(N_bar,2);
Probaility_LDOS_Edge_avg_std(:,1)=mean(Probaility_LDOS_Edge,2);
Probaility_LDOS_Edge_avg_std(:,2)=std(Probaility_LDOS_Edge,1,2);

Probaility_LDOS_Corner_avg_std=zeros(N_bar,2);
Probaility_LDOS_Corner_avg_std(:,1)=mean(Probaility_LDOS_Corner,2);
Probaility_LDOS_Corner_avg_std(:,2)=std(Probaility_LDOS_Corner,1,2);

Probaility_LDOS_Bulk_cylinder_avg_std=zeros(N_bar,2);
Probaility_LDOS_Bulk_cylinder_avg_std(:,1)=mean(Probaility_LDOS_Bulk_cylinder,2);
Probaility_LDOS_Bulk_cylinder_avg_std(:,2)=std(Probaility_LDOS_Bulk_cylinder,1,2);

Probaility_LDOS_Edge_cylinder_avg_std=zeros(N_bar,2);
Probaility_LDOS_Edge_cylinder_avg_std(:,1)=mean(Probaility_LDOS_Edge_cylinder,2);
Probaility_LDOS_Edge_cylinder_avg_std(:,2)=std(Probaility_LDOS_Edge_cylinder,1,2);

% Operation and Store for output
%Store for operations
L_all_num_of_LDOS=[L_2DPBC;L_2DOBC;L_1DPBC];%number of LDOS calculated in the mean.
Mean_Psi_2q={Mean_Psi_2q_all,Mean_Psi_2q_Bulk,Mean_Psi_2q_Edge,Mean_Psi_2q_Corner,Mean_Psi_2q_Bulk_cylinder,Mean_Psi_2q_Edge_cylinder};
Mean_Psi_alpha_q_top={Mean_Psi_alpha_q_top_all,Mean_Psi_alpha_q_top_Bulk,Mean_Psi_alpha_q_top_Edge,...
    Mean_Psi_alpha_q_top_Corner,Mean_Psi_alpha_q_top_Bulk_cylinder,Mean_Psi_alpha_q_top_Edge_cylinder};

%Generate
Tau_Lq_disorder_avg_std=zeros(6,N_q,2);
Delta_Lq_disorder_avg_std=zeros(6,N_q,2);
Gama_Lq_disorder_avg_std=zeros(6,N_q,2);
Alpha_Lq_disorder_avg_std=zeros(6,N_q,2);
f_alpha_Lq_disorder_avg_std=zeros(6,N_q,2);

Tau2Delta=q_line'*ones(1,N_realization);
Delta2Gama=(q_line.*(1-q_line))'*ones(1,N_realization);

L_r=2*L;
%这里的标准差更像是，每次无序下，进行单独计算，然后进行多次统计，这样很可能L增加确实标准差变小，但是很可能是很多次的随机可以认为是完全独立的
%意味着每个小的体系单独计算，事实上，之前的标准差是这么算的。但是很可能不对，事实上，统计的误差，Tau等，应该是拟合误差。用那个做bar，而这个标准差基本上没用
%L or 2L不会影响拟合的截距。
for p=1:6
    Total_mean_2q=(sum(Mean_Psi_2q{p}.*kron(ones(N_q,1),L_all_num_of_LDOS(p,:)),2)/sum(L_all_num_of_LDOS(p,:)));
    Total_mean_Psi_alpha_q_top=sum(Mean_Psi_alpha_q_top{p}.*kron(ones(N_q,1),L_all_num_of_LDOS(p,:)),2)/sum(L_all_num_of_LDOS(p,:));

    Mean_2q_loop=zeros(N_q,1);
    Mean_Psi_2q_temp=Mean_Psi_2q{p};
    for i=1:N_q
        Mean_2q_loop(i)=sum(Mean_Psi_2q_temp(i,:).*L_all_num_of_LDOS(p,:),2)/sum(L_all_num_of_LDOS(p,:));
    end
    Error_check_2q=(Mean_2q_loop-Total_mean_2q)./Mean_2q_loop;

    Total_mean_top_loop=zeros(N_q,1);
    Mean_Psi_2q_top_temp=Mean_Psi_alpha_q_top{p};
    for i=1:N_q
        Total_mean_top_loop(i)=sum(Mean_Psi_2q_top_temp(i,:).*L_all_num_of_LDOS(p,:),2)/sum(L_all_num_of_LDOS(p,:));
    end   
    Error_check_alpha_top=(Total_mean_top_loop-Total_mean_Psi_alpha_q_top)./Total_mean_top_loop;
    %Using loop
    Total_mean_Psi_alpha_q_top=Total_mean_top_loop;
    Total_mean_2q=Mean_2q_loop;

    Tau_Lq_disorder_avg_std(p,:,1)=-log(Total_mean_2q*(L_r)^d_edge(p))/log(L_r);% mean use total mean
    Tau_Lq_disorder_avg_std(p,:,2)=std(-log(Mean_Psi_2q{p}*(L_r)^d_edge(p))/log(L_r),1,2);%Can we get the stand deviation, std
    Tau_Diff_avg_f_seq=(Tau_Lq_disorder_avg_std(p,:,1).'-mean(-log(Mean_Psi_2q{p}*(L_r)^d_edge(p))/log(L_r),2))./(Tau_Lq_disorder_avg_std(p,:,1).');
%This diff check is to check when should average be done.

    Delta_Lq_disorder_avg_std(p,:,1)=Tau_Lq_disorder_avg_std(p,:,1)-2*q_line+d_edge(p);%use total mean
    Delta_Lq_disorder_avg_std(p,:,2)=std(-log(Mean_Psi_2q{p}*(L_r)^d_edge(p))/log(L_r)-2*Tau2Delta+d_edge(p),1,2);%Can we get the stand deviation
    Delta_Diff_avg_f_seq=(Delta_Lq_disorder_avg_std(p,:,1).'-mean(-log(Mean_Psi_2q{p}*(L_r)^d_edge(p))/log(L_r)-2*Tau2Delta+d_edge(p),2))./Delta_Lq_disorder_avg_std(p,:,1).';

    Gama_Lq_disorder_avg_std(p,:,1)=Delta_Lq_disorder_avg_std(p,:,1)./(q_line.*(1-q_line));
    Gama_Lq_disorder_avg_std(p,:,2)=std((-log(Mean_Psi_2q{p}*(L_r)^d_edge(p))/log(L_r)-2*Tau2Delta+d_edge(p))./Delta2Gama,1,2);%Can we get the stand deviation
    Gama_Diff_avg_f_seq=(Gama_Lq_disorder_avg_std(p,:,1).'-mean((-log(Mean_Psi_2q{p}*(L_r)^d_edge(p))/log(L_r)-2*Tau2Delta+d_edge(p))./Delta2Gama,2))./(Gama_Lq_disorder_avg_std(p,:,1).');

    Alpha_Lq_disorder_avg_std(p,:,1)=-(Total_mean_Psi_alpha_q_top./Total_mean_2q)/log(L_r);
    Alpha_Lq_disorder_avg_std(p,:,2)=std(-(Mean_Psi_alpha_q_top{p}./Mean_Psi_2q{p})/log(L_r),1,2);
    Alpha_Diff_avg_f_seq=(Alpha_Lq_disorder_avg_std(p,:,1).'-mean(-(Mean_Psi_alpha_q_top{p}./Mean_Psi_2q{p})/log(L_r),2))./(Alpha_Lq_disorder_avg_std(p,:,1).');

    f_alpha_Lq_disorder_avg_std(p,:,1)=(log(Total_mean_2q)/log(L_r)).'+Alpha_Lq_disorder_avg_std(p,:,1).*q_line+d_edge(p);
    f_alpha_Lq_disorder_avg_std(p,:,2)=std(log(Mean_Psi_2q{p})/log(L_r)+(Mean_Psi_alpha_q_top{p}./Mean_Psi_2q{p}).*Tau2Delta+d_edge(p),1,2);
    f_alpha_Diff_avg_f_seq=(f_alpha_Lq_disorder_avg_std(p,:,1).'-mean(log(Mean_Psi_2q{p})/log(L_r)+(Mean_Psi_alpha_q_top{p}./Mean_Psi_2q{p}).*Tau2Delta+d_edge(p),2))./(f_alpha_Lq_disorder_avg_std(p,:,1).');
end

end

function S_U2=SU2(Theta)
S_U2=[cos(Theta),1i*sin(Theta);1i*sin(Theta),cos(Theta)];
%S_U2=[cos(Theta),sin(Theta);-sin(Theta),cos(Theta)]; %it can be like this.
end
