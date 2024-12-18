function varargout=Single_setting_cal(Boundary_x_botton,Boundary_y_botton,SN,P_all_Bulk,Eps_center,N_q,Disorder_phase_level,q_line,...
    Port_boundary_L,Port_boundary_R,Port_boundary_T,Port_boundary_B,...
    Mx,My,Array_index,N_port,Posi_port_total_in,...
    Num_eigenstate,stream_index,stream,Graph_off,LDOS_Range,N_bar)
%Depend on the boundary condition, we output 1, 2, 3 values for probability
%and MF results
stream.Substream = stream_index;
q_line=q_line.';
L=2*Mx;
%% Boundary, Bulk, Corner recognizaton
% Based on Multifractality and Conformal Invariance at 2D Metal-Insulator Transition in the Spin-Orbit Symmetry Class
% as the eigenstate is for a, we only count a's index
% Output: Port_index_Corner,Port_index_edge, Port_index_Bulk,

%Parameters
W_corner_cell=4;
H_edge_scatter=1;% if even just the whole cell 1~8; if odd in cell index left 1234, right 5678,top 3456, bottom 1278,
L_bulk_cell_x=ceil(Mx/6);L_bulk_cell_y=ceil(My/6);

if (Boundary_x_botton==1)||(Boundary_y_botton==1)%Do not need it for double PBC or TBC
    [Port_index_Corner,Port_index_Edge,Port_index_Bulk]=Port_classification(Boundary_x_botton,Boundary_y_botton,W_corner_cell,H_edge_scatter,L_bulk_cell_x,L_bulk_cell_y,...
        Mx,My,Array_index);

    % For check:
%     figure(2)
%     scatter(Posi_port_total_in(Port_index_Bulk,1),Posi_port_total_in(Port_index_Bulk,2))
%     hold on
%     scatter(Posi_port_total_in(Port_index_Edge,1),Posi_port_total_in(Port_index_Edge,2))
%     hold on
%     scatter(Posi_port_total_in(Port_index_Corner,1),Posi_port_total_in(Port_index_Corner,2))
%     daspect([1 1 1])
end

%% Functionality I: Cavity spectra calculation

%Boundary setting
if (Boundary_x_botton==1)&&(Boundary_y_botton==1)%Double OBC
    % Pair all reflection ports:
    P_all=[Port_boundary_L;Port_boundary_R;Port_boundary_T;Port_boundary_B;P_all_Bulk];
elseif (Boundary_x_botton==1)&&(Boundary_y_botton==0)%X OBC, Y PBC
    Port_boundary_TB=[Port_boundary_T(:,1),Port_boundary_B(:,2);Port_boundary_B(:,1),Port_boundary_T(:,2)];%out from bottom linked with the corresponding port at the top input one and the in from the top linked with corresponding port at the bottom one
    P_all=[Port_boundary_L;Port_boundary_R;Port_boundary_TB;P_all_Bulk];
elseif (Boundary_x_botton==0)&&(Boundary_y_botton==1)%X PBC, Y OBC
    Port_boundary_LR=[Port_boundary_L(:,1),Port_boundary_R(:,2);Port_boundary_R(:,1),Port_boundary_L(:,2)];%out from bottom linked with the corresponding port at the top input one and the in from the top linked with corresponding port at the bottom one
    P_all=[Port_boundary_LR;Port_boundary_T;Port_boundary_B;P_all_Bulk];
elseif (Boundary_x_botton==0)&&(Boundary_y_botton==0)%X PBC, Y PBC
    Port_boundary_LR=[Port_boundary_L(:,1),Port_boundary_R(:,2);Port_boundary_R(:,1),Port_boundary_L(:,2)];%out from bottom linked with the corresponding port at the top input one and the in from the top linked with corresponding port at the bottom one
    Port_boundary_TB=[Port_boundary_T(:,1),Port_boundary_B(:,2);Port_boundary_B(:,1),Port_boundary_T(:,2)];%out from bottom linked with the corresponding port at the top input one and the in from the top linked with corresponding port at the bottom one
    P_all=[Port_boundary_LR;Port_boundary_TB;P_all_Bulk];
end

% P matrix generation
Element_phase_all=exp(1i*(rand(stream,N_port,1)-1/2)*Disorder_phase_level);
P_matrix=sparse(P_all(:,1),P_all(:,2),Element_phase_all,N_port,N_port);

% for debug
% [final,permu]=sort(P_all(:,1),'ascend');
% ZZ=P_all(permu,:);%check rank

U_evolution=P_matrix\SN;

%% I.A Several eigenstates for Multifractality check
[V_sparse,~]=eigs(U_evolution,Num_eigenstate,exp(1i*Eps_center));
%,'largestreal'); eigs(U_evolution,1,'largestreal') not work
% Eigenvalues_sparse=angle(diag(D_sparse));
if Graph_off==0
    figure(5)
    scatter(Posi_port_total_in(:,1),Posi_port_total_in(:,2),20,abs(V_sparse(:,1)),'filled');
end



%% Calculation of Probability
LDOS_bar=(LDOS_Range(2)-LDOS_Range(1))/N_bar*(0:N_bar)+LDOS_Range(1);
LDOS_bar_center=LDOS_bar(1:(end-1))-(LDOS_Range(2)-LDOS_Range(1))/N_bar/2;

if (Boundary_x_botton==1)&&(Boundary_y_botton==1)%Double OBC
    %     Probaility_LDOS=cell(1,3);
    LDOS_all_Bulk=log(abs(V_sparse(Port_index_Bulk,:)).^2*L^2);LDOS_all_Bulk=reshape(LDOS_all_Bulk,[],1);
    LDOS_all_Edge=log(abs(V_sparse(Port_index_Edge,:)).^2*L^2);LDOS_all_Edge=reshape(LDOS_all_Edge,[],1);
    LDOS_all_Corner=log(abs(V_sparse(Port_index_Corner,:)).^2*L^2);LDOS_all_Corner=reshape(LDOS_all_Corner,[],1);

    varargout{1}=P_LDOS(LDOS_all_Bulk,LDOS_bar,LDOS_Range,N_bar);%Probaility_LDOS_Bulk
    varargout{2}=P_LDOS(LDOS_all_Edge,LDOS_bar,LDOS_Range,N_bar);%Probaility_LDOS_Edge
    varargout{3}=P_LDOS(LDOS_all_Corner,LDOS_bar,LDOS_Range,N_bar);%Probaility_LDOS_Corner

    %For debug
    if Graph_off==1
        Probaility_LDOS_Bulk=varargout{1};
        Probaility_LDOS_Edge=varargout{2};
        Probaility_LDOS_Corner=varargout{3};
        figure(1)
        subplot(1,2,1)
        hold on
        subplot(1,2,2)
        %     scatter(LDOS_bar_center,Probaility_LDOS,20,'MarkerEdgeColor',mycolor('RoyalBlue'),...
        %         'MarkerFaceColor',mycolor('LightBlue'));
        %     hold on
        plot(LDOS_bar_center,Probaility_LDOS_Bulk,"-",'LineWidth',1,'Color',mycolor('ForestGreen'))
        hold on
        plot(LDOS_bar_center,Probaility_LDOS_Edge,"-",'LineWidth',1,'Color',mycolor('Pink'))
        hold on
        plot(LDOS_bar_center,Probaility_LDOS_Corner,"-",'LineWidth',1,'Color',mycolor('Gray'))
        legend('Bulk','Edge','Corner','Interpreter','latex')
        xlabel('ln $|\Psi|^2 L^2$', 'interpreter', 'latex','FontName', 'Arial','FontAngle','italic', 'FontSize', 15,'LineWidth', 1);
        ylabel('Probaility $ P$','interpreter', 'latex','FontName', 'Arial', 'FontSize', 15, 'LineWidth', 1);
        set(gca,'FontName', 'Arial', 'FontSize', 15, 'LineWidth', 1);
        title({'Probabilitydistribution functions';[ '$L$', '=', num2str(L),...
            '; $N=$',num2str(Num_eigenstate)]},...
            'interpreter', 'latex','FontName', 'Arial', 'FontSize', 15, 'LineWidth', 1)
        box on
        hold on
    end

elseif (Boundary_x_botton==1)&&(Boundary_y_botton==0)%X OBC, Y PBC
    %     Probaility_LDOS=cell(1,2);
    LDOS_all_Bulk=log(abs(V_sparse(Port_index_Bulk,:)).^2*L^2);LDOS_all_Bulk=reshape(LDOS_all_Bulk,[],1);
    LDOS_all_Edge=log(abs(V_sparse(Port_index_Edge,:)).^2*L^2);LDOS_all_Edge=reshape(LDOS_all_Edge,[],1);

    varargout{1}=P_LDOS(LDOS_all_Bulk,LDOS_bar,LDOS_Range,N_bar);%Probaility_LDOS_Bulk
    varargout{2}=P_LDOS(LDOS_all_Edge,LDOS_bar,LDOS_Range,N_bar);%Probaility_LDOS_Edge

elseif (Boundary_x_botton==0)&&(Boundary_y_botton==1)%X PBC, Y OBC
    %     Probaility_LDOS=cell(1,2);
    LDOS_all_Bulk=log(abs(V_sparse(Port_index_Bulk,:)).^2*L^2);LDOS_all_Bulk=reshape(LDOS_all_Bulk,[],1);
    LDOS_all_Edge=log(abs(V_sparse(Port_index_Edge,:)).^2*L^2);LDOS_all_Edge=reshape(LDOS_all_Edge,[],1);

    varargout{1}=P_LDOS(LDOS_all_Bulk,LDOS_bar,LDOS_Range,N_bar);%Probaility_LDOS_Bulk
    varargout{2}=P_LDOS(LDOS_all_Edge,LDOS_bar,LDOS_Range,N_bar);%Probaility_LDOS_Edge
        %For debug
    if Graph_off==1
        Probaility_LDOS_Bulk=varargout{1};
        Probaility_LDOS_Edge=varargout{2};
        figure(1)
        subplot(1,2,2)
        %     scatter(LDOS_bar_center,Probaility_LDOS,20,'MarkerEdgeColor',mycolor('RoyalBlue'),...
        %         'MarkerFaceColor',mycolor('LightBlue'));
        %     hold on
        plot(LDOS_bar_center,Probaility_LDOS_Bulk,"*",'LineWidth',1,'Color',mycolor('ForestGreen'))
        hold on
        plot(LDOS_bar_center,Probaility_LDOS_Edge,"--",'LineWidth',1,'Color',mycolor('Pink'))
        xlabel('ln $|\Psi|^2 L^2$', 'interpreter', 'latex','FontName', 'Arial','FontAngle','italic', 'FontSize', 15,'LineWidth', 1);
        ylabel('Probaility $ P$','interpreter', 'latex','FontName', 'Arial', 'FontSize', 15, 'LineWidth', 1);
        legend('Bulk Cylinder','Edge Cylinder','Interpreter','latex')
        set(gca,'FontName', 'Arial', 'FontSize', 15, 'LineWidth', 1);
        title({'Probabilitydistribution functions';[ '$L$', '=', num2str(L),...
            '; $N=$',num2str(Num_eigenstate)]},...
            'interpreter', 'latex','FontName', 'Arial', 'FontSize', 15, 'LineWidth', 1)
        box on
        hold on
    end

elseif (Boundary_x_botton==0)&&(Boundary_y_botton==0)%X PBC, Y PBC
    LDOS_all=log(abs(V_sparse).^2*L^2);LDOS_all=reshape(LDOS_all,[],1);
    varargout{1}=P_LDOS(LDOS_all,LDOS_bar,LDOS_Range,N_bar);%Probaility_LDOS_all

    %For debug
    if Graph_off==1
        Probaility_LDOS_all=varargout{1};
        figure(1)
        subplot(1,2,1)
        bar_stat=bar(LDOS_bar_center,Probaility_LDOS_all,'FaceColor','flat','BarWidth', 1);bar_stat.CData=Probaility_LDOS_all.';
        colormap(gca,mymap('Blues'))
        caxis(gca,[0, max(Probaility_LDOS_all)])%min(P_state_statistic), max(P_state_statistic)
        axis([LDOS_Range(1), LDOS_Range(2), 0, max(Probaility_LDOS_all)])
        xlabel('ln $|\Psi|^2 L^2$', 'interpreter', 'latex','FontName', 'Arial','FontAngle','italic', 'FontSize', 10,'LineWidth', 1);
        ylabel('Probaility $ P$','interpreter', 'latex','FontName', 'Arial', 'FontSize', 10, 'LineWidth', 1);
        legend('Total Bulk','Interpreter','latex')
        set(gca,'FontName', 'Arial', 'FontSize', 15, 'LineWidth', 1);
        title({'LDOS statistics Bulk (2D PBC)';[ '$L$', '=', num2str(L),...
            '; $N=$',num2str(Num_eigenstate)]}, 'interpreter', 'latex','FontName', 'Arial', 'FontSize', 15, 'LineWidth', 1)
        colorbar
        box on
        hold on
        subplot(1,2,2)
        %     scatter(LDOS_bar_center,Probaility_LDOS,20,'MarkerEdgeColor',mycolor('RoyalBlue'),...
        %         'MarkerFaceColor',mycolor('LightBlue'));
        %     hold on
        plot(LDOS_bar_center,Probaility_LDOS_all,"--",'LineWidth',1,'Color',mycolor('RoyalBlue'))
        hold on
    end
end


%% Calculation of Multifractality
%MF exponents Tau, MF singularity spectra
Mean_Psi_2q=@(Psi,q)(mean(abs(Psi).^(2*q),'all'));
%Legendre transform MF singularity spectra
Mean_Psi_alpha_q_top=@(Psi,q)(mean(abs(Psi).^(2*q).*log(abs(Psi).^2),'all'));

% (Port_index_Bulk,:)
% (Port_index_Edge,:)
% (Port_index_Corner,:)

if (Boundary_x_botton==1)&&(Boundary_y_botton==1)%Double OBC
 %对tau来说，本次无序平均是在log之前，而不是每次的总平均。
 %Obuse的文章中的是拟合的不确定性。而不是平均值的不确定性。原因是Boundary Multifractality at the Integer Quantum Hall Plateau Transition: Implications for the Critical Theory
 %一文的Fig.1作为每个单独计算的点，并无variation bar，因此，计算中的LDOS平均是对所有无序下的统一平均。
    varargout{4}=arrayfun(@(q) Mean_Psi_2q(V_sparse(Port_index_Bulk,:), q), q_line);%L^-dx(1)* exp(-Tau_q_disorder_Bulk*log(L))
    varargout{5}=arrayfun(@(q) Mean_Psi_2q(V_sparse(Port_index_Edge,:), q), q_line);%L^-dx(1)* exp(-Tau_q_disorder_Edge*log(L))
    varargout{6}=arrayfun(@(q) Mean_Psi_2q(V_sparse(Port_index_Corner,:), q), q_line);%L^-dx(1)* exp(-Tau_q_disorder_Corner*log(L))

    varargout{7}=arrayfun(@(q) Mean_Psi_alpha_q_top(V_sparse(Port_index_Bulk,:), q), q_line);%Numerator of alpha_Lq
    varargout{8}=arrayfun(@(q) Mean_Psi_alpha_q_top(V_sparse(Port_index_Edge,:), q), q_line);%
    varargout{9}=arrayfun(@(q) Mean_Psi_alpha_q_top(V_sparse(Port_index_Corner,:), q), q_line);%

    varargout{10}=[length(Port_index_Bulk),length(Port_index_Edge),length(Port_index_Corner)]';

elseif (Boundary_x_botton==1)&&(Boundary_y_botton==0)%X OBC, Y PBC
    varargout{3}=arrayfun(@(q) Mean_Psi_2q(V_sparse(Port_index_Bulk,:), q), q_line);%L^-dx(1)* exp(-Tau_q_disorder_Bulk*log(L))
    varargout{4}=arrayfun(@(q) Mean_Psi_2q(V_sparse(Port_index_Edge,:), q), q_line);%L^-dx(1)* exp(-Tau_q_disorder_Edge*log(L))

    varargout{5}=arrayfun(@(q) Mean_Psi_alpha_q_top(V_sparse(Port_index_Bulk,:), q), q_line);%Numerator of alpha_Lq
    varargout{6}=arrayfun(@(q) Mean_Psi_alpha_q_top(V_sparse(Port_index_Edge,:), q), q_line);

    varargout{7}=[length(Port_index_Bulk),length(Port_index_Edge)]';
elseif (Boundary_x_botton==0)&&(Boundary_y_botton==1)%X PBC, Y OBC
    varargout{3}=arrayfun(@(q) Mean_Psi_2q(V_sparse(Port_index_Bulk,:), q), q_line);%L^-dx(1)* exp(-Tau_q_disorder_Bulk*log(L))
    varargout{4}=arrayfun(@(q) Mean_Psi_2q(V_sparse(Port_index_Edge,:), q), q_line);%L^-dx(1)* exp(-Tau_q_disorder_Edge*log(L))

    varargout{5}=arrayfun(@(q) Mean_Psi_alpha_q_top(V_sparse(Port_index_Bulk,:), q), q_line);%Numerator of alpha_Lq
    varargout{6}=arrayfun(@(q) Mean_Psi_alpha_q_top(V_sparse(Port_index_Edge,:), q), q_line);

    varargout{7}=[length(Port_index_Bulk),length(Port_index_Edge)]';

elseif (Boundary_x_botton==0)&&(Boundary_y_botton==0)%X PBC, Y PBC
    %      bsxfun(@Mean_Psi_2q,V_sparse,q_line)
    % result = arrayfun(@(q) Mean_Psi_2q(V_sparse, q), q_line);
    varargout{2}=arrayfun(@(q) Mean_Psi_2q(V_sparse, q), q_line);%L^-dx(1)* exp(-Tau_q_disorder_Bulk*log(L)), plot with L check intercept
    %     (varargout{2}-2*q_line+dx(1))./(q_line.*(q_line-1))
    varargout{3}=arrayfun(@(q) Mean_Psi_alpha_q_top(V_sparse, q), q_line);%%Numerator of alpha_Lq
    varargout{4}=N_port;
end

end

function Probaility_LDOS=P_LDOS(LDOS,LDOS_bar,LDOS_Range,N_bar)
Probaility_LDOS=zeros(N_bar,1);
for i=1:N_bar
    [row,~]=find((LDOS<LDOS_bar(i+1))&(LDOS>=LDOS_bar(i)));
    Probaility_LDOS(i)=length(row);
end
Probaility_LDOS=(Probaility_LDOS/length(LDOS))/((LDOS_Range(2)-LDOS_Range(1))/N_bar);
end
