function [Port_index_Corner,Port_index_edge, Port_index_Bulk]=Port_classification(Boundary_x_botton,Boundary_y_botton,W_corner_cell,H_edge_scatter,L_bulk_cell_x,L_bulk_cell_y,...
    Mx,My,Array_index)
% Parameter: Boundary condition(Boundary_x_botton,Boundary_y_botton),W_corner_cell,H_edge_scatter,L_bulk_cell_x,L_bulk_cell_y
% Input : Mx,My,Array_index
% Output: Port_index_Corner,Port_index_edge, Port_index_Bulk,

H_edge_cell=floor(H_edge_scatter/2);

if (Boundary_x_botton==1)&&(Boundary_y_botton==1)%Double OBC
    L_cell_range_x=floor((Mx-L_bulk_cell_x)/2):ceil((Mx+L_bulk_cell_x)/2);%square, for both x and y directions
    L_cell_range_y=floor((My-L_bulk_cell_y)/2):ceil((My+L_bulk_cell_y)/2);%square, for both x and y directions
elseif (Boundary_x_botton==1)&&(Boundary_y_botton==0)%X OBC, Y PBC
    L_cell_range_x=floor((Mx-L_bulk_cell_x)/2):ceil((Mx+L_bulk_cell_x)/2);%square, for both x and y directions
    L_cell_range_y=1:My;%square, for both x and y directions
elseif (Boundary_x_botton==0)&&(Boundary_y_botton==1)%X PBC, Y OBC
    L_cell_range_x=1:Mx;%square, for both x and y directions
    L_cell_range_y=floor((My-L_bulk_cell_y)/2):ceil((My+L_bulk_cell_y)/2);%square, for both x and y directions
end
%Corner_part, only work for Double OBC
if (Boundary_x_botton==1)&&(Boundary_y_botton==1)%Double OBC
    Rc_line_corner_LB=reshape(Array_index(1:W_corner_cell,1:W_corner_cell),[],1);%left edge
    Rc_line_corner_LT=reshape(Array_index(1:W_corner_cell,(My-W_corner_cell+1):My),[],1);
    Rc_line_corner_RB=reshape(Array_index((Mx-W_corner_cell+1):Mx,1:W_corner_cell),[],1);%right edge
    Rc_line_corner_RT=reshape(Array_index((Mx-W_corner_cell+1):Mx,(My-W_corner_cell+1):My),[],1);

    Rc_line_corner=sort([Rc_line_corner_LB;Rc_line_corner_LT;Rc_line_corner_RB;Rc_line_corner_RT],'ascend');%-1 is important
    Port_index_Corner=kron((Rc_line_corner-1)*8,ones(8,1))+kron(ones(length(Rc_line_corner),1),(1:8)');
else
    Port_index_Corner=[];
end

%Edge part
% if even just the whole cell 1~8;
if H_edge_cell>0
    Rc_line_edge_left=reshape(Array_index(1:H_edge_cell,L_cell_range_y),[],1);%left edge
    Rc_line_edge_bot=reshape(Array_index(L_cell_range_x,1:H_edge_cell),[],1);
    Rc_line_edge_right=reshape(Array_index((Mx-H_edge_cell+1):Mx,L_cell_range_y),[],1);%right edge
    Rc_line_edge_top=reshape(Array_index(L_cell_range_x,(My-H_edge_cell+1):My),[],1);
    if (Boundary_x_botton==1)&&(Boundary_y_botton==1)%Double OBC
        Rc_line_edge=sort([Rc_line_edge_left;Rc_line_edge_bot;Rc_line_edge_right;Rc_line_edge_top],'ascend');
    elseif (Boundary_x_botton==1)&&(Boundary_y_botton==0)%X OBC, Y PBC
        Rc_line_edge=sort([Rc_line_edge_left;Rc_line_edge_right],'ascend');
    elseif (Boundary_x_botton==0)&&(Boundary_y_botton==1)%X PBC, Y OBC
        Rc_line_edge=sort([Rc_line_edge_bot;Rc_line_edge_top],'ascend');
    end
    Port_index_edge_1=kron((Rc_line_edge-1)*8,ones(8,1))+kron(ones(length(Rc_line_edge),1),(1:8)');
else
    Port_index_edge_1=[];
end

%if odd in cell index, there are still ports occupy half of the cell with the index: left edge 1234, right edge 5678,top edge 3456, bottom edge 1278, This is the residual part
if mod(H_edge_scatter,2)==1
    Port_line_edge_left_residual=reshape(Array_index(H_edge_cell+1,L_cell_range_y),[],1);
    Port_line_edge_left_residual=kron((Port_line_edge_left_residual-1)*8,ones(4,1))+kron(ones(length(Port_line_edge_left_residual),1),(1:4)');

    Port_line_edge_right_residual=reshape(Array_index(Mx-H_edge_cell,L_cell_range_y),[],1);
    Port_line_edge_right_residual=kron((Port_line_edge_right_residual-1)*8,ones(4,1))+kron(ones(length(Port_line_edge_right_residual),1),(5:8)');

    Port_line_edge_top_residual=reshape(Array_index(L_cell_range_x,My-H_edge_cell),[],1);
    Port_line_edge_top_residual=kron((Port_line_edge_top_residual-1)*8,ones(4,1))+kron(ones(length(Port_line_edge_top_residual),1),(3:6)');

    Port_line_edge_bot_residual=reshape(Array_index(L_cell_range_x,H_edge_cell+1),[],1);
    Port_line_edge_bot_residual=kron((Port_line_edge_bot_residual-1)*8,ones(4,1))+kron(ones(length(Port_line_edge_bot_residual),1),[1 2 7 8]');

    if (Boundary_x_botton==1)&&(Boundary_y_botton==1)%Double OBC
        Port_index_edge_2=[Port_line_edge_left_residual;Port_line_edge_right_residual;Port_line_edge_top_residual;Port_line_edge_bot_residual];
    elseif (Boundary_x_botton==1)&&(Boundary_y_botton==0)%X OBC, Y PBC
        Port_index_edge_2=[Port_line_edge_left_residual;Port_line_edge_right_residual];
    elseif (Boundary_x_botton==0)&&(Boundary_y_botton==1)%X PBC, Y OBC
        Port_index_edge_2=[Port_line_edge_top_residual;Port_line_edge_bot_residual];
    end
else
    Port_index_edge_2=[];
end
Port_index_edge=sort([Port_index_edge_1;Port_index_edge_2],'ascend');

%Bulk part
%For Double PBC, we just take all the ports
if (Boundary_x_botton==1)&&(Boundary_y_botton==1)%Double OBC
    Rc_line_Bulk=sort(reshape(Array_index(L_cell_range_x,L_cell_range_y),[],1),'ascend');
elseif (Boundary_x_botton==1)&&(Boundary_y_botton==0)%X OBC, Y PBC
    Rc_line_Bulk=sort(reshape(Array_index(L_cell_range_x,1:My),[],1),'ascend');
elseif (Boundary_x_botton==0)&&(Boundary_y_botton==1)%X PBC, Y OBC
    Rc_line_Bulk=sort(reshape(Array_index(1:Mx,L_cell_range_y),[],1),'ascend');
end
Port_index_Bulk=kron((Rc_line_Bulk-1)*8,ones(8,1))+kron(ones(length(Rc_line_Bulk),1),(1:8)');

end