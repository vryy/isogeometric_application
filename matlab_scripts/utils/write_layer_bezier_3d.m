function write_layer_bezier_3d(fid, name, nurbs, C, Cxi, Cet, Cze, params)
%
% write the Bezier geometry info for 1 layer for 1 nurbs volume patch
%
% Input:
%   fid     file handler
%   name    layer name
%   nurbs   nurbs geometry
%   C       Bezier extraction operator
%   Cxi     Bezier extraction operator - xi
%   Cet     Bezier extraction operator - eta
%   Cze     Bezier extraction operator - zeta
%   params  additional parameters

%%
sizes = size(nurbs.coefs);
u_dim = sizes(2);
v_dim = sizes(3);
w_dim = sizes(4);

%check if n >= p + 1
n1 = u_dim;
p1 = nurbs.order(1) - 1;
if(n1 < p1 + 1)
    error('n1 < p1 + 1: not enough Bezier elements')
end

n2 = v_dim;
p2 = nurbs.order(2) - 1;
if(n2 < p2 + 1)
    error('n2 < p2 + 1: not enough Bezier elements')
end

n3 = w_dim;
p3 = nurbs.order(3) - 1;
if(n3 < p3 + 1)
    error('n3 < p3 + 1: not enough Bezier elements')
end

ne1 = size(Cxi, 3); %number of element along xi direction
ne2 = size(Cet, 3); %number of element along eta direction
ne3 = size(Cze, 3); %number of element along zeta direction
ne = size(C, 3);

if ne1 * ne2 * ne3 ~= ne
    error('invalid number of elements along 3 directions')
end

fprintf(fid, '\t\t## begin layer_info for layer %s\n', name);
fprintf(fid, '\t\tself.layer_list.append(''%s'')\n\n', name);

%%
fprintf(fid, '\t\t# nodal info\n');
fprintf(fid, '\t\tcurrent_nodal_set = {}\n');
cnt = 1;
for i = 1:w_dim
    for j = 1:v_dim
        for k = 1:u_dim
            point = nurbs.coefs(:, k, j, i);
            point(1:3) = point(1:3) / point(4);
            fprintf(fid, '\t\tcurrent_nodal_set[%d] = [%.10f, %.10f, %.10f]\n', cnt, point(1), point(2), point(3));
            cnt = cnt + 1;
        end
    end
end
fprintf(fid, '\t\tself.layer_nodes_sets[''%s''] = current_nodal_set\n\n', name);
fprintf('number of printed nodes = %d\n', cnt - 1);

%%
fprintf(fid, '\t\t# entity connectivities\n');
fprintf(fid, '\t\tcurrent_entity_set = {}\n');
cnt = 1;
% prop = 1;
b1 = p1 + 1;
sum_mul1 = 0;
for i = 1:ne1
    
    %check the multiplicity
    tmp = b1;
    while b1 <= (n1 + p1 + 1) && nurbs.knots{1}(b1+1) == nurbs.knots{1}(b1)
        b1 = b1 + 1;
    end
    mul1 = b1 - tmp + 1;
    b1 = b1 + 1;
    sum_mul1 = sum_mul1 + (mul1 - 1);
    
    b2 = p2 + 1;
    sum_mul2 = 0;
    for j = 1:ne2
        
        %check the multiplicity
        tmp = b2;
        while b2 <= (n2 + p2 + 1) && nurbs.knots{2}(b2+1) == nurbs.knots{2}(b2)
            b2 = b2 + 1;
        end
        mul2 = b2 - tmp + 1;
        b2 = b2 + 1;
        sum_mul2 = sum_mul2 + (mul2 - 1);
        
        b3 = p3 + 1;
        sum_mul3 = 0;
        for k = 1:ne3
            
            %check the multiplicity
            tmp = b3;
            while b3 <= (n3 + p3 + 1) && nurbs.knots{3}(b3+1) == nurbs.knots{3}(b3)
                b3 = b3 + 1;
            end
            mul3 = b3 - tmp + 1;
            b3 = b3 + 1;
            sum_mul3 = sum_mul3 + (mul3 - 1);
            
            fprintf(fid, '\t\tcurrent_entity_set[%d] = [', cnt);
            for i1 = 1:p1+1
                for j1 = 1:p2+1
                    for k1 = 1:p3+1
                        id1 = i + i1 - 1 + sum_mul1;
                        id2 = j + j1 - 1 + sum_mul2;
                        id3 = k + k1 - 1 + sum_mul3;
                        id = id1 + (id2 + (id3 - 1) * v_dim - 1) * u_dim;
                        fprintf(fid, '%d', id);
                        if (i1==p1+1) && (j1==p2+1) && (k1==p3+1)
                            fprintf(fid, ']\n');
                        else
                            fprintf(fid, ', ');
                        end
                    end
                end
            end
            cnt = cnt + 1;
        end
    end
end
fprintf(fid, '\t\tself.layer_entities_sets[''%s''] = current_entity_set\n\n', name);
fprintf('number of printed entities = %d\n', cnt - 1);

%%
fprintf(fid, '\t\t# entity data\n');
fprintf(fid, '\t\tcurrent_entity_info = {}\n');
fprintf(fid, '\t\ttemp = {}\n');
cnt = 1;
b1 = p1 + 1;
sum_mul1 = 0;
for i = 1:ne1
    
    %check the multiplicity
    tmp = b1;
    while b1 <= (n1 + p1 + 1) && nurbs.knots{1}(b1+1) == nurbs.knots{1}(b1)
        b1 = b1 + 1;
    end
    mul1 = b1 - tmp + 1;
    b1 = b1 + 1;
    sum_mul1 = sum_mul1 + (mul1 - 1);
    
    b2 = p2 + 1;
    sum_mul2 = 0;
    for j = 1:ne2
        
        %check the multiplicity
        tmp = b2;
        while b2 <= (n2 + p2 + 1) && nurbs.knots{2}(b2+1) == nurbs.knots{2}(b2)
            b2 = b2 + 1;
        end
        mul2 = b2 - tmp + 1;
        b2 = b2 + 1;
        sum_mul2 = sum_mul2 + (mul2 - 1);
        
        b3 = p3 + 1;
        sum_mul3 = 0;
        for k = 1:ne3
            
            %check the multiplicity
            tmp = b3;
            while b3 <= (n3 + p3 + 1) && nurbs.knots{3}(b3+1) == nurbs.knots{3}(b3)
                b3 = b3 + 1;
            end
            mul3 = b3 - tmp + 1;
            b3 = b3 + 1;
            sum_mul3 = sum_mul3 + (mul3 - 1);
            
            fprintf(fid, '\t\ttemp[%d] = [', cnt);
            for i1 = 1:p1+1
                for j1 = 1:p2+1
                    for k1 = 1:p3+1
                        id1 = i + i1 - 1 + sum_mul1;
                        id2 = j + j1 - 1 + sum_mul2;
                        id3 = k + k1 - 1 + sum_mul3;
                        point = nurbs.coefs(:, id1, id2, id3);
                        fprintf(fid, '%.10f', point(4));
                        if (i1==p1+1) && (j1==p2+1) && (k1==p3+1)
                            fprintf(fid, ']\n');
                        else
                            fprintf(fid, ', ');
                        end
                    end
                end
            end
            cnt = cnt + 1;
        end
    end
end
fprintf(fid, '\t\tcurrent_entity_info[''NURBS_WEIGHTS''] = temp\n\n');

fprintf(fid, '\t\ttemp = {}\n');
cnt = 1;
for k1 = 1:ne1
    for k2 = 1:ne2
        for k3 = 1:ne3
            %Ce = C(:, :, cnt);
            Ce = mat2mcsr(C(:, :, cnt)); %change to modified compressed sparse row matrix
            sc1 = size(Ce, 1);
            sc2 = size(Ce, 2);
            fprintf(fid, '\t\ttemp[%d] = [', cnt);
            for i = 1:sc1
                fprintf(fid, '[');
                for j = 1:sc2
                    if i == 1 %change to modified compressed sparse row matrix
                        fprintf(fid, ' %d', Ce(i, j));
                    else
                        fprintf(fid, ' %.10f', Ce(i, j));
                    end
                    if (j==sc2)
                        fprintf(fid, ']');
                    else
                        fprintf(fid, ',');
                    end
                end
                if (i==sc1)
                    fprintf(fid, ']\n');
                else
                    fprintf(fid, ',');
                end
            end
            cnt = cnt + 1;
        end
    end
end
fprintf(fid, '\t\tcurrent_entity_info[''EXTRACTION_OPERATOR_MCSR''] = temp\n\n');

fprintf(fid, '\t\ttemp = {}\n');
cnt = 1;
for i = 1:ne1
    for j = 1:ne2
        for k = 1:ne3
            fprintf(fid, '\t\ttemp[%d] = %d\n', cnt, nurbs.order(1) - 1);
            cnt = cnt + 1;
        end
    end
end
fprintf(fid, '\t\tcurrent_entity_info[''NURBS_DEGREE_1''] = temp\n\n');

fprintf(fid, '\t\ttemp = {}\n');
cnt = 1;
for i = 1:ne1
    for j = 1:ne2
        for k = 1:ne3
            fprintf(fid, '\t\ttemp[%d] = %d\n', cnt, nurbs.order(2) - 1);
            cnt = cnt + 1;
        end
    end
end
fprintf(fid, '\t\tcurrent_entity_info[''NURBS_DEGREE_2''] = temp\n\n');

fprintf(fid, '\t\ttemp = {}\n');
cnt = 1;
for i = 1:ne1
    for j = 1:ne2
        for k = 1:ne3
            fprintf(fid, '\t\ttemp[%d] = %d\n', cnt, nurbs.order(3) - 1);
            cnt = cnt + 1;
        end
    end
end
fprintf(fid, '\t\tcurrent_entity_info[''NURBS_DEGREE_3''] = temp\n\n');

fprintf(fid, '\t\ttemp = {}\n');
cnt = 1;
div_u = params.num_division_u;
for i = 1:ne1
    for j = 1:ne2
        for k = 1:ne3
            fprintf(fid, '\t\ttemp[%d] = %d\n', cnt, div_u);
            cnt = cnt + 1;
        end
    end
end
fprintf(fid, '\t\tcurrent_entity_info[''NUM_DIVISION_1''] = temp\n\n');

fprintf(fid, '\t\ttemp = {}\n');
cnt = 1;
div_v = params.num_division_v;
for i = 1:ne1
    for j = 1:ne2
        for k = 1:ne3
            fprintf(fid, '\t\ttemp[%d] = %d\n', cnt, div_v);
            cnt = cnt + 1;
        end
    end
end
fprintf(fid, '\t\tcurrent_entity_info[''NUM_DIVISION_2''] = temp\n\n');

fprintf(fid, '\t\ttemp = {}\n');
cnt = 1;
div_w = params.num_division_w;
for i = 1:ne1
    for j = 1:ne2
        for k = 1:ne3
            fprintf(fid, '\t\ttemp[%d] = %d\n', cnt, div_w);
            cnt = cnt + 1;
        end
    end
end
fprintf(fid, '\t\tcurrent_entity_info[''NUM_DIVISION_3''] = temp\n\n');

if isfield(params,'activation_level')
    fprintf(fid, '\t\ttemp = {}\n');
    cnt = 1;
    alevel = params.activation_level;
    for i = 1:ne1
        for j = 1:ne2
            for k = 1:ne3
                fprintf(fid, '\t\ttemp[%d] = %d\n', cnt, alevel);
                cnt = cnt + 1;
            end
        end
    end
    fprintf(fid, '\t\tcurrent_entity_info[''ACTIVATION_LEVEL''] = temp\n\n');
end

fprintf(fid, '\t\tself.layer_entity_info_sets[''%s''] = current_entity_info\n\n', name);

%%
fprintf(fid, '\t\t# layer boundary marker\n');
fprintf(fid, '\t\tboundary_marker = {}\n');

% extract the boundary control point indices
boundary_u0 = zeros(v_dim,w_dim);
boundary_u1 = zeros(v_dim,w_dim);
boundary_v0 = zeros(u_dim,w_dim);
boundary_v1 = zeros(u_dim,w_dim);
boundary_w0 = zeros(u_dim,v_dim);
boundary_w1 = zeros(u_dim,v_dim);
boundary_u0v0 = zeros(1,w_dim);
boundary_u0v1 = zeros(1,w_dim);
boundary_u1v0 = zeros(1,w_dim);
boundary_u1v1 = zeros(1,w_dim);
boundary_u0w0 = zeros(1,v_dim);
boundary_u0w1 = zeros(1,v_dim);
boundary_u1w0 = zeros(1,v_dim);
boundary_u1w1 = zeros(1,v_dim);
boundary_v0w0 = zeros(1,u_dim);
boundary_v0w1 = zeros(1,u_dim);
boundary_v1w0 = zeros(1,u_dim);
boundary_v1w1 = zeros(1,u_dim);
boundary_u0v0w0 = zeros(1,1);
boundary_u1v0w0 = zeros(1,1);
boundary_u0v1w0 = zeros(1,1);
boundary_u1v1w0 = zeros(1,1);
boundary_u0v0w1 = zeros(1,1);
boundary_u1v0w1 = zeros(1,1);
boundary_u0v1w1 = zeros(1,1);
boundary_u1v1w1 = zeros(1,1);
cnt = 1;
for i = 1:w_dim
    for j = 1:v_dim
        for k = 1:u_dim
%             point = nurbs.coefs(:, k, j, i);
            if k==1
                boundary_u0(j,i) = cnt;
            elseif k==u_dim
                boundary_u1(j,i) = cnt;
            end

            if j==1
                boundary_v0(k,i) = cnt;
            elseif j==v_dim
                boundary_v1(k,i) = cnt;
            end

            if i==1
                boundary_w0(k,j) = cnt;
            elseif i==w_dim
                boundary_w1(k,j) = cnt;
            end

            if k==1 && j==1
                boundary_u0v0(i) = cnt;
            elseif k==1 && j==v_dim
                boundary_u0v1(i) = cnt;
            elseif k==u_dim && j==1
                boundary_u1v0(i) = cnt;
            elseif k==u_dim && j==v_dim
                boundary_u1v1(i) = cnt;
            end

            if k==1 && i==1
                boundary_u0w0(j) = cnt;
            elseif k==1 && i==w_dim
                boundary_u0w1(j) = cnt;
            elseif k==u_dim && i==1
                boundary_u1w0(j) = cnt;
            elseif k==u_dim && i==w_dim
                boundary_u1w1(j) = cnt;
            end

            if j==1 && i==1
                boundary_v0w0(k) = cnt;
            elseif j==1 && i==w_dim
                boundary_v0w1(k) = cnt;
            elseif j==v_dim && i==1
                boundary_v1w0(k) = cnt;
            elseif j==v_dim && i==w_dim
                boundary_v1w1(k) = cnt;
            end

            if k==1 && j==1 && i==1
                boundary_u0v0w0 = cnt;
            elseif k==1 && j==1 && i==w_dim
                boundary_u0v0w1 = cnt;
            elseif k==1 && j==v_dim && i==1
                boundary_u0v1w0 = cnt;
            elseif k==1 && j==v_dim && i==w_dim
                boundary_u0v1w1 = cnt;
            elseif k==u_dim && j==1 && i==1
                boundary_u1v0w0 = cnt;
            elseif k==u_dim && j==1 && i==w_dim
                boundary_u1v0w1 = cnt;
            elseif k==u_dim && j==v_dim && i==1
                boundary_u1v1w0 = cnt;
            elseif k==u_dim && j==v_dim && i==w_dim
                boundary_u1v1w1 = cnt;
            end

            cnt = cnt + 1;
        end
    end
end

fprintf(fid, '\t\tboundary_marker[''u0''] = ');
print_integer_matrix(fid,boundary_u0);
fprintf(fid, '\t\tboundary_marker[''u1''] = ');
print_integer_matrix(fid,boundary_u1);
fprintf(fid, '\t\tboundary_marker[''v0''] = ');
print_integer_matrix(fid,boundary_v0);
fprintf(fid, '\t\tboundary_marker[''v1''] = ');
print_integer_matrix(fid,boundary_v1);
fprintf(fid, '\t\tboundary_marker[''w0''] = ');
print_integer_matrix(fid,boundary_w0);
fprintf(fid, '\t\tboundary_marker[''w1''] = ');
print_integer_matrix(fid,boundary_w1);

fprintf(fid, '\t\tboundary_marker[''u0v0''] = ');
print_integer_vector(fid,boundary_u0v0);
fprintf(fid, '\t\tboundary_marker[''u0v1''] = ');
print_integer_vector(fid,boundary_u0v1);
fprintf(fid, '\t\tboundary_marker[''u1v0''] = ');
print_integer_vector(fid,boundary_u1v0);
fprintf(fid, '\t\tboundary_marker[''u1v1''] = ');
print_integer_vector(fid,boundary_u1v1);

fprintf(fid, '\t\tboundary_marker[''u0w0''] = ');
print_integer_vector(fid,boundary_u0w0);
fprintf(fid, '\t\tboundary_marker[''u0w1''] = ');
print_integer_vector(fid,boundary_u0w1);
fprintf(fid, '\t\tboundary_marker[''u1w0''] = ');
print_integer_vector(fid,boundary_u1w0);
fprintf(fid, '\t\tboundary_marker[''u1w1''] = ');
print_integer_vector(fid,boundary_u1w1);

fprintf(fid, '\t\tboundary_marker[''v0w0''] = ');
print_integer_vector(fid,boundary_v0w0);
fprintf(fid, '\t\tboundary_marker[''v0w1''] = ');
print_integer_vector(fid,boundary_v0w1);
fprintf(fid, '\t\tboundary_marker[''v1w0''] = ');
print_integer_vector(fid,boundary_v1w0);
fprintf(fid, '\t\tboundary_marker[''v1w1''] = ');
print_integer_vector(fid,boundary_v1w1);

fprintf(fid, '\t\tboundary_marker[''u0v0w0''] = ');
print_integer_vector(fid,boundary_u0v0w0);
fprintf(fid, '\t\tboundary_marker[''u1v0w0''] = ');
print_integer_vector(fid,boundary_u1v0w0);
fprintf(fid, '\t\tboundary_marker[''u0v1w0''] = ');
print_integer_vector(fid,boundary_u0v1w0);
fprintf(fid, '\t\tboundary_marker[''u1v1w0''] = ');
print_integer_vector(fid,boundary_u1v1w0);
fprintf(fid, '\t\tboundary_marker[''u0v0w1''] = ');
print_integer_vector(fid,boundary_u0v0w1);
fprintf(fid, '\t\tboundary_marker[''u1v0w1''] = ');
print_integer_vector(fid,boundary_u1v0w1);
fprintf(fid, '\t\tboundary_marker[''u0v1w1''] = ');
print_integer_vector(fid,boundary_u0v1w1);
fprintf(fid, '\t\tboundary_marker[''u1v1w1''] = ');
print_integer_vector(fid,boundary_u1v1w1);

fprintf(fid, '\t\tself.layer_boundary_marker[''%s''] = boundary_marker\n\n', name);

%%
fprintf(fid, '\t\t# layer attributes\n');
fprintf(fid, '\t\tself.layer_attributes[''%s''] = {}\n', name);

fprintf(fid, '\t\t################### end layer_info %s ###################\n\n', name);
fprintf('write layer %s completed\n', name);
disp('------------------------');

end

% function print_integer_matrix(fid,M)
%     fprintf(fid,'[');
%     for i = 1:size(M,1)
%         fprintf(fid,'[');
%         for j = 1:size(M,2)
%             fprintf(fid,' %d',M(i,j));
%             if j < size(M,2)
%                 fprintf(fid,',');
%             end
%         end
%         if i < size(M,1)
%             fprintf(fid,'], ');
%         else
%             fprintf(fid,']');
%         end
%     end
%     fprintf(fid, ']\n');
% end

