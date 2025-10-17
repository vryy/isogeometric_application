function write_layer_bezier_2d(fid, name, nurbs, C, Cxi, Cet, params)
%
% write the Bezier geometry info for 1 layer for 1 nurbs surface patch
%
% Input:
%   fid     file handler
%   name    layer name
%   nurbs   nurbs geometry
%   C       Bezier extraction operator
%   Cxi     Bezier extraction operator - xi
%   Cet     Bezier extraction operator - eta
%   params  additional parameters

%%
sizes = size(nurbs.coefs);
u_dim = sizes(2);
v_dim = sizes(3);

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

ne1 = size(Cxi, 3); %number of element along xi direction
ne2 = size(Cet, 3); %number of element along eta direction
ne = size(C, 3);

if ne1 * ne2 ~= ne
    error('invalid number of elements along 2 directions')
end

fprintf(fid, '\t\t## begin layer_info for layer %s\n', name);
fprintf(fid, '\t\tself.layer_list.append(''%s'')\n\n', name);

%%
fprintf(fid, '\t\t# nodal info\n');
fprintf(fid, '\t\tcurrent_nodal_set = {}\n');
cnt = 1;
for i = 1:v_dim
    for j = 1:u_dim
        point = nurbs.coefs(:, j, i);
        point(1:3) = point(1:3) / point(4);
        fprintf(fid, '\t\tcurrent_nodal_set[%d] = [%.10f, %.10f, %.10f]\n', cnt, point(1), point(2), point(3));
        cnt = cnt + 1;
    end
end
fprintf(fid, '\t\tself.layer_nodes_sets[''%s''] = current_nodal_set\n\n', name);
fprintf('number of printed nodes = %d\n', cnt - 1);

%%
fprintf(fid, '\t\t# entity connectivities\n');
fprintf(fid, '\t\tcurrent_entity_set = {}\n');
cnt = 1;
prop = 1;
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

        fprintf(fid, '\t\tcurrent_entity_set[%d] = [', cnt);
        for i1 = 1:p1+1
            for j1 = 1:p2+1
                id1 = i + i1 - 1 + sum_mul1;
                id2 = j + j1 - 1 + sum_mul2;
                id = id1 + (id2 - 1) * u_dim;
                fprintf(fid, '%d', id);
                if (i1==p1+1) && (j1==p2+1)
                    fprintf(fid, ']\n');
                else
                    fprintf(fid, ', ');
                end
            end
        end
        cnt = cnt + 1;
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

        fprintf(fid, '\t\ttemp[%d] = [', cnt);
        for i1 = 1:p1+1
            for j1 = 1:p2+1
               id1 = i + i1 - 1 + sum_mul1;
               id2 = j + j1 - 1 + sum_mul2;
                point = nurbs.coefs(:, id1, id2);
                fprintf(fid, '%.10f', point(4));
                if (i1==p1+1) && (j1==p2+1)
                    fprintf(fid, ']\n');
                else
                    fprintf(fid, ', ');
                end
            end
        end
        cnt = cnt + 1;
    end
end
fprintf(fid, '\t\tcurrent_entity_info[''NURBS_WEIGHTS''] = temp\n\n');

fprintf(fid, '\t\ttemp = {}\n');
cnt = 1;
for k1 = 1:ne1
    for k2 = 1:ne2
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
fprintf(fid, '\t\tcurrent_entity_info[''EXTRACTION_OPERATOR_MCSR''] = temp\n\n');

fprintf(fid, '\t\ttemp = {}\n');
cnt = 1;
for i = 1:ne1
    for j = 1:ne2
        fprintf(fid, '\t\ttemp[%d] = %d\n', cnt, nurbs.order(1) - 1);
        cnt = cnt + 1;
    end
end
fprintf(fid, '\t\tcurrent_entity_info[''NURBS_DEGREE_1''] = temp\n\n');

fprintf(fid, '\t\ttemp = {}\n');
cnt = 1;
for i = 1:ne1
    for j = 1:ne2
        fprintf(fid, '\t\ttemp[%d] = %d\n', cnt, nurbs.order(2) - 1);
        cnt = cnt + 1;
    end
end
fprintf(fid, '\t\tcurrent_entity_info[''NURBS_DEGREE_2''] = temp\n\n');

fprintf(fid, '\t\ttemp = {}\n');
cnt = 1;
div_u = params.num_division_u;
for i = 1:ne1
    for j = 1:ne2
        fprintf(fid, '\t\ttemp[%d] = %d\n', cnt, div_u);
        cnt = cnt + 1;
    end
end
fprintf(fid, '\t\tcurrent_entity_info[''NUM_DIVISION_1''] = temp\n\n');

fprintf(fid, '\t\ttemp = {}\n');
cnt = 1;
div_v = params.num_division_v;
for i = 1:ne1
    for j = 1:ne2
        fprintf(fid, '\t\ttemp[%d] = %d\n', cnt, div_v);
        cnt = cnt + 1;
    end
end
fprintf(fid, '\t\tcurrent_entity_info[''NUM_DIVISION_2''] = temp\n\n');

if isfield(params,'master_index')
    fprintf(fid, '\t\ttemp = {}\n');
    cnt = 1;
    mi = params.master_index;
    for i = 1:ne1
        for j = 1:ne2
            fprintf(fid, '\t\ttemp[%d] = %d\n', cnt, mi);
            cnt = cnt + 1;
        end
    end
    fprintf(fid, '\t\tcurrent_entity_info[''MASTER_INDEX''] = temp\n\n');
end

if isfield(params,'slave_index')
    fprintf(fid, '\t\ttemp = {}\n');
    cnt = 1;
    si = params.slave_index;
    for i = 1:ne1
        for j = 1:ne2
            fprintf(fid, '\t\ttemp[%d] = %d\n', cnt, si);
            cnt = cnt + 1;
        end
    end
    fprintf(fid, '\t\tcurrent_entity_info[''SLAVE_INDEX''] = temp\n\n');
end

fprintf(fid, '\t\tself.layer_entity_info_sets[''%s''] = current_entity_info\n\n', name);

%%
fprintf(fid, '\t\t# layer attributes\n');
fprintf(fid, '\t\tself.layer_attributes[''%s''] = {}\n', name);
fprintf(fid, '\t\tself.layer_attributes[''%s''][''%s''] = %d\n', name, 'NURBS_NUMBER_1', nurbs.number(1));
fprintf(fid, '\t\tself.layer_attributes[''%s''][''%s''] = %d\n', name, 'NURBS_NUMBER_2', nurbs.number(2));
fprintf(fid, '\t\tself.layer_attributes[''%s''][''%s''] = %d\n', name, 'NURBS_DEGREE_1', nurbs.order(1)-1);
fprintf(fid, '\t\tself.layer_attributes[''%s''][''%s''] = %d\n', name, 'NURBS_DEGREE_2', nurbs.order(2)-1);

fprintf(fid, '\t\t################### end layer_info %s ###################\n\n', name);
fprintf('write layer %s completed\n', name);
disp('------------------------');
