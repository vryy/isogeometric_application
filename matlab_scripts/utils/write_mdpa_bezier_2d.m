%generates mdpa bezier data for single nurbs patch with open knot vector,
%with knot multiplicity
function write_mdpa_bezier_2d(nurbs, C, Cxi, Cet, fn, params)

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

fid = fopen(fn, 'wt');

fprintf(fid, '//KRATOS isogeometric application data file\n');
fprintf(fid, '//(c) 2014 Hoang Giang Bui, Ruhr-University Bochum\n\n');

fprintf(fid, 'Begin ModelPartData\n');
fprintf(fid, 'End ModelPartData\n\n');

fprintf(fid, 'Begin Properties 1\n');
fprintf(fid, 'End Properties\n\n');

fprintf(fid, 'Begin Nodes\n');
cnt = 1;
for i = 1:v_dim
    for j = 1:u_dim
        point = nurbs.coefs(:, j, i);
        point(1:3) = point(1:3) / point(4);
        fprintf(fid, '%d %f %f %f\n', cnt, point(1), point(2), point(3));
        cnt = cnt + 1;
    end
end
fprintf(fid, 'End Nodes\n\n');

fprintf(fid, 'Begin Elements KinematicLinearGeo2dBezier\n');
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
        fprintf(fid, '%d %d', cnt, prop);
        
        %check the multiplicity
        tmp = b2;
        while b2 <= (n2 + p2 + 1) && nurbs.knots{2}(b2+1) == nurbs.knots{2}(b2)
            b2 = b2 + 1;
        end
        mul2 = b2 - tmp + 1;
        b2 = b2 + 1;
        sum_mul2 = sum_mul2 + (mul2 - 1);
        
        for k = 1:p1+1
            for l = 1:p2+1
                id1 = i + k - 1 + sum_mul1;
                id2 = j + l - 1 + sum_mul2;
                id = id1 + (id2 - 1) * u_dim;
                fprintf(fid, ' %d', id);
            end
        end
        cnt = cnt + 1;
        fprintf(fid, '\n');
    end
end
fprintf(fid, 'End Elements\n\n');

fprintf(fid, 'Begin ElementalData NURBS_WEIGHTS\n');
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
        
        fprintf(fid, '%d [%d] (', cnt, (p1 + 1) * (p2 + 1));
        for k = 1:p1+1
            for l = 1:p2+1
                id1 = i + k - 1 + sum_mul1;
                id2 = j + l - 1 + sum_mul2;
                point = nurbs.coefs(:, id1, id2);
                fprintf(fid, ' %f', point(4));
                if (k==p1+1) && (l==p2+1)
                    fprintf(fid, ')\n');
                else
                    fprintf(fid, ',');
                end
            end
        end
        cnt = cnt + 1;
    end
end
fprintf(fid, 'End ElementalData\n\n');

fprintf(fid, 'Begin ElementalData EXTRACTION_OPERATOR_MCSR\n');
cnt = 1;
for k = 1:ne1
    for l = 1:ne2
        %Ce = C(:, :, cnt);
        Ce = mat2mcsr(C(:, :, cnt)); %change to modified compressed sparse row matrix
        sc1 = size(Ce, 1);
        sc2 = size(Ce, 2);
        fprintf(fid, '%d [%d, %d] (', cnt, sc1, sc2);
        for i = 1:sc1
            fprintf(fid, '(');
            for j = 1:sc2
                if i == 1 %change to modified compressed sparse row matrix
                    fprintf(fid, ' %d', Ce(i, j));
                else
                    fprintf(fid, ' %f', Ce(i, j));
                end
                if (j==sc2)
                    fprintf(fid, ')');
                else
                    fprintf(fid, ',');
                end
            end
            if (i==sc1)
                fprintf(fid, ')\n');
            else
                fprintf(fid, ',');
            end
        end
        cnt = cnt + 1;
    end
end
fprintf(fid, 'End ElementalData\n\n');

fprintf(fid, 'Begin ElementalData NURBS_DEGREE_1\n');
cnt = 1;
for i = 1:ne1
    for j = 1:ne2
        fprintf(fid, '%d %d', cnt, nurbs.order(1) - 1);
        cnt = cnt + 1;
        fprintf(fid, '\n');
    end
end
fprintf(fid, 'End ElementalData\n\n');

fprintf(fid, 'Begin ElementalData NURBS_DEGREE_2\n');
cnt = 1;
for i = 1:ne1
    for j = 1:ne2
        fprintf(fid, '%d %d', cnt, nurbs.order(2) - 1);
        cnt = cnt + 1;
        fprintf(fid, '\n');
    end
end
fprintf(fid, 'End ElementalData\n\n');

fprintf(fid, 'Begin ElementalData NUM_DIVISION_1\n');
cnt = 1;
div_u = params.num_division_u;
for i = 1:ne1
    for j = 1:ne2
        fprintf(fid, '%d %d', cnt, div_u);
        cnt = cnt + 1;
        fprintf(fid, '\n');
    end
end
fprintf(fid, 'End ElementalData\n\n');

fprintf(fid, 'Begin ElementalData NUM_DIVISION_2\n');
cnt = 1;
div_v = params.num_division_v;
for i = 1:ne1
    for j = 1:ne2
        fprintf(fid, '%d %d', cnt, div_v);
        cnt = cnt + 1;
        fprintf(fid, '\n');
    end
end
fprintf(fid, 'End ElementalData\n\n');

fclose(fid);
