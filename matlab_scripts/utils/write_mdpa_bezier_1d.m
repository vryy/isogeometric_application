%generates mdpa bezier data for single nurbs patch with open knot vector,
%no knot multiplicity
function write_mdpa_bezier_1d(nurbs, C)

sizes = size(nurbs.coefs);
u_dim = sizes(2);

%check if n >= p + 1
n = u_dim;
p = nurbs.order - 1;
if(n < p + 1)
    error('n < p + 1: not enough Bezier elements')
end

ne = size(C, 3); %number of Bezier element

fid = fopen('nurbs.mdpa', 'wt');

fprintf(fid, 'Begin Nodes\n');
cnt = 1;
weights = zeros(1, u_dim);
for j = 1:u_dim
    point = nurbs.coefs(:, j);
    weights(j) = point(4);
    point(1:3) = point(1:3) / point(4);
    fprintf(fid, '%d %f %f %f\n', cnt, point(1), point(2), point(3));
    cnt = cnt + 1;
end
fprintf(fid, 'End Nodes\n\n');

fprintf(fid, 'Begin Elements KinematicLinearGeo1dBezier\n');
cnt = 1;
prop = 1;
b = p + 1;
sum_mul = 0;
for i = 1:ne
    fprintf(fid, '%d %d', cnt, prop);
    
    %check the multiplicity
    tmp = b;
    while b <= (n + p + 1) && nurbs.knots(b+1) == nurbs.knots(b)
        b = b + 1;
    end
    mul = b - tmp + 1;
    b = b + 1;
    sum_mul = sum_mul + mul - 1;
    
    for j = 1:p+1
        id = i + j - 1 + sum_mul;
        fprintf(fid, ' %d', id);
    end
    cnt = cnt + 1;
    fprintf(fid, '\n');
end
fprintf(fid, 'End Elements\n\n');

fprintf(fid, 'Begin ElementalData NURBS_WEIGHTS\n');
cnt = 1;
b = p + 1;
sum_mul = 0;
for i = 1:ne
    fprintf(fid, '%d [%d] (', cnt, p + 1);
    
    %check the multiplicity
    tmp = b;
    while b <= (n + p + 1) && nurbs.knots(b+1) == nurbs.knots(b)
        b = b + 1;
    end
    mul = b - tmp + 1;
    b = b + 1;
    sum_mul = sum_mul + mul - 1;
    
    for j = 1:p+1
        id = i + j - 1 + sum_mul;
        point = nurbs.coefs(:, id);
        fprintf(fid, ' %f', point(4));
        if (j==p + 1)
            fprintf(fid, ')\n');
        else
            fprintf(fid, ',');
        end
    end
    cnt = cnt + 1;
end
fprintf(fid, 'End ElementalData\n\n');

fprintf(fid, 'Begin ElementalData NURBS_DEGREE_1\n');
cnt = 1;
for i = 1:ne
    fprintf(fid, '%d %d\n', cnt, p);
    cnt = cnt + 1;
end
fprintf(fid, 'End ElementalData\n\n');

fprintf(fid, 'Begin ElementalData EXTRACTION_OPERATOR\n');
cnt = 1;
for k = 1:ne
    %Ce = C(:, :, k);
    Ce = mat2mcsr(C(:, :, k)); %change to modified compressed sparse row matrix
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
fprintf(fid, 'End ElementalData\n\n');



fclose(fid);
