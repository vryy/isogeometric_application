%generate mdpa for single nurbs patch
function write_mdpa_2d(geometry)

nurbs = geometry.nurbs;

nurbs

fid = fopen('nurbs.mdpa', 'wt');

sizes = size(nurbs.coefs);
u_dim = sizes(2);
v_dim = sizes(3);

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

fprintf(fid, 'Begin Elements KinematicLinearGeo2dNURBS\n');
cnt = 1;
prop = 1;
fprintf(fid, '%d %d', cnt, prop);
for i = 1:u_dim
    for j = 1:v_dim
        id = i + (j-1) * u_dim;
        fprintf(fid, ' %d', id);
    end
end
fprintf(fid, '\nEnd Elements\n\n');

fprintf(fid, 'Begin ElementalData NURBS_WEIGHTS\n');
cnt = 1;
fprintf(fid, '%d [%d] (', cnt, u_dim * v_dim);
for i = 1:v_dim
    for j = 1:u_dim
        point = nurbs.coefs(:, j, i);
        fprintf(fid, ' %f', point(4));
        if (i==v_dim) && (j==u_dim)
            fprintf(fid, ')\n');
        else
            fprintf(fid, ',');
        end
    end
end
fprintf(fid, 'End ElementalData\n\n');

fprintf(fid, 'Begin ElementalData NURBS_KNOTS_1\n');
cnt = 1;
u_knots = nurbs.knots{1};
u_len = size(u_knots, 2);
fprintf(fid, '%d [%d] (', cnt, u_len);
for i = 1:u_len
    fprintf(fid, ' %f', u_knots(i));
    if (i==u_len)
        fprintf(fid, ')\n');
    else
        fprintf(fid, ',');
    end
end
fprintf(fid, 'End ElementalData\n\n');

fprintf(fid, 'Begin ElementalData NURBS_KNOTS_2\n');
cnt = 1;
v_knots = nurbs.knots{2};
v_len = size(v_knots, 2);
fprintf(fid, '%d [%d] (', cnt, v_len);
for i = 1:v_len
    fprintf(fid, ' %f', v_knots(i));
    if (i==v_len)
        fprintf(fid, ')\n');
    else
        fprintf(fid, ',');
    end
end
fprintf(fid, 'End ElementalData\n\n');

fprintf(fid, 'Begin ElementalData NURBS_DIMENSION_1\n');
cnt = 1;
fprintf(fid, '%d %d\n', cnt, u_dim);
fprintf(fid, 'End ElementalData\n\n');

fprintf(fid, 'Begin ElementalData NURBS_DIMENSION_2\n');
cnt = 1;
fprintf(fid, '%d %d\n', cnt, v_dim);
fprintf(fid, 'End ElementalData\n\n');

fprintf(fid, 'Begin ElementalData NURBS_DEGREE_1\n');
cnt = 1;
order = nurbs.order;
fprintf(fid, '%d %d\n', cnt, order(1) - 1);
fprintf(fid, 'End ElementalData\n\n');

fprintf(fid, 'Begin ElementalData NURBS_DEGREE_2\n');
cnt = 1;
order = nurbs.order;
fprintf(fid, '%d %d\n', cnt, order(2) - 1);
fprintf(fid, 'End ElementalData\n\n');

fprintf(fid, 'Begin ElementalData NUM_DIVISION_1\n');
cnt = 1;
div_u = 30;
fprintf(fid, '%d %d\n', cnt, div_u);
fprintf(fid, 'End ElementalData\n\n');

fprintf(fid, 'Begin ElementalData NUM_DIVISION_2\n');
cnt = 1;
div_v = 30;
fprintf(fid, '%d %d\n', cnt, div_v);
fprintf(fid, 'End ElementalData\n\n');

fclose(fid);
