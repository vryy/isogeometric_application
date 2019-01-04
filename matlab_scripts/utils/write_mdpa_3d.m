function write_mdpa_3d(geometry)

nurbs = geometry.nurbs;

nurbs

fid = fopen('nurbs.mdpa', 'wt');

sizes = size(nurbs.coefs);
u_dim = sizes(2);
v_dim = sizes(3);
w_dim = sizes(4);

fprintf(fid, 'Begin Nodes\n');
cnt = 1;
for i = 1:w_dim
    for j = 1:v_dim
        for k = 1:u_dim
            point = nurbs.coefs(:, k, j, i);
            point(1:3) = point(1:3) / point(4);
            fprintf(fid, '%d %f %f %f\n', cnt, point(1), point(2), point(3));
            cnt = cnt + 1;
        end
    end
end
fprintf(fid, 'End Nodes\n\n');

fprintf(fid, 'Begin Elements KinematicLinearGeo3dNURBS\n');
cnt = 1;
prop = 1;
fprintf(fid, '%d %d', cnt, prop);
for i = 1:u_dim
    for j = 1:v_dim
        for k = 1:w_dim
            id = i + (j + (k - 1) * v_dim - 1) * u_dim;
            fprintf(fid, ' %d', id);
        end
    end
end
fprintf(fid, '\nEnd Elements\n\n');

fprintf(fid, 'Begin ElementalData NURBS_WEIGHTS\n');
cnt = 1;
fprintf(fid, '%d [%d] (', cnt, u_dim * v_dim * w_dim);
for i = 1:w_dim
    for j = 1:v_dim
        for k = 1:u_dim
            point = nurbs.coefs(:, k, j, i);
            fprintf(fid, ' %f', point(4));
            if (i==w_dim) && (j==v_dim) && (k==u_dim)
                fprintf(fid, ')\n');
            else
                fprintf(fid, ',');
            end
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

fprintf(fid, 'Begin ElementalData NURBS_KNOTS_3\n');
cnt = 1;
w_knots = nurbs.knots{3};
w_len = size(w_knots, 2);
fprintf(fid, '%d [%d] (', cnt, w_len);
for i = 1:w_len
    fprintf(fid, ' %f', w_knots(i));
    if (i==w_len)
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

fprintf(fid, 'Begin ElementalData NURBS_DIMENSION_3\n');
cnt = 1;
fprintf(fid, '%d %d\n', cnt, w_dim);
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

fprintf(fid, 'Begin ElementalData NURBS_DEGREE_3\n');
cnt = 1;
order = nurbs.order;
fprintf(fid, '%d %d\n', cnt, order(3) - 1);
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

fprintf(fid, 'Begin ElementalData NUM_DIVISION_3\n');
cnt = 1;
div_w = 30;
fprintf(fid, '%d %d\n', cnt, div_w);
fprintf(fid, 'End ElementalData\n\n');

fclose(fid);
