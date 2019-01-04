%generates mdpa bezier data for single tspline
function write_mdpa_bezier_tspline_2d(tspline, fn)

fid = fopen(fn, 'wt');

fprintf(fid, '//KRATOS isogeometric application data file\n');
fprintf(fid, '//(c) 2014 Hoang Giang Bui, Ruhr-University Bochum\n\n');

fprintf(fid, 'Begin ModelPartData\n');
fprintf(fid, 'End ModelPartData\n\n');

fprintf(fid, 'Begin Properties 1\n');
fprintf(fid, 'End Properties\n\n');

fprintf(fid, 'Begin Nodes\n');
cnt = 1;
for i = 1:tspline.ndof
	point = tspline.control_points(:, i);
    %point(1:3) = point(1:3) / point(4);
    fprintf(fid, '%d %f %f %f\n', cnt, point(1), point(2), point(3));
    cnt = cnt + 1;
end
fprintf(fid, 'End Nodes\n\n');

fprintf(fid, 'Begin Elements KinematicLinearGeo2dBezier\n');
cnt = 1;
prop = 1;
for i = 1:tspline.nel
    fprintf(fid, '%d %d', cnt, prop);
    elem = tspline.elements(i);
    con = elem.connectivity;
    for j = 1:elem.nsh
        fprintf(fid, ' %d', con(j));
    end
    cnt = cnt + 1;
    fprintf(fid, '\n');
end
fprintf(fid, 'End Elements\n\n');

fprintf(fid, 'Begin ElementalData NURBS_WEIGHTS\n');
cnt = 1;
for i = 1:tspline.nel
    elem = tspline.elements(i);
    con = elem.connectivity;
    fprintf(fid, '%d [%d] (', cnt, elem.nsh);
    points = tspline.control_points(:, con);
    for j = 1:elem.nsh
        fprintf(fid, ' %d', points(4, j));
        if j == elem.nsh
            fprintf(fid, ')\n');
        else
            fprintf(fid, ',');
        end
    end
    cnt = cnt + 1;
end
fprintf(fid, 'End ElementalData\n\n');

fprintf(fid, 'Begin ElementalData EXTRACTION_OPERATOR\n');
cnt = 1;
for k = 1:tspline.nel
    elem = tspline.elements(k);
    %Ce = elem.extraction;
    Ce = mat2mcsr(elem.extraction); %change to modified compressed sparse row matrix
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

fprintf(fid, 'Begin ElementalData NURBS_DEGREE_1\n');
cnt = 1;
for k = 1:tspline.nel
    elem = tspline.elements(k);
    p = elem.degree;
    fprintf(fid, '%d %d\n', cnt, p(1));
    cnt = cnt + 1;
end
fprintf(fid, 'End ElementalData\n\n');

fprintf(fid, 'Begin ElementalData NURBS_DEGREE_2\n');
cnt = 1;
for k = 1:tspline.nel
    elem = tspline.elements(k);
    p = elem.degree;
    fprintf(fid, '%d %d\n', cnt, p(2));
    cnt = cnt + 1;
end
fprintf(fid, 'End ElementalData\n\n');

fprintf(fid, 'Begin ElementalData NUM_DIVISION_1\n');
cnt = 1;
div_u = 1;
for k = 1:tspline.nel
    fprintf(fid, '%d %d\n', cnt, div_u);
    cnt = cnt + 1;
end
fprintf(fid, 'End ElementalData\n\n');

fprintf(fid, 'Begin ElementalData NUM_DIVISION_2\n');
cnt = 1;
div_v = 1;
for k = 1:tspline.nel
    fprintf(fid, '%d %d\n', cnt, div_v);
    cnt = cnt + 1;
end
fprintf(fid, 'End ElementalData\n\n');

fclose(fid);
