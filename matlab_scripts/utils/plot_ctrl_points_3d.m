function plot_ctrl_points_3d(nurbs,params)
% Plot the control points of a single 2D patch
%   params.label: turn the id of each control point on/off
%   params.axis: turn the axis on/off
%   params.legend: turn the legend for u/v on/off
%   params.number: if it is given, the id of each control points will be
%     derived by the given value; otherwise, incremental values will be used.
%   params.patch_id: if this is given, the patch id will be plotted at the
%     center of gravity of all control points
% Sometimes quiver3 does not show the arrow head correctly, due to some arrow
% head drawing face is the same as screen plane. This problem is reported here:
%   https://savannah.gnu.org/bugs/?47185
% Also, quiver3 will not display arrowheads for vectors which exist only
% in the z axis (ex. <0, 0, 1>) despite set(vector, 'ShowArrowHead', 'on'), see
%   https://savannah.gnu.org/bugs/?53187
axis equal;
hold on;

if nargin==1
    params.label='off';
    params.axis='off';
end

sizes = size(nurbs.coefs);
u_dim = sizes(2);
v_dim = sizes(3);
w_dim = sizes(4);

u_color = 'black';
v_color = 'cyan';
w_color = 'red';

if ~isfield(params,'point_color')
    params.point_color = 'blue';
end
rgbColor = colorNameToRGB(params.point_color); % Convert to RGB

if ~isfield(params,'label')
    params.label = 'off';
end

if ~isfield(params,'legend')
    params.legend = 'on';
end

if ~isfield(params,'axis')
    params.axis = 'off';
end

if ~isfield(params,'number')
    params.number = 1:(u_dim*v_dim*w_dim);
end

if ~isfield(params,'text_dc')
    params.text_dc = 1.01;
end

if ~isfield(params,'arrow_size')
    params.arrow_size = 0.33;
end

if ~isfield(params,'font_size')
  params.font_size = 20;
end

if ~isfield(params,'point_style')
  params.point_style = 'o';
end

%%
cnt = 1;
for i = 1:w_dim
    for j = 1:v_dim
        for k = 1:u_dim
            point = nurbs.coefs(:, k, j, i);
            point(1:3) = point(1:3) / point(4);
            S = scatter3(point(1),point(2),point(3));
            set(S,'Marker',params.point_style);
            set(S,'CData',rgbColor);
            if strcmp(params.label,'on')
                text(point(1)*params.text_dc, point(2)*params.text_dc, point(3)*params.text_dc, num2str(params.number(cnt)), "fontsize", params.font_size);
            end
            cnt = cnt + 1;
            if k > 1
                if k == 2
                    if strcmp(params.legend,'on')
                        L = quiver3(old_point_u(1),old_point_u(2),old_point_u(3),point(1)-old_point_u(1),point(2)-old_point_u(2),point(3)-old_point_u(3));
                        %set(L, 'markersize', 10.0);
                        set(L, "maxheadsize", params.arrow_size, 'AutoScale','on');
                        u_plot_for_legend = L;
                    else
                        L = line([old_point_u(1) point(1)], [old_point_u(2) point(2)], [old_point_u(3) point(3)]);
                    end
                else
                    L = line([old_point_u(1) point(1)], [old_point_u(2) point(2)], [old_point_u(3) point(3)]);
                end
                set(L, 'color', u_color);
            end
            old_point_u = point;
        end
    end
end

%%
for i = 1:w_dim
    for k = 1:u_dim
        for j = 1:v_dim
            point = nurbs.coefs(:, k, j, i);
            point(1:3) = point(1:3) / point(4);
            if j > 1
                if j == 2
                    if strcmp(params.legend,'on')
                        L = quiver3(old_point_v(1),old_point_v(2),old_point_v(3),point(1)-old_point_v(1),point(2)-old_point_v(2),point(3)-old_point_v(3));
                        set(L, "maxheadsize", params.arrow_size, 'AutoScale','on');
                        v_plot_for_legend = L;
                    else
                        L = line([old_point_v(1) point(1)], [old_point_v(2) point(2)], [old_point_v(3) point(3)]);
                    end
                else
                    L = line([old_point_v(1) point(1)], [old_point_v(2) point(2)], [old_point_v(3) point(3)]);
                end
                set(L, 'color', v_color);
            end
            old_point_v = point;
        end
    end
end

%%
for j = 1:v_dim
    for k = 1:u_dim
        for i = 1:w_dim
            point = nurbs.coefs(:, k, j, i);
            point(1:3) = point(1:3) / point(4);
            if i > 1
                if i == 2
                    if strcmp(params.legend,'on')
                        L = quiver3(old_point_w(1),old_point_w(2),old_point_w(3),point(1)-old_point_w(1),point(2)-old_point_w(2),point(3)-old_point_w(3));
                        set (L, "maxheadsize", params.arrow_size);
                        w_plot_for_legend = L;
                    else
                        L = line([old_point_w(1) point(1)], [old_point_w(2) point(2)], [old_point_w(3) point(3)]);
                    end
                else
                    L = line([old_point_w(1) point(1)], [old_point_w(2) point(2)], [old_point_w(3) point(3)]);
                end
                set(L, 'color', w_color);
            end
            old_point_w = point;
        end
    end
end

if strcmp(params.legend,'on')
    legend([u_plot_for_legend,v_plot_for_legend,w_plot_for_legend], 'u-dim', 'v-dim', 'w-dim');
end

if isfield(params,'patch_id')
    cnt = 1;
    cen = [0 0 0];
    for i = 1:w_dim
        for j = 1:v_dim
            for k = 1:u_dim
                point = nurbs.coefs(:, k, j, i);
                point(1:3) = point(1:3) / point(4);
                cen = cen + point;
                cnt = cnt + 1;
            end
        end
    end
    cen = cen / (w_dim*v_dim*u_dim);
    text(cen(1), cen(2), cen(3), num2str(params.patch_id));
%     [x y z] = sphere(10);
%     h = surfl(x+cen(1), y+cen(2), z+cen(3));
%     set(h, 'FaceAlpha', 0.5)
%     shading interp
end

xlabel('x');
ylabel('y');
zlabel('z');

if strcmp(params.axis,'on')
    Lx = line([0 1], [0 0], [0 0]);
    text(1, 0, 0, 'X');
    set(Lx, 'color', 'magenta');

    Ly = line([0 0], [0 1], [0 0]);
    text(0, 1, 0, 'Y');
    set(Ly, 'color', 'magenta');

    Lz = line([0 0], [0 0], [0 1]);
    text(0, 0, 1, 'Z');
    set(Lz, 'color', 'magenta');
end

