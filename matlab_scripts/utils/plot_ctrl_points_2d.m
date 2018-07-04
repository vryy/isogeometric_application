function plot_ctrl_points_2d(nurbs,params)
% Plot the control points of a single 2D patch
%   params.label: turn the id of each control point on/off
%   params.axis: turn the axis on/off
%   params.legend: turn the legend for u/v on/off
%   params.number: if it is given, the id of each control points will be
%   derived by the given value; otherwise, incremental values will be used.
axis equal
hold on

if nargin==1
    params.label='off';
    params.axis='off';
end

sizes = size(nurbs.coefs);
u_dim = sizes(2);
v_dim = sizes(3);

u_color = 'black';
v_color = 'cyan';

if ~isfield(params,'point_color')
    params.point_color = 'blue';
end

if ~isfield(params,'label')
    params.label = 'off';
end

if ~isfield(params,'axis')
    params.axis = 'off';
end

if ~isfield(params,'legend')
    params.legend = 'on';
end

if ~isfield(params,'number')
    params.number = 1:(u_dim*v_dim);
end

%%
cnt = 1;
for j = 1:v_dim
    for k = 1:u_dim
        point = nurbs.coefs(:, k, j);
        point(1:3) = point(1:3) / point(4);
        scatter3(point(1), point(2), point(3), params.point_color);
        if strcmp(params.label,'on')
            text(point(1), point(2), point(3), num2str(params.number(cnt)), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'cap');
        end
        cnt = cnt + 1;
        if k > 1
            if k == 2
                if strcmp(params.legend,'on')
                    L = quiver3(old_point_u(1), old_point_u(2), old_point_u(3), point(1)-old_point_u(1), point(2)-old_point_u(2), point(3)-old_point_u(3));
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


%%
for k = 1:u_dim
    for j = 1:v_dim
        point = nurbs.coefs(:, k, j);
        point(1:3) = point(1:3) / point(4);
        if j > 1
            if j == 2
                if strcmp(params.legend,'on')
                    L = quiver3(old_point_v(1), old_point_v(2), old_point_v(3), point(1)-old_point_v(1), point(2)-old_point_v(2), point(3)-old_point_v(3));
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

if strcmp(params.legend,'on')
    legend([u_plot_for_legend,v_plot_for_legend], 'u-dim', 'v-dim');
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

