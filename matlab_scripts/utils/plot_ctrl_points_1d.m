function plot_ctrl_points_1d(nurbs,params)

axis equal;
hold on;

sizes = size(nurbs.coefs);
u_dim = sizes(2);

if nargin==1
    params.label='off';
    params.axis='off';
end

if ~isfield(params,'u_color')
    params.u_color = 'black';
end

% params.label
% params.axis

if ~isfield(params,'point_color')
    params.point_color = 'blue';
end
rgbColor = colorNameToRGB(params.point_color); % Convert to RGB

if ~isfield(params,'line_style')
    params.line_style = '-';
end

if ~isfield(params,'point_style')
    params.point_style = 'o';
end

if ~isfield(params,'axis')
    params.axis = 'off';
end

if ~isfield(params,'label')
    params.label = 'off';
end

%%
cnt = 1;
        for k = 1:u_dim
            point = nurbs.coefs(:, k);
            point(1:3) = point(1:3) / point(4);
            S = scatter3(point(1),point(2),point(3));
            set(S,'Marker',params.point_style);
            set(S,'CData',rgbColor);
            if strcmp(params.label,'on')
                text(point(1),point(2),point(3),num2str(cnt));
            end
            cnt = cnt + 1;
            if k > 1
                L = line([old_point_u(1) point(1)],[old_point_u(2) point(2)],[old_point_u(3) point(3)],'LineStyle',params.line_style);
                set(L,'color',params.u_color);
            end
            old_point_u = point;
        end

% u = plot(0,0,u_color);
% legend([u], 'u-dim');

if strcmp(params.axis,'on')
    Lx = line([0 1],[0 0],[0 0]);
    text(1,0,0,'X');
    set(Lx,'color','magenta');

    Ly = line([0 0],[0 1],[0 0]);
    text(0,1,0,'Y');
    set(Ly,'color','magenta');

    Lz = line([0 0],[0 0],[0 1]);
    text(0,0,1,'Z');
    set(Lz,'color','magenta');
end
