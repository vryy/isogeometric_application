%% plot the control points of the hierarchical B-Splines
function plot_ctrl_points_hbsplines_2d(P,W,EqId,params)

if ~isfield(params, 'color')
    params.color = 'red';
end

axis off

Px = P(:,1) ./ W';
Py = P(:,2) ./ W';
plot(Px,Py,'o','color',params.color);

for i = 1:length(W)
    text(Px(i),Py(i),num2str(EqId(i)), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'cap');
end

