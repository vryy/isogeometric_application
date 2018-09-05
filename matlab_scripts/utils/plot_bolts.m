function plot_bolts(bolts,params)

axis equal;
hold on;

nbolts = length(bolts);
for i = 1:nbolts
    L = line(bolts{i}(:,1), bolts{i}(:,2), bolts{i}(:,3));
end

