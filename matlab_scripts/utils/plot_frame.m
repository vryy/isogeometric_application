%% Plot the Frenet frame
function plot_frame(C,B,N,T,s)

quiver3(C(1),C(2),C(3),s*B(1),s*B(2),s*B(3),'r','LineWidth',2,'MaxHeadSize',1.0);
quiver3(C(1),C(2),C(3),s*N(1),s*N(2),s*N(3),'g','LineWidth',2,'MaxHeadSize',1.0);
quiver3(C(1),C(2),C(3),s*T(1),s*T(2),s*T(3),'b','LineWidth',2,'MaxHeadSize',1.0);

end

