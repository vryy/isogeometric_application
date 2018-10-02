%% create a sample NURBS surface for testing
function srf = create_sample_nurbs_surface

srf = nrb4surf([0.0 0.0 0.5],[1.0 0.0 -0.5],[0.0 1.0 -0.5],[1.0 1.0 0.5]);
srf = nrbdegelev(srf, [1 1]);

sizes = size(srf.coefs);

for i = 1:sizes(2)
    for j = 1:sizes(3)
        srf.coefs(4,i,j) = 1.0/(1.0 + sin(i)*cos(j)/37.0);
    end
end

srf.coefs(3,2,1) = 0.3;
srf.coefs(3,2,3) = 0.3;
srf.coefs(3,3,2) = 0.3;
srf.coefs(3,1,2) = 0.3;

end

