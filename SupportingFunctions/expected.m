function E = expected (x,y,params)
% The expected counts per pixel.
%E = params (10).*twoDGauss (x,y,params (1),params (3),params (5),params (7)) + params (11).*twoDGauss (x,y,params (2),params (4),params (6),params (8)) + params (9).^2;
E = params (10).*twoDGauss (x,y,params (1),params (3),params (5),params (7)) + params (11).*twoDGauss (x,y,params (2),params (4),params (6),params (8)) + params (9);
end