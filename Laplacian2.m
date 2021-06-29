function lap=Laplacian2(vtg)
k_i = sum(sum(vtg,3),2)/2;
k_ij = sum(vtg,3);
lap = 2*diag(k_i) - k_ij;
end