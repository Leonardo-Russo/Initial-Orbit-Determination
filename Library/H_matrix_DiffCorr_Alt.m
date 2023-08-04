function H = H_matrix_DiffCorr_Alt(X, Rgs)
% WARNING: this function may contain unfixed errors, please consult the
% H_matrix function for a correct application of the differential
% correction theory for the construction of the H matrix.

delta = 1e-10;
dX = delta * X;

H = zeros(2,6);

[Rax, Decx] = ECI2TopoRaDec(Rgs, X(1:3));
[Rax_p, Decx_p] = ECI2TopoRaDec(Rgs, X(1:3)+[dX(1) 0 0]');

dRa_dx = deg2rad(Rax_p-Rax)/delta;
dDec_dx = deg2rad(Decx_p-Decx)/delta;


[Ray, Decy] = ECI2TopoRaDec(Rgs, X(1:3));
[Ray_p, Decy_p] = ECI2TopoRaDec(Rgs, X(1:3)+[0 dX(2) 0]');

dRa_dy = deg2rad(Ray_p-Ray)/delta;
dDec_dy = deg2rad(Decy_p-Decy)/delta;


[Raz, Decz] = ECI2TopoRaDec(Rgs, X(1:3));
[Raz_p, Decz_p] = ECI2TopoRaDec(Rgs, X(1:3)+[0 0 dX(3)]');

dRa_dz = deg2rad(Raz_p-Raz)/delta;
dDec_dz = deg2rad(Decz_p-Decz)/delta;


H(1, 1) = dRa_dx;
H(1, 2) = dRa_dy;
H(1, 3) = dRa_dz;
H(2, 1) = dDec_dx;
H(2, 2) = dDec_dy;
H(2, 3) = dDec_dz;


end