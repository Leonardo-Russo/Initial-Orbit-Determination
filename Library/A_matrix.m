function A = A_matrix(r_vect)

global mu

x = r_vect(1);
y = r_vect(2);
z = r_vect(3);

r = norm(r_vect);

A = [0    0    0    1    0    0;
     0    0    0    0    1    0;
     0    0    0    0    0    1;
     3*mu*x^2/r^5 - mu/r^3  3*mu*x*y/r^5  3*mu*x*z/r^5  0   0   0;
     3*mu*x*y/r^5  3*mu*y^2/r^5 - mu/r^3  3*mu*y*z/r^5  0   0   0;
     3*mu*x*z/r^5  3*mu*y*z/r^5  3*mu*z^2/r^5 - mu/r^3  0   0   0];

end