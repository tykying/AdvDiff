m = 4*4;

cb1 = zeros(1, m);
cb2 = zeros(1, m);
v = 1:m;
for i = 1: m
    if (mod(i, 2) == 0) 
        cb1(i) = 1;
        cb2(i) = 0;
    else 
        cb1(i) = 0;
        cb2(i) = 1;
    end
end

v_b = 0.5*(v(1:end-1) + v(2:end))

err = 5*cb1 - 2*cb2;

v = v + err;

mode1 = dot(v, cb1);
mode2 = dot(v, cb2);

v_b = 0.5*(v(1:end-1) + v(2:end))
