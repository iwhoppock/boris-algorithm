function [x,v] = boris_rotaion(x,v,charge,mass,vAc,dt,B,E)
    t = charge ./ mass .* B .* 0.5 .* dt;
    s = 2 .* t ./ (1 + t.*t);
    v_minus = v + charge ./ (mass * vAc) .* E .* 0.5 .* dt; 
    v_prime = v_minus + cross(v_minus,t);
    v_plus = v_minus + cross(v_prime,s);
    v = v_plus + charge ./ (mass * vAc) .* E .* 0.5 .* dt;
    x = x + v .* dt;
end