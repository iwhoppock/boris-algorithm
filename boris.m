clear
close

dt = 1e-2;
mass = 1.0;
charge = 1.0;
vAc = 3e-4;

duration = 5000;

v = [0, 1, 0];
x = [-1, 0, 0];

B = [0, 0, 1];
E = [0, 0, 0];

X = zeros(duration,3);
V = zeros(duration,3);


for time = 1:1:duration
    t = charge ./ mass .* B .* 0.5 .* dt;
    s = 2 .* t ./ (1 + t.*t);
    v_minus = v + charge ./ (mass * vAc) .* E .* 0.5 .* dt; 
    v_prime = v_minus + cross(v_minus,t);
    v_plus = v_minus + cross(v_prime,s);
    v = v_plus + charge ./ (mass * vAc) .* E .* 0.5 .* dt;
    x = x + v .* dt;
    X(time,:) = x;
    V(time,:) = v;
end    


plot(X(:,1),X(:,2),'k','Linewidth',2); hold on;
%plot3(X(:,1),X(:,2),X(:,3),'k','Linewidth',2); hold on;
set(gca,'TickLabelInterpreter','latex','Fontsize',14)

ylabel('$ y / d_{\rm p} $','Interpreter','latex','Fontsize',20);
xlabel('$ x / d_{\rm p} $','Interpreter','latex','Fontsize',20);
%zlabel('$ z / d_{\rm p} $','Interpreter','latex','Fontsize',20);

%title('$ title $','Interpreter','latex','Fontsize',20);


