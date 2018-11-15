filename=('particle_tracer.csv');
M = csvread(filename);
M = sortrows(M,1);
number_of_particles = 10;
timesteps = 20000;
if_statement_modulo = 1;
datapoints = timesteps / if_statement_modulo;


figure
for i = 0:1:number_of_particles-1
    plot3(M(1+i*datapoints:datapoints*i+datapoints,3),M(1+i*datapoints:datapoints*i+datapoints,4),M(1+i*datapoints:datapoints*i+datapoints,5)); hold on; 
    %plot(M(1+i*datapoints:datapoints*i+datapoints,3),M(1+i*datapoints:datapoints*i+datapoints,4)); hold on;
end
set(gca,'TickLabelInterpreter','latex','Fontsize',14)

ylabel('$ y / d_{\rm p} $','Interpreter','latex','Fontsize',20);
xlabel('$ x / d_{\rm p} $','Interpreter','latex','Fontsize',20);
zlabel('$ z / d_{\rm p} $','Interpreter','latex','Fontsize',20);

