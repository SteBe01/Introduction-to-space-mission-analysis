% This script is used to create and export the graph 
% used for maneuver comparison

clear, clc
close all

StandardDv=[7.6009 33329.2874];
StandardDt=[7.7407 25639.3816];
InvStandardDv=[9.0096 23114.5090];
InvStandardDt=[9.0384 22630.7839];
Bielliptic=[7.1918 66417.6263];
Dv_min=[6.6324 7971.6116];
Direct=[14.3652 6700.1589];
Intersection=[7.0185 7816.0531];

h=figure;
plot(StandardDv(2),StandardDv(1),'xr',LineWidth=1.5);
text(StandardDv(2)+1500,StandardDv(1),"Standard (v)","FontSize",13)
hold on, grid minor
xlabel("Time [s]",FontSize=20)
ylabel("\Deltav [km/s]",FontSize=20)

plot(StandardDt(2),StandardDt(1),'xr',LineWidth=1.5);
text(StandardDt(2)-16000,StandardDt(1)+0.05,"Standard (t)","FontSize",13)

plot(InvStandardDv(2),InvStandardDv(1),'xr',LineWidth=1.5);
text(InvStandardDv(2)+2000,InvStandardDv(1),"InvStandard (v)","FontSize",13)

plot(InvStandardDt(2),InvStandardDt(1),'xr',LineWidth=1.5);
text(InvStandardDt(2)-20000,InvStandardDt(1),"InvStandard (t)","FontSize",13)

plot(Bielliptic(2),Bielliptic(1),'xr',LineWidth=1.5);
text(Bielliptic(2)-12000,Bielliptic(1),"Bielliptic","FontSize",13)

plot(Dv_min(2),Dv_min(1),'xr',LineWidth=1.5);
text(Dv_min(2)+1000,Dv_min(1)-0.08,"Direct \Deltav_{min}","FontSize",13)

plot(Direct(2),Direct(1),'xr',LineWidth=1.5);
text(Direct(2)+1000,Direct(1),"Direct","FontSize",13)

plot(Intersection(2),Intersection(1),'xr',LineWidth=1.5);
text(Intersection(2)+1000,Intersection(1),"Intersection","FontSize",13)


%%

clear, clc
close all

StandardDv=[7.6009 33329.2874];
StandardDt=[7.7407 25639.3816];
InvStandardDv=[9.0096 23114.5090];
InvStandardDt=[9.0384 22630.7839];
Bielliptic=[7.1918 66417.6263];
Dv_min=[6.6324 7971.6116];
Direct=[14.3652 6700.1589];
Intersection=[7.0185 7816.0531];

h=figure;
plot(StandardDv(2),StandardDv(1),'xr',LineWidth=1.5);
text(StandardDv(2)+1000,StandardDv(1),"Standard (v)","FontSize",13)
hold on, grid minor
xlabel("Time [s]",FontSize=20)
ylabel("\Deltav [km/s]",FontSize=20)

plot(StandardDt(2),StandardDt(1),'xr',LineWidth=1.5);
text(StandardDt(2)-10000,StandardDt(1),"Standard (t)","FontSize",13)

plot(InvStandardDv(2),InvStandardDv(1),'xr',LineWidth=1.5);
text(InvStandardDv(2)+1000,InvStandardDv(1),"InvStandard (v)","FontSize",13)

plot(InvStandardDt(2),InvStandardDt(1),'xr',LineWidth=1.5);
text(InvStandardDt(2)-12000,InvStandardDt(1),"InvStandard (t)","FontSize",13)

plot(Bielliptic(2),Bielliptic(1),'xr',LineWidth=1.5);
text(Bielliptic(2)-8000,Bielliptic(1),"Bielliptic","FontSize",13)

plot(Dv_min(2),Dv_min(1),'xr',LineWidth=1.5);
text(Dv_min(2)+1000,Dv_min(1),"Direct \Deltav_{min}","FontSize",13)

plot(Direct(2),Direct(1),'xr',LineWidth=1.5);
text(Direct(2)+1000,Direct(1),"Direct","FontSize",13)

plot(Intersection(2),Intersection(1),'xr',LineWidth=1.5);
text(Intersection(2)+1000,Intersection(1),"Intersection","FontSize",13)


%%
set(h,'PaperOrientation','landscape');
set(h,'PaperPosition', [1 1 28 19]);
print(gcf, 'Graph.pdf', '-dpdf', '-fillpage')

%%
print(gcf, 'Graph1.pdf', '-dpdf')