%x0=[0.01 0.001 0.1]; % p r s
x0=[0.1 0.1 0.1];
p=1;
p2=p;
much={'Low','Moderate','High'};

subplot(2,3,1)
xx=trajectories2(x0, p);

subplot(2,3,1)
triplot2(xx, p2)
title(['Density'],'FontSize',30)


figure
names={'traj1','traj2b','traj2','traj3'};
for i=1:4
    
    load(names{i})

    subplot(4,4,i)
    diffplot(tt,pp,rr,ss,10)
    subplot(4,4,4+i)
    diffplot(tt,pp,rr,ss,100)
    subplot(4,4,8+i)
    diffplot(tt,pp,rr,ss,300)
    subplot(4,4,12+i)
    diffplot(tt,pp,rr,ss,600)

end