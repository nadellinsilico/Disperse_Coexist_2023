function xx=trajectories2(x0,p)

CPA = [1 1 0];
CEF = [0 1 1];
CEC = [1 0 1];

opts = odeset('Events',@(t,x) events(t,x,p));
[tt,xx]=ode113(@(t,x) equations(t,x,p),[0 800],x0,opts);

tt=tt/log(2);

hold on
plot(tt,xx(:,1),'Color',CPA,'LineWidth',3)
plot(tt,xx(:,2),'Color',CEF,'LineWidth',3)
plot(tt,xx(:,3),'Color',CEC,'LineWidth',3)
xlabel('Time (doublings)','FontSize',20)
ylabel('Fraction of population','FontSize',20)
legend('Producer','Resistant','Sensitive','FontSize',15)
axis([0 tt(end) 0 1.1])
set(gca,'FontSize',15)

function dx=equations(t,x,p)

% populations rescaled by highest carrying capacity (Kpa=1)

% growth of each species
gpa=0.042;
%gef=0.005;
gef=0.109;
gec=0.036;

% carrying capacity for each species
Kpa=0.71981117;
Kef=0.194193739;
Kec=0.607803022;


% dispersal of each species
spa=2;
sef=0;
sec=0;
s=p(1); % triggers dispersal

pa=x(1);
ef=x(2);
ec=x(3);

dpa=gpa*pa.*(1-(pa+(gef/gpa)*ef+(gec/gpa)*ec)/Kpa)-s*spa*pa;
def=gef*ef.*(1-(ef+(gpa/gef)*pa+(gec/gef)*ec)/Kef)-s*sef*ef;
dec=gec*ec.*(1-(ec+(gpa/gec)*pa+(gef/gec)*ef)/Kec)-s*sec*ec;

dx=[dpa def dec]';

end

function [value,isterminal,direction]=events(t,x,p)
    dx=equations(t,x,p);
    value=norm(dx./x)-1e-3;
    isterminal=1;
    direction=0;
end
end