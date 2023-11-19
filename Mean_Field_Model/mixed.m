CPA = [1 1 0];
CEF = [0 1 1];
CEC = [1 0 1];

% PA thresholds for turning dispersal on/off
pa_thresh_up=0.9*0.71981117;
pa_thresh_down=0.1*0.71981117;

T=800;
ttt=[];
xxx=[];
p=[0 pa_thresh_up];
t0=0;
x0=[0.02 0.01 0.01];

while t0<T

    opts = odeset('Events',@(t,x) events(t,x,p));
    [tt,xx]=ode15s(@(t,x) FF(t,x,p),[t0 T],x0,opts);
    t0=tt(end);
    x0=xx(end,:);
    ttt=[ttt; tt];
    xxx=[xxx; xx];
    
    if p==[0 pa_thresh_up]
        p=[1 pa_thresh_down];
    elseif p==[1 pa_thresh_down]
        p=[0 pa_thresh_up];
    end
    
end

hold on
plot(ttt,xxx(:,1),'Color',CPA,'LineWidth',3)
plot(ttt,xxx(:,2),'Color',CEF,'LineWidth',3)
plot(ttt,xxx(:,3),'Color',CEC,'LineWidth',3)
xlabel('Time (doublings)','FontSize',20)
ylabel('Fraction of population','FontSize',20)
legend('PA','EF','EC','FontSize',15)
axis([0 tt(end) 0 1.1])
set(gca,'FontSize',15)
ch1 = xxx(:,1,:);
ch2 = xxx(:,2,:);
ch3 = xxx(:,3,:);

function dx=FF(t,x,p)

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
spa=0;
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

function [position,isterminal,direction] = events(t,x,p)
  position = x(1)-p(2); % The value that we want to be zero
  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
end

