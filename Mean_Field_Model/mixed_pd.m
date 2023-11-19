function xxx=mixed_pd(v)

% PA thresholds for turning dispersal on/off
pa_thresh_up=v(1);
% disp(pa_thresh_up)
pa_thresh_down=v(2);

T=4000;
xxx=[NaN];
p=[0 pa_thresh_up v];
t0=0;
x0=[0.01 0.01 0.01];

while t0<T
    
    opts = odeset('Events',@(t,x) events(t,x,p));
    [tt,xx]=ode15s(@(t,x) FF(t,x,p),[t0 T],x0,opts);
    t0=tt(end);
    x0=xx(end,:);
    
    if p==[0 pa_thresh_up v]
        p=[1 pa_thresh_down v];
        
    elseif p==[1 pa_thresh_down v]
        p=[0 pa_thresh_up v];        
    xxx=xx(end,1)/(xx(end,1)+xx(end,2)+xx(end,3));    
        
    
    end
    
end

function dx=FF(t,x,p)

% populations rescaled by highest carrying capacity (Kpa=1)

% growth of each species
gpa=p(5);
gef=p(6);
gec=p(7);

% carrying capacity for each species
Kpa=p(8);
Kef=p(9);
Kec=p(10);

% dispersal of each species
spa=p(11);
sef=p(12);
sec=p(13);
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

end