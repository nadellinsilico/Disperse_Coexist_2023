function triplot_heatmap_line(xx, p2)

%run mixed.m, then call triplot_heatmap_line(xxx, dispersal on/ff (1/0))

CPA = [1 1 0];
CEF = [0 1 1];
CEC = [1 0 1];

n=10;
t = [[0, 0.5, 1]; [0, sqrt(3)/2, 0]];

ng = ( ( n + 1 ) * ( n + 2 ) ) / 2;
tg = zeros ( 2, ng );

p = 0;

for i = 0 : n
for j = 0 : n - i
  k = n - i - j;
  p = p + 1;
  tg(1:2,p) = ( i * t(1:2,1) + j * t(1:2,2) + k * t(1:2,3) ) / n;
end
end

abundances=@(x,y) [1-x-y/sqrt(3), x-y/sqrt(3), 2*y/sqrt(3)];
components=@(a,b,c) [3*(b-a)/4, sqrt(3)*(2*c-b-a)/4];

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
s=p2(1); % triggers dispersal

X=@(pa,ef,ec) pa*(1-(pa+(gef/gpa)*ef+(gec/gpa)*ec)/Kpa)-s*spa*pa;
dp =@(pa,ef,ec) gpa*pa.*(1-(pa+(gef/gpa)*ef+(gec/gpa)*ec)/Kpa)-s*spa*pa;
dr =@(pa,ef,ec) gef*ef.*(1-(ef+(gpa/gef)*pa+(gec/gef)*ec)/Kef)-s*sef*ef;
ds =@(pa,ef,ec) gec*ec.*(1-(ec+(gpa/gec)*pa+(gef/gec)*ef)/Kec)-s*sec*ec;

N=length(tg);
aa=zeros(3,1);
cc=zeros(2,N);
for n=1:N
    aa=abundances(tg(1,n),tg(2,n));
    
    pa=aa(3);
    ef=aa(2);
    ec=aa(1);
    
    pp=dp(aa(3),aa(2),aa(1));
    rr=dr(aa(3),aa(2),aa(1));
    ss=ds(aa(3),aa(2),aa(1));

    ppp=pp-pa*(pp+rr+ss);
    rrr=rr-ef*(pp+rr+ss);
    sss=ss-ec*(pp+rr+ss);
    
    cc(:,n)=components(sss,rrr,ppp);

end

w=sum(xx,2);
for i=1:3
xx(:,i)=xx(:,i)./w;
end
position=@(a,b,c) [b+c/2, sqrt(3)*c/2];
xy=position(xx(:,3),xx(:,2),xx(:,1));

hold on
plot([0 0.5],[0 sqrt(3)/2],'k')
plot([0 1],[0 0],'k')
plot([0.5 1],[sqrt(3)/2 0],'k')

plot(0.5,sqrt(3)/2,'o','Color',CPA,'MarkerSize',15,'LineWidth',5)
plot(1,0,'o','Color',CEF,'MarkerSize',15,'LineWidth',5)
plot(0,0,'o','Color',CEC,'MarkerSize',15,'LineWidth',5)

time = uint64(1):uint64(length(xy));

plot(xy(:,1),xy(:,2),'Color','black','LineWidth',1)
scatter(xy(:,1),xy(:,2), 30, time, 'filled');
alpha(.5)

colormap(summer);
colorbar;

axis off
set(gcf, 'renderer', 'painters')
