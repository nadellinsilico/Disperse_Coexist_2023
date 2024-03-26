% growth of each species
gpa=0.042;
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

% PA thresholds for turning dispersal on/off
pa_thresh_up=0.9*Kpa;
pa_thresh_down=0.1*Kpa;

v=[pa_thresh_up pa_thresh_down gpa gef gec Kpa Kef Kec spa sef sec];

bias=[0.01:0.05:2];
KK = [0.031:0.00075:0.06075];

for i=1:length(bias)
    for j=1:length(KK)
        vv=v;
        vv(10:11)=bias(i);
        vv(3)=KK(j);
        M(i,j,:)=mixed_pd_phase(vv);
        [i j]
    end
end

imagesc(KK,bias,M)
colorbar()
caxis([0,6]);
set(gca,'YDir','normal')
xlabel('gPA')
ylabel('Bias')