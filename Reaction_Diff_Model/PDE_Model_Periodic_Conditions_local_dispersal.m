
N = 51;
Tend = 40000;
D = 0.000025;
rng(1769)

t = [];
ch1 = [];
ch2 = [];
ch3 = [];

%
dx = 1 / (N - 1);
dt = 0.1 * (0.5 * dx^2 / D);
k = dt / (2 * dx^2);

%
x = [0:dx:1]; y = [0:dx:1];


pa = 0.3 * (rand(N, N) < 0.1);
ef = 0.3 * (rand(N, N) < 0.1);
ec = 0.3 * (rand(N, N) < 0.1);

pa = reshape(pa, N^2, 1);
ef = reshape(ef, N^2, 1);
ec = reshape(ec, N^2, 1);
s = zeros(N^2, 1);
time = 0;

%
A1D = diffmat(N); Imat = eye(N);
A2D = kron(A1D, Imat) + kron(Imat, A1D);
Amat = eye(N^2) - D * k * A2D;

hfig = figure;
v = VideoWriter('pa_dispersal_local.mp4', 'MPEG-4');
v.FrameRate = 5;
open(v);
Nt = int16(Tend / dt);
F(Nt) = struct('cdata', [], 'colormap', []);

CPA = [1 1 0];
CEF = [0 1 1];
CEC = [1 0 1];

%
for i = 1:Nt
    t(end + 1) = time;
    ch1(end + 1) = [sum(pa)];
    ch2(end + 1) = [sum(ef)];
    ch3(end + 1) = [sum(ec)];

    pa = Amat \ (pa + D * k * A2D * pa);
    ef = Amat \ (ef + D * k * A2D * ef);
    ec = Amat \ (ec + D * k * A2D * ec);

    [dxx, ss] = FF(0, [pa ef ec], s);
    pa = pa + dt * dxx(:, 1);
    ef = ef + dt * dxx(:, 2);
    ec = ec + dt * dxx(:, 3);

    time = time + dt;
    s = ss;

    % carrying capacity for each species
    Kpa = 0.71981117;
    Kef = 0.194193739;
    Kec = 0.607803022;
    xpa = reshape(pa, N, N);
    xef = reshape(ef, N, N);
    xec = reshape(ec, N, N);
    C(:, :, 1) = min(1, xpa * CPA(1) / Kpa + xef * CEF(1) / Kef + xec * CEC(1) / Kec);
    C(:, :, 2) = min(1, xpa * CPA(2) / Kpa + xef * CEF(2) / Kef + xec * CEC(2) / Kec);
    C(:, :, 3) = min(1, xpa * CPA(3) / Kpa + xef * CEF(3) / Kef + xec * CEC(3) / Kec);
    image(x, y, C);
    text(0.8, 0.1, ['time ' num2str(time)], 'Color', 'w')
    drawnow;

    F(i) = getframe(hfig);
    writeVideo(v, F(i));

end
close(v);

function [A] = diffmat(n)

A = zeros(n);

for i = 1:n
    A(i, i) = -2.0;
end

for i = 2:n
    A(i, i - 1) = 1.0;
    A(i - 1, i) = 1.0;
end

% Add periodic boundary conditions 
A(1, end) = 1.0;
A(end, 1) = 1.0;

return

end

function [dx ss] = FF(t, x, s)

% growth of each species
gpa = 0.042;
gef = 0.109;
gec = 0.036;

% carrying capacity for each species
Kpa = 0.71981117;
Kef = 0.194193739;
Kec = 0.607803022;

% dispersal of each species
spa = 0;
sef = 0;
sec = 0;

pa_thresh_up = 0.9 * Kpa;
pa_thresh_down = 0.1 * Kpa;

pa = x(:, 1);
ef = x(:, 2);
ec = x(:, 3);

dpa = gpa * pa .* (1 - (pa + (gef / gpa) * ef + (gec / gpa) * ec) / Kpa) - spa * (s .* pa);
def = gef * ef .* (1 - (ef + (gpa / gef) * pa + (gec / gef) * ec) / Kef) - sef * (s .* ef);
dec = gec * ec .* (1 - (ec + (gpa / gec) * pa + (gef / gec) * ef) / Kec) - sec * (s .* ec);

dx = [dpa def dec];

for i = 1:length(pa)
    if pa(i) >= pa_thresh_up && s(i) == 0
        ss(i) = 1;
    elseif pa(i) <= pa_thresh_down && s(i) == 1
        ss(i) = 0;
    else
        ss(i) = s(i);
    end
end
ss = ss';

end