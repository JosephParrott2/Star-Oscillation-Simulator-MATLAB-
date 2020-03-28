%% Part A.1
r=1; % R/Rsun cancels out
m=1; % M/Msun cancels out
conv=1e5; % converting the frequency into hearing range
teff=5777;
dv=134.9*(m^-0.5)*(r^-1.5)*1e-6*conv;
Vmax=3090*m*(r^-2)*(teff/5777)^-0.5*1e-6*conv;

fs = 44100; % standard sampling rate
T = 1/fs; % sampling period
t = 0:T:5; % time vector
Xtotal=0; % 0 is the starting value
for n=-2:1:2 % does the equation for each frequency
    phi=2*pi*rand(1,1); % random phase
X=sin(2*pi*(Vmax+n*dv)*t+phi);
Xtotal=Xtotal+X; % adds all of the frequencies together since they're interfering
end


sound(Xtotal, fs); % play the signal

%% Part A.2

load star_evolving.mat
m=1;
teff=5777;
conv=1e5;

for r=[radius] % Doing the process for each radius
    figure % makes it so each set of subplots for each radius is on a new figure
    m=1;
Vmax=3090*m*(r^-2)*(teff/5777)^-0.5*1e-6;
dv=134.9*(m^-0.5)*(r^-1.5)*1e-6; % Defining Vmax and dv when the mass cancels but there are different radii
t=linspace(0,0.16e5,100000) % Performs the action between t=0 and 1.6e4
Xtotal=0

Xtotalconverted=0
fs = 44100; % standard sampling rate
T = 1/fs; % sampling period
tC = 0:T:5; % This is the time, with a conversion 1e5

for n=-2:1:2 % Much like A.1 we're using a loop to cycle through each frequency
    phi=2*pi*rand(1,1); % random phase
X=sin(2*pi*(Vmax+n*dv)*t+phi);
Xtotal=Xtotal+X;
subplot(2,2,1)
plot(t,X)
xlabel('Time (s)')
ylabel('Amplitude: Arbitrary units')
title('Time series: Individual Modes')
hold on

Xconverted=sin(2*pi*(Vmax*conv+n*dv*conv)*tC+phi);
Xtotalconverted=Xtotalconverted+Xconverted;
end

subplot(2,2,2)
plot(t,Xtotal)
xlabel('Time (s)')
ylabel('Amplitude: Arbitrary units')
title('Time series: Collective Behaviour')

subplot(2,2,3)
scatter(age, radius, 'g','fill')
xlim([0 12e9])
ylim([0 10])
xlabel('Time(s)')
ylabel('R/Rs')
title('Radius plotted against age')
grid on

subplot(2,2,4)
[x,y,z] = sphere;
surf(r*x,r*y,r*z) % sphere centred at the origin
axis square
rscale=5;
xlim([-rscale rscale])
ylim([-rscale rscale])
zlim([-rscale rscale])
shading interp
xlabel('X(Rs)')
ylabel('Y(Rs)')
zlabel('Z(Rs)')

sound(Xtotalconverted,fs)
pause(5) % The pause is so you don't have all of the different sounds playing at once

end
% QuestionA2: As the star gets older, the frequency of the sound decreases,
% therefore the pitch decreases as well. This is because the radius
% increases, therefore the standing waves created within the star having a
% longer wavelength, hence a lower frequency.
%% Part B.1

syms r m % Making it so r and m are variables rather than values

[solm solr]=solve(1e-6*(1954-3090*m*(r^-2))==0, 1e-6*(95.66-134.9*(m^-0.5)*(r^-1.5))==0); % Solving the equation for Vmax and dv
RatioM=double(solm) % This is the ratio of a star's mass compared to the sun's
RatioR=double(solr) % This is the ratio of a star's radius compared to the sun's

Rsun=6.9959e5 
Msun=1.989e30

Rstar=double(solr*Rsun) % This is using the ratio and the value of Rsun to find the radius of the star
Mstar=double(solm*Msun) % This is using the ratio and the value of Rsun to find the mass of the star

Rearth=6.371e3
d=1.5e-5
Rplanet=(sqrt(d)*Rstar)/Rearth % Using the radius of the star to work out the radius of the planet

%% Part B.2


load star_evolving.mat
y=transpose(age) % Transpose changes matrix dimensions from 1x6 to 6x1
x=transpose(radius) % Transpose changes matrix dimensions from 1x6 to 6x1


n = length(x);
s = (1:n)';
t = (1:.05:n)';
u = splinetx(s,x,t); % Using splinetx to create a smooth line interpolation connecting points
v = splinetx(s,y,t);
clf reset
plot(x,y,'.',u,v,'-');


