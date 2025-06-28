% 加载数据
D2R=pi/180.0;
R2D=180.0/pi;
a=6378137.0;
e=0.08181919104;

data = importdata('Result.txt');
result=data.data;
data2=importdata('Res-17.pos');
res=data2.data;

time=result(:,1);
lat=result(:,2);
lon=result(:,3);
h=result(:,4);
vn=result(:,5);
ve=result(:,6);
vd=result(:,7);
roll=result(:,8);
pitch=result(:,9);
heading=result(:,10);

t=result(:,1);
latstd=result(:,11);
lonstd=result(:,12);
hstd=result(:,13);
vnstd=result(:,14);
vestd=result(:,15);
vdstd=result(:,16);
rollstd=result(:,17);
pitchstd=result(:,18);
headingstd=result(:,19);

time2=res(:,1);
lat2=res(:,2);
lon2=res(:,3);
h2=res(:,4);
vn2=res(:,5);
ve2=res(:,6);
vd2=-res(:,7);
roll2=res(:,8);
pitch2=res(:,9);
heading2=res(:,10);

 % 找出共同的时间范围
min_time = max(min(time), min(time2));
max_time = min(max(time), max(time2));
% 筛选在共同时间范围内的数据
idx = (time >= min_time) & (time <= max_time);
time=time(idx);
lat=lat(idx);
lon=lon(idx);
h=h(idx);
vn=vn(idx);
ve=ve(idx);
vd=vd(idx);
roll=roll(idx);
pitch=pitch(idx);
heading=heading(idx);
idx2=(time2>=min_time)&(time<=max_time);
time2=time2(idx2);
lat2=lat2(idx2);
lon2=lon2(idx2);
h2=h2(idx2);
vn2=vn2(idx2);
ve2=ve2(idx2);
vd2=vd2(idx2);
roll2=roll2(idx2);
pitch2=pitch2(idx2);
heading2=heading2(idx2);

%转NED
lat=lat*D2R;
lon=lon*D2R;
lat2=lat2*D2R;
lon2=lon2*D2R;
firstlat=lat(1);
firstlon=lon(1);
firsth=h(1);
firstlat2=lat2(1);
firstlon2=lon2(1);
firsth2=h2(1);
RM1 = a*(1-e^2)/sqrt((1-e^2*(sin(firstlat))^2)^3);
RN1 = a/sqrt(1-e^2*(sin(firstlat))^2);
DR1 = diag([RM1 + firsth, (RN1 + firsth)*cos(firsth), -1]);
RM2 = a*(1-e^2)/sqrt((1-e^2*(sin(firstlat2))^2)^3);
RN2 = a/sqrt(1-e^2*(sin(firstlat2))^2);
DR2 = diag([RM2 + firsth2, (RN2 + firsth2)*cos(firsth2), -1]);
% blh to ned
num=size(lat,1);
pos1 = zeros(num,3);
pos2=zeros(num,3);
for i = 1:size(pos1, 1)
    delta_blh1 = [lat(i)-firstlat,lon(i)-firstlon,h(i)-firsth];
    delta_pos1 = DR1 * delta_blh1';
    pos1(i, :) = delta_pos1';
    
    delta_blh2 = [lat2(i)-firstlat2,lon2(i)-firstlon2,h2(i)-firsth2];
    delta_pos2 = DR2 * delta_blh2';
    pos2(i, :) = delta_pos2';
end

dh=h-h2;
dpos=pos1-pos2;
dvn=vn-vn2;
dve=ve-ve2;
dvd=vd-vd2;
droll=roll-roll2;
dpitch=pitch-pitch2;
dheading=heading-heading2;
dheading=mod(dheading + 180, 360) - 180;

figure;
plot(pos1(:,2),pos1(:,1),pos2(:,2),pos2(:,1));
grid on;
legend('编程结果','PosMind结果');
xlabel('E(m)');
ylabel('N(m)');

figure;
subplot(3,1,1);
plot(time,dpos(:,1));
grid on;
xlabel('Time(s)');
ylabel('dN(m)');

subplot(3,1,2);
plot(time,dpos(:,2));
grid on;
xlabel('Time(s)');
ylabel('dE(m)');

subplot(3,1,3);
plot(time,dpos(:,3));
grid on;
xlabel('Time(s)');
ylabel('dH(m)');

figure;
subplot(3,1,1);
plot(time,dvn);
grid on;
xlabel('Time(s)');
ylabel('dVn(m/s)');

subplot(3,1,2);
plot(time,dve);
grid on;
xlabel('Time(s)');
ylabel('dVe(m/s)');

subplot(3,1,3);
plot(time,dvd);
grid on;
xlabel('Time(s)');
ylabel('dVd(m/s)');

figure;
subplot(3,1,1);
plot(time,droll);
grid on;
xlabel('Time(s)');
ylabel('dRoll(deg)');

subplot(3,1,2);
plot(time,dpitch);
grid on;
xlabel('Time(s)');
ylabel('dPitch(deg)');

subplot(3,1,3);
plot(time,dheading);
grid on;
xlabel('Time(s)');
ylabel('dHeading(deg)');

%%
%std
figure;
subplot(3,1,1);
grid on;
plot(t,latstd);
ylabel('latSTD(rad)');

subplot(3,1,2);
grid on;
plot(t,lonstd);
ylabel('lonSTD(rad)');

subplot(3,1,3);
grid on;
plot(t,hstd);
ylabel('hSTD(m)');
xlabel('Time(s)');

figure;
subplot(3,1,1);
grid on;
plot(t,vnstd);
ylabel('VnSTD(m/s)');

subplot(3,1,2);
grid on;
plot(t,vestd);
ylabel('VeSTD(m/s)');

subplot(3,1,3);
grid on;
plot(t,vdstd);
ylabel('VdSTD(m/s)');
xlabel('Time(s)');

figure;
subplot(3,1,1);
grid on;
plot(t,rollstd);
ylabel('RollSTD(Deg)');

subplot(3,1,2);
grid on;
plot(t,pitchstd);
ylabel('PitchSTD(Deg)');

subplot(3,1,3);
grid on;
plot(t,headingstd);
ylabel('HeadingSTD(Deg)');
xlabel('Time(s)');

%%计算RMS
dN_rms=rms(dpos(:,1));
dE_rms=rms(dpos(:,2));
dH_rms=rms(dpos(:,3));
dvn_rms=rms(dvn);
dve_rms=rms(dve);
dvd_rms=rms(dvd);
droll_rms=rms(droll);
dpitch_rms=rms(dpitch);
dheading_rms=rms(dheading);