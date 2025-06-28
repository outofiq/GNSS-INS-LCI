D2R=pi/180.0;
a=6378137.0;
e=0.08181919104;

dataL = importdata('LCI.pos');
resultL=dataL.data;
dataT=importdata('TCI.pos');
resultT=dataT.data;

timeL=resultL(:,1);
latL=resultL(:,2);
lonL=resultL(:,3);
hL=resultL(:,4);
vnL=resultL(:,5);
veL=resultL(:,6);
vuL=resultL(:,7);
rollL=resultL(:,8);
pitchL=resultL(:,9);
headingL=resultL(:,10);

timeT=resultT(:,1);
latT=resultT(:,2);
lonT=resultT(:,3);
hT=resultT(:,4);
vnT=resultT(:,5);
veT=resultT(:,6);
vuT=resultT(:,7);
rollT=resultT(:,8);
pitchT=resultT(:,9);
headingT=resultT(:,10);

dh=hL-hT;
dvn=vnL-vnT;
dve=veL-veT;
dvu=vuL-vuT;
droll=rollL-rollT;
dpitch=pitchL-pitchT;
dheading=headingL-headingT;
dheading=mod(dheading+180,360)-180;

%%
%转NED
latL=latL*D2R;
lonL=lonL*D2R;
latT=latT*D2R;
lonT=lonT*D2R;
firstlatL=latL(1);
firstlonL=lonL(1);
firsthL=hL(1);
firstlatT=latT(1);
firstlonT=lonT(1);
firsthT=hT(1);
RML = a*(1-e^2)/sqrt((1-e^2*(sin(firstlatL))^2)^3);
RNL = a/sqrt(1-e^2*(sin(firstlatL))^2);
DRL = diag([RML + firsthL, (RNL + firsthL)*cos(firsthL), -1]);
RMT = a*(1-e^2)/sqrt((1-e^2*(sin(firstlatT))^2)^3);
RNT = a/sqrt(1-e^2*(sin(firstlatT))^2);
DRT = diag([RMT + firsthT, (RNT + firsthT)*cos(firsthT), -1]);
% blh to ned
num=size(latL,1);
posL = zeros(num,3);
posT=zeros(num,3);
for i = 1:size(posL, 1)
    delta_blhL = [latL(i)-firstlatL,lonL(i)-firstlonL,hL(i)-firsthL];
    delta_posL = DRL * delta_blhL';
    posL(i, :) = delta_posL';
    
    delta_blhT = [latT(i)-firstlatT,lonT(i)-firstlonT,hT(i)-firsthT];
    delta_posT = DRT * delta_blhT';
    posT(i, :) = delta_posT';
end
dN=posL(:,1)-posT(:,1);
dE=posL(:,2)-posT(:,2);

figure;
grid on;
plot(posL(:,2),posL(:,1),posT(:,2),posT(:,1));
legend('松组合轨迹','紧组合轨迹');
xlabel('E(m)');
ylabel('N(m)');

figure;
subplot(3,1,1);
grid on;
plot(timeL,dN);
ylabel('dN(m)');

subplot(3,1,2);
grid on;
plot(timeL,dE);
ylabel('dE(m)');

subplot(3,1,3);
grid on;
plot(timeL,dh);
ylabel('dZ(m)');
xlabel('Time(s)');

figure;
subplot(3,1,1);
grid on;
plot(timeL,dvn);
ylabel('dVn(m/s)');

subplot(3,1,2);
grid on;
plot(timeL,dve);
ylabel('dVe(m/s)');

subplot(3,1,3);
grid on;
plot(timeL,dvu);
ylabel('dVu(m/s)');
xlabel('Time(s)');

figure;
subplot(3,1,1);
grid on;
plot(timeL,droll);
ylabel('dRoll(deg)');

subplot(3,1,2);
grid on;
plot(timeL,dpitch);
ylabel('dPitch(deg)');

subplot(3,1,3);
grid on;
plot(timeL,dheading);
ylabel('dHeading(deg)');
xlabel('Time(s)');