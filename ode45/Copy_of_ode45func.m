function ydot = ode45func(t,y)
L=-5*[-1.5,1.5,0,0,0,0;2,-2,0,0,0,0;0.9,0,-2.8,0,1.9,0;0,1.2,0,-2.5,0,1.3;0,0,1.4,1.8,-3.2,0;0,0,0,0,0.7,-0.7];
ydot=-L*y+1/5*abs(1/4*cos(2*pi*t)); %΢�ַ���