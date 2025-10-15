
% generate map and points
load('usborder.mat','x','y');
%rng(3,'twister') % Makes stops in Maine & Florida, and is reproducible
rng('shuffle');
nStops = 50; % You can use any number, but the problem size scales as N^2

stopsLon = zeros(nStops,1); % Allocate x-coordinates of nStops
stopsLat = stopsLon; % Allocate y-coordinates
n = 1;
while (n <= nStops)
    xp = rand*1.5;
    yp = rand;
    if inpolygon(xp,yp,x,y) % Test if inside the border
        stopsLon(n) = xp;
        stopsLat(n) = yp;
        n = n+1;
    end
end

traveling_salesman(stopsLon,stopsLat);
