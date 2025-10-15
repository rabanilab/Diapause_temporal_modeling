function h = big_figure()

h = figure;
scrsz = get(0,'ScreenSize'); %[left, bottom, width, height]
set(h, 'OuterPosition',[1 scrsz(4) 0.5*scrsz(3) 0.5*scrsz(4)]);
