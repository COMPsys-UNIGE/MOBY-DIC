[x y] = meshgrid(linspace(-10,8,100),linspace(-10,8,100));
p = [x(:) y(:)];

global usemex

usemex = 0;

tic
[u allu] = ctrl.eval(p);
toc

usemex = 1;

tic
[umex allumex] = ctrl.eval(p);
toc

figure
hold on
plot(allu{1},'.')
plot(allumex{1},'rx')

figure
hold on
plot(allu{2},'.')
plot(allumex{2},'rx')

figure
hold on
plot(allu{3},'.')
plot(allumex{3},'rx')
