function [Si,Sd,h] = traveling_salesman(stopsLon,stopsLat,stopsIds)

if (nargin < 3)
    stopsIds = [];
end

% calculate distances between points
nStops = max(size(stopsLon));
idxs = nchoosek(1:nStops,2);
dist = hypot(stopsLat(idxs(:,1)) - stopsLat(idxs(:,2)), ...
             stopsLon(idxs(:,1)) - stopsLon(idxs(:,2)));
lendist = length(dist);
G = graph(idxs(:,1),idxs(:,2));

% plot initial map

h = figure;
scrsz = get(0,'ScreenSize');
set(h, 'OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)]);
subplot(2,2,1);
hold on;
hGraph = plot_base(G,stopsLon,stopsLat,nStops,stopsIds);
hold off;
title(sprintf('Data (# of stops: %d)',nStops));

%%%% solution with subtours
% constraints
tsp = optimproblem;
trips = optimvar('trips',lendist,1,'Type','integer','LowerBound',0,'UpperBound',1);
tsp.Objective = dist'*trips;

constr2trips = optimconstr(nStops,1);
for stop = 1:nStops
    whichIdxs = outedges(G,stop); % Identify trips associated with the stop
    if (stop == 1)
        constr2trips(stop) = sum(trips(whichIdxs)) == 1;
    elseif (stop == nStops)
        constr2trips(stop) = sum(trips(whichIdxs)) == 1;
    else
        constr2trips(stop) = sum(trips(whichIdxs)) == 2;
    end
end
tsp.Constraints.constr2trips = constr2trips;

% optimizatoin
opts = optimoptions('intlinprog','Display','off');
tspsol = solve(tsp,'options',opts);
tspsol.trips = logical(round(tspsol.trips));
Gsol = graph(idxs(tspsol.trips,1),idxs(tspsol.trips,2),[],numnodes(G));
tourIdxs = conncomp(Gsol);
numtours = max(tourIdxs); % Number of subtours
fprintf('# of subtours: %d\n',numtours);

% plot final map
subplot(2,2,2);
hold on;
hGraph = plot_base(G,stopsLon,stopsLat,nStops,stopsIds);
% Gsol = graph(idxs(tspsol.trips,1),idxs(tspsol.trips,2)); % Also works in most cases
highlight(hGraph,Gsol,'LineStyle','-')
hold off;
title(sprintf('Solution with Subtours (# of subtours: %d)',numtours));


%%%% solution without subtours
k = 1;
while numtours > 1 % Repeat until there is just one subtour
    
    % Add the subtour constraints
    for ii = 1:numtours
        inSubTour = (tourIdxs == ii); % Edges in current subtour
        a = all(inSubTour(idxs),2); % Complete graph indices with both ends in subtour
        constrname = "subtourconstr" + num2str(k);
        tsp.Constraints.(constrname) = sum(trips(a)) <= (nnz(inSubTour) - 1);
        k = k + 1;        
    end
    
    % Try to optimize again
    [tspsol,~,~,output] = solve(tsp,'options',opts);
    tspsol.trips = logical(round(tspsol.trips));
    Gsol = graph(idxs(tspsol.trips,1),idxs(tspsol.trips,2),[],numnodes(G));
    % Gsol = graph(idxs(tspsol.trips,1),idxs(tspsol.trips,2)); % Also works in most cases
    
    % How many subtours this time?
    tourIdxs = conncomp(Gsol);
    numtours = max(tourIdxs); % Number of subtours
    fprintf('# of subtours: %d\n',numtours);
    disp(output.absolutegap)
end

% Plot new solution
subplot(2,2,3);
hold on;
hGraph = plot_base(G,stopsLon,stopsLat,nStops,stopsIds);
highlight(hGraph,Gsol,'LineStyle','-');
hold off;
title('Solution with Subtours Eliminated');

fprintf('solution:\n');
I = [idxs(tspsol.trips,:) (1:nStops-1)'];
D = dist(tspsol.trips);
i = 1;
p = I(i,1);
q = setdiff(I(i,1:2),p);
j = setdiff(1:nStops-1,i);

Si = p;
Sd = 0;
while (p<nStops)
    num2cell([p q D(i,:)])
    Si = [Si q];
    Sd = [Sd D(i,:)];
    
    J = I(j,:);
    i = J(J(:,1)==q,3);
    if (isempty(i))
        i = J(J(:,2)==q,3);
        p = I(i,2);
    else
        p = I(i,1);
    end
    q = setdiff(I(i,1:2),p);
    j = setdiff(j,i);
end

function hGraph = plot_base(G,stopsLon,stopsLat,nStops,stopsIds)

hGraph = plot(G,'XData',stopsLon,'YData',stopsLat,'LineStyle','none','NodeLabel',{});
for i = 1:nStops
    if (isempty(stopsIds))
        text(stopsLon(i),stopsLat(i),num2str(i),'FontSize',14);
    else
        text(stopsLon(i),stopsLat(i),stopsIds{i},'FontSize',14);
    end
end
plot(stopsLon(1),stopsLat(1),'.r','markersize',20);
plot(stopsLon(end),stopsLat(end),'.r','markersize',20);
