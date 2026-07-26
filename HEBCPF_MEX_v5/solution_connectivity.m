function [G, comp, start_idx] = solution_connectivity(VBook, out_png, varargin)
%SOLUTION_CONNECTIVITY Solution connectivity diagram from the VBook record.
%
%   [G,comp] = solution_connectivity(VBook)
%   [G,comp] = solution_connectivity(VBook, 'case39_conn.png')
%   [G,comp,si] = solution_connectivity(VBook, 'case39_conn3d.png','Dim',3, 'Restarts',12, 'Start',solu, 'Zsave',Zsave)
%
% VBook(s,e) is the id of the trace that covered solution s while deforming
% equation e (0 = pair not traced). Solutions stamped with the same trace id
% lie on ONE continuation curve, hence are mutually reachable:
%   nodes  : solutions (rows of VBook), coloured by DISCOVERY ORDER (the
%            smallest trace id in each row: dark = found early, bright = late)
%   edges  : each multi-solution trace group -> star to its lowest member
%   weight : number of distinct traces linking the same pair
% A single connected component is expected for a complete single-session run
% (every solution is found on a trace through an already-known one); more than
% one component signals merged or inconsistent bookkeeping.
%
% Options (name-value):
%   'Dim'      2 or 3 (default 3). 3 uses the force3 layout.
%   'Restarts' N random layout restarts (default 12); the layout kept is the
%              "most stretched" one: max over restarts of the minimum
%              nearest-neighbour node distance (normalized scale), with mean
%              pairwise distance as tie-break. For 3D the view angle is then
%              chosen to minimize projected overlap the same way.
%   'Start'    the starting (operating-point) solution: either its column
%              INDEX in Zsave, or the solution VECTOR (e.g. solu from main.m)
%              -- the vector form requires 'Zsave' and is matched
%              antipodal-aware. Highlighted as a large BLACK marker (no label).
%   'Zsave'    solution matrix, only needed for the vector form of 'Start'.
%
% Returns the MATLAB graph G, the connected-component index comp of each
% solution, and start_idx (the matched starting-solution column). If out_png
% is given, writes the PNG plus a companion <out_png>.fig with a rotation-
% static header, so the 3D view can be rotated interactively (openfig).

opt = struct('Dim',3, 'Restarts',12, 'Start',[], 'Zsave',[]);
for k = 1:2:numel(varargin)
    opt.(varargin{k}) = varargin{k+1};
end
if nargin < 2, out_png = ''; end

%% ---- build the graph from VBook ----
n = size(VBook,1);
[r, ~, v] = find(double(VBook));
[vs, ord] = sort(v);
rs = r(ord);
bnd = [1; find(diff(vs) ~= 0) + 1; numel(vs) + 1];
src = cell(numel(bnd)-1,1); dst = cell(numel(bnd)-1,1);
n_traces = numel(bnd)-1; n_link = 0; max_group = 0;
for g = 1:n_traces
    m = rs(bnd(g):bnd(g+1)-1);
    max_group = max(max_group, numel(m));
    if numel(m) > 1
        n_link = n_link + 1;
        src{g} = repmat(m(1), numel(m)-1, 1);
        dst{g} = m(2:end);
    end
end
src = cell2mat(src); dst = cell2mat(dst);
if isempty(src)
    G = graph([], [], [], n);
else
    E = [min(src,dst) max(src,dst)];
    [Eu, ~, ic] = unique(E, 'rows');
    G = graph(Eu(:,1), Eu(:,2), accumarray(ic,1), n);
end
comp = conncomp(G)';
ncomp = max(comp);

% discovery order: a solution is first stamped by the trace that found it,
% so the smallest positive trace id in its VBook row is its discovery time.
dmin = accumarray(r, v, [n 1], @min, inf);
[~, dord] = sort(dmin);
disc_rank = zeros(n,1); disc_rank(dord) = 1:n;

%% ---- locate the starting solution ----
start_idx = [];
if ~isempty(opt.Start)
    if isscalar(opt.Start)
        start_idx = opt.Start;
    else
        if isempty(opt.Zsave)
            error('solution_connectivity:needZsave', ...
                'Vector ''Start'' needs ''Zsave'' to locate the column.');
        end
        s0 = opt.Start(:);
        d = min(vecnorm(opt.Zsave - s0, inf, 1), vecnorm(opt.Zsave + s0, inf, 1));
        [dmin, start_idx] = min(d);
        fprintf('CONN|start_match|idx=%d|dist=%.3e\n', start_idx, dmin);
    end
end

fprintf('CONN|solutions=%d|traces=%d|linking_traces=%d|edges=%d|components=%d|max_trace_group=%d\n', ...
    n, n_traces, n_link, numedges(G), ncomp, max_group);

%% ---- layout restarts: keep the most stretched embedding ----
dim = opt.Dim;
layoutname = 'force'; if dim == 3, layoutname = 'force3'; end
fl = figure('Visible','off');
bestP = []; bestScore = -inf;
for k = 1:max(1, opt.Restarts)
    rng(k, 'twister');
    ht = plot(G, 'Layout', layoutname, 'Iterations', 300, 'UseGravity', true);
    if dim == 3
        P = [ht.XData' ht.YData' ht.ZData'];
    else
        P = [ht.XData' ht.YData'];
    end
    delete(ht);
    Pn = (P - mean(P,1)) / max(eps, max(range(P,1)));   % normalized scale
    D = pdist(Pn);
    score = min(D) + 0.05*mean(D);
    if score > bestScore
        bestScore = score; bestP = P;
    end
end
close(fl);
fprintf('CONN|layout|dim=%d|restarts=%d|best_minNN=%.4f\n', dim, opt.Restarts, bestScore);

%% ---- render ----
f = figure('Visible','off', 'Position',[0 0 1500 1100], 'Color','w');
if dim == 3
    h = plot(G, 'XData',bestP(:,1), 'YData',bestP(:,2), 'ZData',bestP(:,3));
else
    h = plot(G, 'XData',bestP(:,1), 'YData',bestP(:,2));
end
h.EdgeAlpha = 0.35;
h.EdgeColor = [0.45 0.45 0.45];
if numedges(G) > 0
    w = G.Edges.Weight;
    h.LineWidth = 0.5 + 2.5*(w - min(w)) / max(1, max(w) - min(w));
end
deg = degree(G);
h.MarkerSize = 3 + 7*deg/max(1,max(deg));
% shade nodes by the order the solutions were collected: dark = early,
% bright = late (discovery rank from the smallest trace id in each row)
h.NodeCData = disc_rank;
colormap(parula(n));
cb = colorbar;
cb.Label.String = 'solution discovery order (early \rightarrow late)';

% pick the least-overlapping view for 3D
if dim == 3
    bestv = [-37.5 30]; bestm = -inf;
    Pn = (bestP - mean(bestP,1)) / max(eps, max(range(bestP,1)));
    for az = 0:30:330
        for el = [10 25 40 60]
            R1 = [cosd(az) -sind(az) 0; sind(az) cosd(az) 0; 0 0 1];
            R2 = [1 0 0; 0 cosd(el) -sind(el); 0 sind(el) cosd(el)];
            Q = Pn*R1'*R2';
            m2 = min(pdist(Q(:,[1 3])));
            if m2 > bestm, bestm = m2; bestv = [az el]; end
        end
    end
    view(bestv(1), bestv(2));
    fprintf('CONN|view|az=%g|el=%g\n', bestv(1), bestv(2));
end

% highlight the starting solution in black (no label). NodeCData shading is
% active, so overlay a marker instead of recoloring through highlight().
if ~isempty(start_idx)
    hold on
    if dim == 3
        scatter3(bestP(start_idx,1), bestP(start_idx,2), bestP(start_idx,3), ...
            120, 'k', 'filled');
    else
        scatter(bestP(start_idx,1), bestP(start_idx,2), 120, 'k', 'filled');
    end
    hold off
end

axis off; if dim == 3, axis vis3d; end

% fill the window: axes over almost the whole figure, slim colorbar at the
% right edge, tight limits. (No camzoom: a zoomed camera re-anchors
% overlays, can clip nodes after interactive rotation, and fights the
% header placement.)
axis tight
set(gca, 'Position', [0.02 0.02 0.88 0.92]);
cb.Location = 'eastoutside';
cb.Position = [0.92 0.15 0.018 0.70];

% header in its OWN invisible 2D overlay axes at the top of the figure.
% A separate axes has its own camera, so rotating the 3D plot cannot move
% it; and unlike figure annotations (which exportgraphics drops in R2022a)
% or in-axes normalized text (which the 3D camera re-anchors), plain text
% in a 2D axes is static, always visible, and always exported.
ax_main = gca;
ax_hdr = axes(f, 'Position',[0 0.955 1 0.04], 'Visible','off', ...
    'HitTest','off', 'PickableParts','none', 'Tag','conn_header');
text(ax_hdr, 0.5, 0.5, sprintf(['Solution connectivity (%dD): %d solutions, ' ...
    '%d edges, %d component(s); node size ~ degree, line width ~ #traces'], ...
    dim, n, numedges(G), ncomp), ...
    'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize', 11, ...
    'Interpreter','none');
set(f, 'CurrentAxes', ax_main);   % rotate tool targets the graph axes

if ~isempty(out_png)
    exportgraphics(f, out_png, 'Resolution', 140);
    [pp, nn] = fileparts(out_png);
    % flip visibility before saving, otherwise the .fig inherits
    % Visible='off' and opens as an invisible window in the desktop
    set(f, 'Visible', 'on');
    savefig(f, fullfile(pp, [nn '.fig']));
    fprintf('CONN|saved=%s (+.fig for interactive rotation)\n', out_png);
end
close(f);
end
