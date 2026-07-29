function generate_reference_connectivity_v5(snapshot_root, reference_root)
%GENERATE_REFERENCE_CONNECTIVITY_V5 Offline case30/case57 reference graphs.
%
% This report-only utility never launches a solver. It reconstructs graph
% edges from complete VBook/Zsave snapshots from the released-v4 benchmark,
% validates those solution matrices against saved V5 scan references, and
% exports two-dimensional connectivity figures. Those raw v4 snapshots are
% not distributed with V5; pass their local directory as SNAPSHOT_ROOT.

release_root = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(snapshot_root)
    snapshot_root = fullfile(release_root,'benchmark_results_20260715', ...
        '20260715_194046','MEX_v4_20260715');
end
if nargin < 2 || isempty(reference_root)
    reference_root = fullfile(release_root,'scheduler_benchmark_v5', ...
        '20260724_3x23','references');
end
if ~isfolder(snapshot_root)
    error('HEBCPF:MissingV4ConnectivityInputs', ...
        ['Historical v4 benchmark snapshots are not distributed with V5. ' ...
         'Pass the local MEX_v4_20260715 snapshot directory as snapshot_root.']);
end
if ~isfolder(reference_root)
    error('HEBCPF:MissingV5ReferenceInputs', ...
        ['V5 reference MAT files are not present. Extract the benchmark raw ' ...
         'release asset or pass its references directory as reference_root.']);
end
addpath(fullfile(release_root,'HEBCPF_MEX_v5.2'));

cases = {'case30','case57'};
restarts = [6 4];
for k = 1:numel(cases)
    case_name = cases{k};
    old = load(fullfile(snapshot_root,[case_name '_snapshot.mat']),'VBook','Zsave');
    current = load(fullfile(reference_root,[case_name '_reference.mat']), ...
        'reference_solutions');
    distance = symmetric_set_distance(old.Zsave,current.reference_solutions);
    fprintf('REFERENCE_CONN|case=%s|solutions=%d|v4_v5_set_distance=%.6g\n', ...
        case_name,size(old.Zsave,2),distance);
    if ~isfinite(distance) || distance > 4e-7
        error('HEBCPF:ReferenceConnectivityMismatch', ...
            '%s v4/v5 solution-set distance %.6g exceeds 4e-7.', ...
            case_name,distance);
    end
    out_png = fullfile(release_root,sprintf( ...
        'connectivity_%s_v4_reference.png',case_name));
    solution_connectivity(old.VBook,out_png,'Dim',2,'Restarts',restarts(k));
end
end

function distance = symmetric_set_distance(a,b)
if ~isequal(size(a),size(b))
    distance = Inf;
    return
end
distance = 0;
for direction = 1:2
    if direction == 1
        from = a; to = b;
    else
        from = b; to = a;
    end
    for j = 1:size(from,2)
        delta_pos = max(abs(to-from(:,j)),[],1);
        delta_neg = max(abs(to+from(:,j)),[],1);
        distance = max(distance,min(min(delta_pos),min(delta_neg)));
    end
end
end
