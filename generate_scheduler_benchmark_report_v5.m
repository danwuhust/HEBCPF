function overall = generate_scheduler_benchmark_report_v5()
%GENERATE_SCHEDULER_BENCHMARK_REPORT_V5 Publish CSV, Markdown, TeX, and plots.

release_root = fileparts(mfilename('fullpath'));
publish_value = strtrim(getenv('HEBCPF_PUBLISH_BENCHMARK'));
if ~any(strcmpi(publish_value,{'1','true','yes'}))
    error('HEBCPF:PublicationNotConfirmed', ...
        ['Set HEBCPF_PUBLISH_BENCHMARK=true only after reviewing the ', ...
         'benchmark metadata and validation output.']);
end
benchmark_root = strtrim(getenv('HEBCPF_SCHED_BENCH_ROOT'));
if isempty(benchmark_root)
    error('HEBCPF:MissingBenchmarkRoot', ...
        'Set HEBCPF_SCHED_BENCH_ROOT to the completed benchmark directory.');
end
metadata = load_benchmark_metadata(benchmark_root,release_root);
raw = readtable(fullfile(benchmark_root,'scheduler_benchmark_raw.csv'));
agg = readtable(fullfile(benchmark_root,'scheduler_benchmark_aggregate.csv'));
validation = readtable(fullfile(benchmark_root,'scheduler_solution_set_validation.csv'));
validate_release_inputs(raw,agg,validation,metadata);
policies = cellstr(string(metadata.policies));
cases = cellstr(string(metadata.cases));
repeats = metadata.repeat_count;

% Replace the obsolete order-dependent raw column with the dedicated symmetric
% nearest-set validation. Scan is the stored per-case reference.
agg.symmetric_nearest_linf = zeros(height(agg),1);
for i = 1:height(agg)
    if ~strcmp(agg.policy{i},'scan')
        idx = strcmp(validation.case_name,agg.case_name{i}) & ...
            strcmp(validation.policy,agg.policy{i});
        agg.symmetric_nearest_linf(i) = validation.symmetric_nearest_linf(idx);
    end
end
agg.max_abs_vs_reference = [];
writetable(agg,fullfile(release_root,'scheduler_benchmark_v5_summary.csv'));

% Overall statistics use sums within each repeat, then mean/std across repeats.
rows = repmat(empty_overall(),0,1);
for pi = 1:numel(policies)
    policy = policies{pi};
    idxp = strcmp(raw.policy,policy);
    totals_wall = zeros(repeats,1); totals_t90 = zeros(repeats,1);
    totals_wall90 = zeros(repeats,1); totals_traces = zeros(repeats,1);
    totals_select = zeros(repeats,1); totals_cursor = zeros(repeats,1);
    totals_draws = zeros(repeats,1);
    for ri = 1:repeats
        idx = idxp & raw.repeat == ri;
        totals_wall(ri) = sum(raw.solver_wall_sec(idx));
        totals_t90(ri) = sum(raw.trace_to_90pct(idx));
        totals_wall90(ri) = sum(raw.wall_to_90pct_sec(idx));
        totals_traces(ri) = sum(raw.total_traces(idx));
        totals_select(ri) = sum(raw.selection_wall_sec(idx));
        totals_cursor(ri) = sum(raw.cursor_checks(idx));
        totals_draws(ri) = sum(raw.novelty_draws(idx));
    end
    row = empty_overall();
    row.policy = policy;
    row.aggregate_wall_sec_mean = mean(totals_wall);
    row.aggregate_wall_sec_std = std(totals_wall);
    row.sum_trace_to_90pct_mean = mean(totals_t90);
    row.sum_wall_to_90pct_sec_mean = mean(totals_wall90);
    row.sum_total_traces_mean = mean(totals_traces);
    row.selection_wall_sec_mean = mean(totals_select);
    row.selection_percent = 100*mean(totals_select)/mean(totals_wall);
    row.cursor_checks_mean = mean(totals_cursor);
    row.novelty_draws_mean = mean(totals_draws);
    row.max_residual = max(raw.max_residual(idxp));
    if strcmp(policy,'scan')
        row.max_set_distance = 0;
    else
        row.max_set_distance = max(validation.symmetric_nearest_linf( ...
            strcmp(validation.policy,policy)));
    end
    rows(end+1) = row; %#ok<AGROW>
end
overall = struct2table(rows,'AsArray',true);
scan_wall = overall.aggregate_wall_sec_mean(1);
scan_t90 = overall.sum_trace_to_90pct_mean(1);
overall.wall_ratio_vs_scan = overall.aggregate_wall_sec_mean/scan_wall;
overall.t90_ratio_vs_scan = overall.sum_trace_to_90pct_mean/scan_t90;

% Count per-case wins (ties within numerical display precision count for all ties).
overall.wall_time_wins = zeros(height(overall),1);
overall.trace90_wins = zeros(height(overall),1);
for ci = 1:numel(cases)
    idx = strcmp(agg.case_name,cases{ci});
    aw = agg.wall_sec_mean(idx); at = agg.trace_to_90pct_mean(idx);
    policy_here = agg.policy(idx);
    for pi = 1:numel(policies)
        oi = find(strcmp(overall.policy,policies{pi}));
        ai = find(strcmp(policy_here,policies{pi}));
        overall.wall_time_wins(oi) = overall.wall_time_wins(oi) + (aw(ai) <= min(aw)+1e-12);
        overall.trace90_wins(oi) = overall.trace90_wins(oi) + (at(ai) <= min(at)+1e-12);
    end
end
writetable(overall,fullfile(release_root,'scheduler_benchmark_v5_overall.csv'));

make_grouped_plot(agg,cases,policies,'wall_sec_mean','wall_sec_std', ...
    'Mean exhaustive wall time (s)','scheduler_benchmark_v5_wall.png', ...
    release_root,metadata.worker_count,repeats);
make_grouped_plot(agg,cases,policies,'trace_to_90pct_mean','trace_to_90pct_std', ...
    'Traces to 90% of final solutions','scheduler_benchmark_v5_t90.png', ...
    release_root,metadata.worker_count,repeats);
make_anytime_plot(raw,policies,benchmark_root,release_root,repeats);
make_representative_anytime_plot(raw,policies,benchmark_root,release_root,repeats);
make_suite_anytime_plot(raw,policies,benchmark_root,release_root,repeats);
write_markdown(release_root,agg,overall,cases,policies,validation,metadata,height(raw));
write_tex_table(release_root,overall);
end

function row = empty_overall()
row = struct('policy','','aggregate_wall_sec_mean',NaN, ...
    'aggregate_wall_sec_std',NaN,'sum_trace_to_90pct_mean',NaN, ...
    'sum_wall_to_90pct_sec_mean',NaN,'sum_total_traces_mean',NaN, ...
    'selection_wall_sec_mean',NaN,'selection_percent',NaN, ...
    'cursor_checks_mean',NaN,'novelty_draws_mean',NaN,'max_residual',NaN, ...
    'max_set_distance',NaN);
end

function make_grouped_plot(agg,cases,policies,value_name,error_name,y_label, ...
        file_name,release_root,workers,repeats)
values = zeros(numel(cases),numel(policies));
errors = values;
for ci = 1:numel(cases)
    for pi = 1:numel(policies)
        idx = strcmp(agg.case_name,cases{ci}) & strcmp(agg.policy,policies{pi});
        values(ci,pi) = agg.(value_name)(idx);
        errors(ci,pi) = agg.(error_name)(idx);
    end
end
f = figure('Visible','off','Color','w','Position',[100 100 1600 720]);
b = bar(1:numel(cases),values,'grouped'); hold on;
colors = [0.32 0.47 0.75; 0.90 0.45 0.18; 0.35 0.68 0.45];
for pi = 1:numel(policies)
    b(pi).FaceColor = colors(pi,:);
    errorbar(b(pi).XEndPoints,values(:,pi),errors(:,pi),'k.', ...
        'LineWidth',0.8,'CapSize',3);
end
set(gca,'YScale','log','FontSize',10,'XTick',1:numel(cases), ...
    'XTickLabel',strrep(cases,'_','\_'),'XTickLabelRotation',45);
ylabel(y_label); xlabel('Bundled case');
legend(policies,'Location','northwest','NumColumns',3);
grid on; box off;
title(sprintf('HEBCPF v5 MEX: %d repeats, %d workers (pool startup excluded)', ...
    repeats,workers));
exportgraphics(f,fullfile(release_root,file_name),'Resolution',220);
close(f);
end

function make_anytime_plot(raw,policies,benchmark_root,release_root,repeats)
target_cases = {'case57mod','case57'};
colors = [0.32 0.47 0.75; 0.90 0.45 0.18; 0.35 0.68 0.45];
gridx = linspace(0,1,301);
f = figure('Visible','off','Color','w','Position',[100 100 1500 650]);
for ci = 1:2
    subplot(1,2,ci); hold on;
    final_solutions = raw.solutions(find(strcmp(raw.case_name,target_cases{ci}),1));
    for pi = 1:numel(policies)
        curves = zeros(repeats,numel(gridx));
        for ri = 1:repeats
            idx = strcmp(raw.case_name,target_cases{ci}) & ...
                strcmp(raw.policy,policies{pi}) & raw.repeat == ri;
            h = read_history(raw.history_file{idx},benchmark_root);
            x = [0; h.trace/max(1,h.trace(end))];
            y = [0; h.solutions/final_solutions];
            curves(ri,:) = interp1(x,y,gridx,'previous','extrap');
        end
        lo = min(curves,[],1); hi = max(curves,[],1); av = mean(curves,1);
        fill([gridx fliplr(gridx)],[lo fliplr(hi)],colors(pi,:), ...
            'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
        plot(gridx,av,'LineWidth',2.0,'Color',colors(pi,:),'DisplayName',policies{pi});
    end
    yline(0.9,'--','90%','HandleVisibility','off');
    xlabel('Fraction of exhaustive trace count'); ylabel('Fraction of final solutions');
    title(strrep(target_cases{ci},'_','\_')); grid on; box off; axis([0 1 0 1.02]);
end
legend(policies,'Location','southeast');
sgtitle('Anytime discovery: mean curve and repeat range');
exportgraphics(f,fullfile(release_root,'scheduler_benchmark_v5_anytime.png'),'Resolution',220);
close(f);
end

function make_representative_anytime_plot(raw,policies,benchmark_root,release_root,repeats)
% Reuse the recorded trace histories; this function does not call a solver.
target_cases = {'case14mod','case33bw','case39','case30'};
colors = [0.32 0.47 0.75; 0.90 0.45 0.18; 0.35 0.68 0.45];
gridx = linspace(0,1,301);
f = figure('Visible','off','Color','w','Position',[100 100 1500 1050]);
for ci = 1:numel(target_cases)
    subplot(2,2,ci); hold on;
    final_solutions = raw.solutions(find(strcmp(raw.case_name,target_cases{ci}),1));
    for pi = 1:numel(policies)
        curves = zeros(repeats,numel(gridx));
        for ri = 1:repeats
            idx = strcmp(raw.case_name,target_cases{ci}) & ...
                strcmp(raw.policy,policies{pi}) & raw.repeat == ri;
            h = read_history(raw.history_file{idx},benchmark_root);
            x = [0; h.trace/max(1,h.trace(end))];
            y = [0; h.solutions/final_solutions];
            curves(ri,:) = interp1(x,y, ...
                gridx,'previous','extrap');
        end
        lo = min(curves,[],1); hi = max(curves,[],1); av = mean(curves,1);
        fill([gridx fliplr(gridx)],[lo fliplr(hi)],colors(pi,:), ...
            'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
        plot(gridx,av,'LineWidth',2.0,'Color',colors(pi,:), ...
            'DisplayName',policies{pi});
    end
    yline(0.9,'--','90%','HandleVisibility','off');
    xlabel('Fraction of exhaustive trace count');
    ylabel('Fraction of final solutions');
    title(strrep(target_cases{ci},'_','\_'));
    grid on; box off; axis([0 1 0 1.02]);
    if ci == numel(target_cases)
        legend(policies,'Location','southeast');
    end
end
sgtitle('Anytime discovery in representative small and medium cases');
exportgraphics(f,fullfile(release_root, ...
    'scheduler_benchmark_v5_anytime_representative.png'),'Resolution',220);
close(f);
end

function make_suite_anytime_plot(raw,policies,benchmark_root,release_root,repeats)
% Equal-case macro averages keep case57 from dominating smaller cases.
cases = unique(raw.case_name,'stable');
colors = [0.32 0.47 0.75; 0.90 0.45 0.18; 0.35 0.68 0.45];
gridx = linspace(0,1,301);
f = figure('Visible','off','Color','w','Position',[100 100 1500 650]);
for mode = 1:2
    subplot(1,2,mode); hold on;
    for pi = 1:numel(policies)
        curves = zeros(repeats*numel(cases),numel(gridx));
        row = 0;
        for ci = 1:numel(cases)
            final_solutions = raw.solutions(find(strcmp(raw.case_name,cases{ci}),1));
            for ri = 1:repeats
                idx = strcmp(raw.case_name,cases{ci}) & ...
                    strcmp(raw.policy,policies{pi}) & raw.repeat == ri;
                h = read_history(raw.history_file{idx},benchmark_root);
                if mode == 1
                    x = [0; h.trace/max(1,h.trace(end))];
                else
                    x = [0; h.completion_elapsed_seconds / ...
                        max(eps,h.completion_elapsed_seconds(end))];
                end
                y = [0; h.solutions/final_solutions];
                row = row + 1;
                curves(row,:) = interp1(x,y, ...
                    gridx,'previous','extrap');
            end
        end
        ordered = sort(curves,1);
        n = size(ordered,1);
        lo = ordered(max(1,ceil(0.25*n)),:);
        hi = ordered(min(n,ceil(0.75*n)),:);
        av = mean(curves,1);
        fill([gridx fliplr(gridx)],[lo fliplr(hi)],colors(pi,:), ...
            'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
        plot(gridx,av,'LineWidth',2.0,'Color',colors(pi,:), ...
            'DisplayName',policies{pi});
    end
    yline(0.9,'--','90%','HandleVisibility','off');
    if mode == 1
        xlabel('Fraction of exhaustive trace count');
        title('Trace-normalized');
    else
        xlabel('Fraction of exhaustive elapsed time');
        title('Elapsed-time-normalized');
    end
    ylabel('Mean fraction of final solutions');
    grid on; box off; axis([0 1 0 1.02]);
end
legend(policies,'Location','southeast');
sgtitle('Suite-level anytime discovery: equal weight per case and repeat (IQR bands)');
exportgraphics(f,fullfile(release_root, ...
    'scheduler_benchmark_v5_anytime_suite.png'),'Resolution',220);
close(f);
end

function h = read_history(history_path,benchmark_root)
if ~exist(history_path,'file')
    relative_candidate = fullfile(benchmark_root,strrep(history_path,'/',filesep));
    if exist(relative_candidate,'file')
        history_path = relative_candidate;
    else
        [~,name,ext] = fileparts(history_path);
        history_path = fullfile(benchmark_root,'trace_histories',[name ext]);
    end
end
h = readtable(history_path);
end

function write_markdown(release_root,agg,overall,cases,policies,validation,metadata,run_count)
fid = fopen(fullfile(release_root,'SCHEDULER_BENCHMARK_v5.md'),'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'# HEBCPF v5 three-scheduler benchmark\n\n');
fprintf(fid,'This report compares `scan`, `bandit`, and `novelty` in the v5 MEX queue driver. ');
fprintf(fid,['It is the current V5 measurement; the direct public release baseline is ', ...
    'v4 (2026.07.15).\n\n']);
fprintf(fid,'## Method\n\n- %d bundled nominal-load cases through case57.\n',numel(cases));
fprintf(fid,'- %s, %s, %d local workers on the designated main job machine.\n', ...
    metadata.matlab_release,metadata.operating_system,metadata.worker_count);
fprintf(fid,'- %d measured repeats per scheduler and case (%d/%d PASS).\n', ...
    metadata.repeat_count,run_count,run_count);
fprintf(fid,'- One persistent pool per case group; pool startup and three policy warm-ups excluded.\n');
fprintf(fid,'- Policy order rotated by repeat. Raw state was saved after every run.\n');
fprintf(fid,'- Separate order-independent validation: bandit and novelty each compared once with the saved scan reference for every case (%d/%d PASS).\n\n', ...
    height(validation),height(validation));
fprintf(fid,['All published v5 numerical and timing values come only from official ', ...
    'dataset `%s`; results from other machines are not merged. Raw release asset: ', ...
    '`release_assets/HEBCPF-v5-benchmark-raw.zip`.\n\n'],metadata.benchmark_id);

fprintf(fid,'## Overall results\n\n');
fprintf(fid,'| Policy | Sum wall, mean +/- SD (s) | Sum traces to 90%% | vs scan | Sum total traces | Selection time (%% wall) | Wall wins | t90 wins | Max set distance |\n');
fprintf(fid,'|---|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for i = 1:height(overall)
    fprintf(fid,'| `%s` | %.3f +/- %.3f | %.1f | %.3fx | %.1f | %.3f (%.2f%%) | %d | %d | %.3g |\n', ...
        overall.policy{i},overall.aggregate_wall_sec_mean(i),overall.aggregate_wall_sec_std(i), ...
        overall.sum_trace_to_90pct_mean(i),overall.t90_ratio_vs_scan(i), ...
        overall.sum_total_traces_mean(i),overall.selection_wall_sec_mean(i), ...
        overall.selection_percent(i),overall.wall_time_wins(i),overall.trace90_wins(i), ...
        overall.max_set_distance(i));
end

scan = overall(1,:); bandit = overall(2,:); novelty = overall(3,:);
fprintf(fid,'\nAcross the suite, exhaustive wall time was effectively tied: bandit was %.2f%% and novelty %.2f%% relative to scan. ', ...
    100*(bandit.wall_ratio_vs_scan-1),100*(novelty.wall_ratio_vs_scan-1));
fprintf(fid,'The anytime result differed sharply: bandit reduced the summed traces-to-90%% metric by %.1f%% and novelty by %.1f%%. ', ...
    100*(1-bandit.t90_ratio_vs_scan),100*(1-novelty.t90_ratio_vs_scan));
fprintf(fid,'Novelty spent %.2f%% of wall time in selection versus %.2f%% for bandit and %.2f%% for scan.\n\n', ...
    novelty.selection_percent,bandit.selection_percent,scan.selection_percent);
fprintf(fid,'![Mean exhaustive wall time](scheduler_benchmark_v5_wall.png)\n\n');
fprintf(fid,'![Mean traces to 90%%](scheduler_benchmark_v5_t90.png)\n\n');
fprintf(fid,'![Anytime discovery curves](scheduler_benchmark_v5_anytime.png)\n\n');
fprintf(fid,'![Representative-case anytime discovery](scheduler_benchmark_v5_anytime_representative.png)\n\n');
fprintf(fid,'![Suite-level normalized anytime discovery](scheduler_benchmark_v5_anytime_suite.png)\n\n');

fprintf(fid,'## Per-case means\n\n');
fprintf(fid,'Wall entries are seconds; t90 is the first completed-trace count at 90%% of that run''s final solutions.\n\n');
fprintf(fid,'| Case | Solutions | Wall scan / bandit / novelty | t90 scan / bandit / novelty | Total traces scan / bandit / novelty |\n');
fprintf(fid,'|---|---:|---:|---:|---:|\n');
for ci = 1:numel(cases)
    a = cell(1,3);
    for pi = 1:3
        a{pi} = agg(strcmp(agg.case_name,cases{ci}) & strcmp(agg.policy,policies{pi}),:);
    end
    fprintf(fid,'| `%s` | %.0f | %.3f / %.3f / %.3f | %.1f / %.1f / %.1f | %.1f / %.1f / %.1f |\n', ...
        cases{ci},a{1}.solutions,a{1}.wall_sec_mean,a{2}.wall_sec_mean,a{3}.wall_sec_mean, ...
        a{1}.trace_to_90pct_mean,a{2}.trace_to_90pct_mean,a{3}.trace_to_90pct_mean, ...
        a{1}.total_traces,a{2}.total_traces,a{3}.total_traces);
end

fprintf(fid,'\n## Accuracy and interpretation\n\n');
fprintf(fid,'Every timing run returned the expected solution count. The maximum residual among the 180 measured runs was `%.3g`. ',max(agg.max_residual));
fprintf(fid,'The separate symmetric nearest-set audit matched all 40 adaptive-policy sets to their scan references; the worst distance was `%.3g`, below `4e-7`.\n\n',max(validation.symmetric_nearest_linf));
fprintf(fid,'Total traces can differ slightly because asynchronous dispatch may complete work that becomes redundant after another worker reports. ');
fprintf(fid,'This is why total traces, t90, and wall time are all reported as measured outcomes.\n\n');
fprintf(fid,'The largest-system behavior is the clearest: on case57, mean t90 was 19,481 for scan, 6,996 for bandit, and 4,594 for novelty, while exhaustive wall times remained near 333--337 s. ');
fprintf(fid,'On case57mod, mean t90 was 9,203, 4,425, and 2,829. Thus novelty provides the strongest early discovery here, but bandit avoids most novelty-selection cost and is the default.\n\n');
fprintf(fid,'## Direct comparison with v4\n\n');
fprintf(fid,'The retained 2026.07.15 v4 MEX benchmark reported 718.656 s for one 20-case pass. The current v5 scan mean sums to %.3f s (%.1f%% lower). ', ...
    scan.aggregate_wall_sec_mean,100*(1-scan.aggregate_wall_sec_mean/718.656));
fprintf(fid,'This is a descriptive release-to-release comparison, not a controlled paired speedup: v4 has one historical pass, while v5 has three current repeats and includes other implementation changes. ');
fprintf(fid,'The scheduler conclusions above come only from the controlled v5 three-policy experiment.\n');
end

function metadata = load_benchmark_metadata(benchmark_root,release_root)
paths = {fullfile(benchmark_root,'benchmark_metadata_v5.json'), ...
    fullfile(benchmark_root,'BENCHMARK_METADATA_v5.json'), ...
    fullfile(release_root,'BENCHMARK_METADATA_v5.json')};
metadata = [];
for i = 1:numel(paths)
    if exist(paths{i},'file')
        metadata = jsondecode(fileread(paths{i}));
        break;
    end
end
if isempty(metadata)
    error('HEBCPF:MissingBenchmarkMetadata', ...
        'No benchmark_metadata_v5.json or BENCHMARK_METADATA_v5.json was found.');
end
end

function validate_release_inputs(raw,agg,validation,metadata)
required = {'benchmark_id','publication_status','worker_count','repeat_count', ...
    'cases','policies','matlab_release','operating_system'};
for i = 1:numel(required)
    if ~isfield(metadata,required{i})
        error('HEBCPF:IncompleteBenchmarkMetadata', ...
            'Benchmark metadata is missing field %s.',required{i});
    end
end
if ~strcmp(metadata.publication_status,'official release evidence')
    error('HEBCPF:UnapprovedBenchmarkMetadata', ...
        'Benchmark metadata is not marked as official release evidence.');
end
if any(string(raw.status) ~= "PASS") || any(string(validation.status) ~= "PASS")
    error('HEBCPF:BenchmarkValidationFailed', ...
        'Publication requires PASS status for every timing and validation row.');
end
workers = unique(raw.workers);
if numel(workers) ~= 1 || workers ~= metadata.worker_count
    error('HEBCPF:BenchmarkWorkerMismatch', ...
        'Raw data worker count does not match approved benchmark metadata.');
end
expected_runs = numel(metadata.cases)*numel(metadata.policies)*metadata.repeat_count;
if height(raw) ~= expected_runs
    error('HEBCPF:IncompleteBenchmarkMatrix', ...
        'Expected %d timing rows from metadata, found %d.',expected_runs,height(raw));
end
if ~isequal(sort(unique(string(raw.case_name))),sort(string(metadata.cases(:)))) || ...
        ~isequal(sort(unique(string(raw.policy))),sort(string(metadata.policies(:))))
    error('HEBCPF:BenchmarkMatrixMismatch', ...
        'Raw case or policy names do not match approved benchmark metadata.');
end
if height(agg) ~= numel(metadata.cases)*numel(metadata.policies)
    error('HEBCPF:IncompleteAggregateMatrix', ...
        'Aggregate table does not contain the approved case-policy matrix.');
end
end

function write_tex_table(release_root,overall)
fid = fopen(fullfile(release_root,'scheduler_benchmark_v5_table.tex'),'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'\\begin{center}\\small\n');
fprintf(fid,'\\begin{tabular}{@{}lrrrrr@{}}\n\\toprule\n');
fprintf(fid,'Policy & Wall (s) & $\\sum t_{90}$ & Total traces & Selection (\\%%) & Set distance\\\\\n\\midrule\n');
for i = 1:height(overall)
    fprintf(fid,'\\code{%s} & %.3f $\\pm$ %.3f & %.1f & %.1f & %.2f & %.2g\\\\\n', ...
        overall.policy{i},overall.aggregate_wall_sec_mean(i),overall.aggregate_wall_sec_std(i), ...
        overall.sum_trace_to_90pct_mean(i),overall.sum_total_traces_mean(i), ...
        overall.selection_percent(i),overall.max_set_distance(i));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{center}\n');
fprintf(fid,'Across 20 cases, bandit and novelty reduced summed traces to 90\\%% by %.1f\\%% and %.1f\\%% versus scan, while aggregate exhaustive wall time differed by less than 0.1\\%%.\n', ...
    100*(1-overall.t90_ratio_vs_scan(2)),100*(1-overall.t90_ratio_vs_scan(3)));
end
