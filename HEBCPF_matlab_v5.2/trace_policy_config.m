function [policy,cfg] = trace_policy_config(policy_in)
%TRACE_POLICY_CONFIG Normalize HEBCPF queue-scheduler policy and parameters.
% Canonical policies are scan, bandit (default), and novelty. The v5 name
% diverse remains a deprecated alias for bandit so existing scripts retain
% their original equation-gain behavior.

if nargin < 1 || isempty(policy_in)
    policy = 'bandit';
else
    if isstring(policy_in) && isscalar(policy_in)
        policy_in = char(policy_in);
    end
    if ~ischar(policy_in)
        error('HEBCPF:InvalidPolicyType', ...
            'HEBCPOLICY must be a character vector or scalar string.');
    end
    policy = lower(strtrim(policy_in));
end

persistent diverse_warning_shown
if strcmp(policy,'diverse')
    if isempty(diverse_warning_shown) || ~diverse_warning_shown
        warning('HEBCPF:DeprecatedDiversePolicy', ...
            ['HEBCPOLICY=''diverse'' is the v5 name for equation-gain ', ...
             'ordering and now maps to ''bandit''. Use ''novelty'' for ', ...
             'bounded novelty-based start selection.']);
        diverse_warning_shown = true;
    end
    policy = 'bandit';
end

valid = {'scan','bandit','novelty'};
if ~any(strcmp(policy,valid))
    error('HEBCPF:UnknownPolicy', ...
        'Unknown HEBCPOLICY ''%s''. Use scan, bandit, or novelty.', policy);
end

cfg = struct();
cfg.default_policy = 'bandit';
cfg.novelty_candidate_cap = 96;
cfg.novelty_reference_cap = 96;
cfg.novelty_draw_budget = 400;
cfg.novelty_projection_dim = 8;
cfg.novelty_projection_seed = 11;
cfg.novelty_sampling_seed = 29;
cfg.bandit_discount_old = 0.7;
cfg.bandit_discount_new = 0.3;
cfg.bandit_optimistic_gain = 5;
cfg.stall_sweeps = 2;
end
