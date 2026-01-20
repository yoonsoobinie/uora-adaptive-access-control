%% UORA sweep (fairness / throughput / attempt rate / idle RU) for Origin
clear; clc;

%% ===================== Common simulation parameters =====================
sim_time  = 60;          % [s]
WIN       = 100;         % sliding window length (TF count)
N_RU      = 9;           % number of RA-RUs
numRuns   = 10;          % average over multiple runs
use_sigmoid_alpha = true; % proposed scheme (sigmoid on)

% STA sweep
N_STA_vec = 5:5:50;      % 5,10,...,50
nNsta     = numel(N_STA_vec);

% OCWmin / OCWmax configurations
cfg = [ ...
  struct('OCWmin',31,'OCWmax',511,  'name','min31max511')
% struct('OCWmin',63,'OCWmax',1023, 'name','min63max1023')
];
nCfg = numel(cfg);

% Result matrices (rows: N_STA, cols: cfg)
fairness_mat = zeros(nNsta, nCfg);
thr_mat      = zeros(nNsta, nCfg);
att_mat      = zeros(nNsta, nCfg);
idleRU_mat   = zeros(nNsta, nCfg);

%% ===================== Main sweep loop =====================
for c = 1:nCfg
    fprintf("=== Config %d/%d: %s (OCWmin=%d, OCWmax=%d)\n", ...
        c, nCfg, cfg(c).name, cfg(c).OCWmin, cfg(c).OCWmax);

    for k = 1:nNsta
        N_STA = N_STA_vec(k);

        fair_runs = zeros(1, numRuns);
        thr_runs  = zeros(1, numRuns);
        att_runs  = zeros(1, numRuns);
        idle_runs = zeros(1, numRuns);

        for r = 1:numRuns
            seed = 1000*c + 100*k + r;

            [thr_Mbps, fairness_succ, attempt_rate, avg_idle_ru] = simulate_uora_once( ...
                N_STA, N_RU, cfg(c).OCWmin, cfg(c).OCWmax, ...
                use_sigmoid_alpha, sim_time, WIN, seed);

            thr_runs(r)  = thr_Mbps;
            fair_runs(r) = fairness_succ;
            att_runs(r)  = attempt_rate;
            idle_runs(r) = avg_idle_ru;
        end

        thr_mat(k, c)      = mean(thr_runs);
        fairness_mat(k, c) = mean(fair_runs);
        att_mat(k, c)      = mean(att_runs);
        idleRU_mat(k, c)   = mean(idle_runs);

        fprintf("  N_STA=%2d → thr=%.3f Mbps, fair=%.5f, att=%.4f, idleRU=%.2f\n", ...
            N_STA, thr_mat(k,c), fairness_mat(k,c), att_mat(k,c), idleRU_mat(k,c));
    end
end

%% ===================== Origin-friendly output tables =====================
STA_col = N_STA_vec(:);
result_fairness_for_origin   = [STA_col, fairness_mat];
result_throughput_for_origin = [STA_col, thr_mat];
result_attempt_for_origin    = [STA_col, att_mat];
result_idleRU_for_origin     = [STA_col, idleRU_mat];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simulate_uora_once (includes idle RU metric)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [thr_Mbps, fairness_succ, attempt_rate, avg_idle_ru] = simulate_uora_once( ...
    N_STA, N_RU, OCWmin, OCWmax, use_sigmoid_alpha, sim_time, WIN, seed)

%% -------------------- RNG (reproducibility) --------------------
% IMPORTANT: rng must be called BEFORE any randi() initialization.
rng(seed, "twister");

%% -------------------- PHY/MAC-ish timing params --------------------
SIFS = 16e-6;
PHY_header = 40e-6;
TF = 100e-6;
BlockAck = 68e-6;

FrameSize = 2000;     % [bytes]
DataRate  = 6.67e6;   % [bps]

T_data  = PHY_header + (FrameSize*8 / DataRate);
T_total = TF + T_data + SIFS + BlockAck;

%% -------------------- Control parameters --------------------
b    = 0.1;
amin = -0.5 * N_RU;
amax =  2.0 * N_RU;

tau  = 0.8;
kmax = 3;

% sigmoid shaping
s  = 5;
p0 = 0.15;

S0 = 1/(1+exp(-s*(0 - p0)));
S1 = 1/(1+exp(-s*(1 - p0)));

%% -------------------- State initialization --------------------
t = 0;

OCWi = OCWmin * ones(1, N_STA);
OBOi = randi([0, OCWmin], 1, N_STA);
a    = zeros(1, N_STA);

fail_hist = false(N_STA, WIN);
idle_hist = false(N_STA, WIN);
filled_len = 0;

%% -------------------- Counters / stats --------------------
att_cnt        = 0;
thr_bits       = 0;
tf_cnt         = 0;
succ_per_STA   = zeros(1, N_STA);
total_idle_RUs = 0;

%% -------------------- Main simulation loop --------------------
while t < sim_time
    tf_cnt = tf_cnt + 1;

    % Decrement OBO by number of RUs per TF
    OBOi = OBOi - N_RU;

    % RU selection (0: no attempt)
    RUi = zeros(1, N_STA);
    ready = (OBOi <= a);
    if any(ready)
        RUi(ready) = randi([1, N_RU], 1, nnz(ready));
    end

    % Evaluate RU outcomes
    coll_idx = [];
    succ_idx = [];
    idle_ru_tf = 0;

    for ru = 1:N_RU
        idx = find(RUi == ru);
        k = numel(idx);

        if k == 0
            idle_ru_tf = idle_ru_tf + 1;
        elseif k == 1
            succ_idx = [succ_idx, idx]; %#ok<AGROW>
            att_cnt = att_cnt + 1;
            thr_bits = thr_bits + FrameSize*8;
            succ_per_STA(idx) = succ_per_STA(idx) + 1;
        else
            coll_idx = [coll_idx, idx]; %#ok<AGROW>
            att_cnt = att_cnt + k;
        end
    end

    total_idle_RUs = total_idle_RUs + idle_ru_tf;
    idle_idx = find(RUi == 0);

    % Sliding window bookkeeping
    y_fail = false(1, N_STA); y_fail(coll_idx) = true;
    y_idle = false(1, N_STA); y_idle(idle_idx) = true;

    col = mod(tf_cnt-1, WIN) + 1;
    fail_hist(:, col) = y_fail.';
    idle_hist(:, col) = y_idle.';

    filled_len = min(WIN, filled_len + 1);

    %% ----- Collision update (a and OCW) -----
    if ~isempty(coll_idx) && use_sigmoid_alpha
        p_fail_sw = sum(fail_hist(:, 1:filled_len), 2).' / filled_len;

        Sp = 1 ./ (1 + exp(-s * (p_fail_sw(coll_idx) - p0)));
        Sn = (Sp - S0) ./ (S1 - S0);
        Sn = max(0, min(1, Sn));

        K = 1 + (kmax - 1) .* Sn;

        a(coll_idx) = clamp(a(coll_idx) - b, amin, amax);
        OCWi(coll_idx) = min(OCWmax, floor(K .* OCWi(coll_idx) + 1));
    end

    % Resample OBO for collided STAs
    for j = 1:numel(coll_idx)
        idx = coll_idx(j);
        OBOi(idx) = randi([0, OCWi(idx)], 1, 1);
    end

    %% ----- Success update -----
    if ~isempty(succ_idx)
        OCWi(succ_idx) = OCWmin;
        a(succ_idx) = clamp(a(succ_idx) + b, amin, amax);
        OBOi(succ_idx) = randi([0, OCWmin], 1, numel(succ_idx));
    end

    %% ----- Idle update (alpha only) -----
    if ~isempty(idle_idx) && use_sigmoid_alpha
        p_idle_sw = sum(idle_hist(:, 1:filled_len), 2).' / filled_len;

        Sp_idle = 1 ./ (1 + exp(-s * (p_idle_sw(idle_idx) - p0)));
        Sn_idle = (Sp_idle - S0) ./ (S1 - S0);
        Sn_idle = max(0, min(1, Sn_idle));

        Sn_eff = Sn_idle;
        Sn_eff(Sn_eff <= tau) = 0;

        % Keep original behavior: only upper bound enforced here
        a(idle_idx) = min(amax, a(idle_idx) + b .* Sn_eff);
    end

    % Advance time
    t = t + T_total;
end

%% -------------------- Metrics --------------------
thr_Mbps     = (thr_bits / t) / 1e6;
attempt_rate = att_cnt / (tf_cnt * N_STA);
avg_idle_ru  = total_idle_RUs / tf_cnt;

% Jain's Fairness Index on success counts
if sum(succ_per_STA) == 0
    fairness_succ = 0;
else
    fairness_succ = (sum(succ_per_STA)^2) / (numel(succ_per_STA) * sum(succ_per_STA.^2));
end

end

%% ===================== helper =====================
function y = clamp(x, lo, hi)
y = min(hi, max(lo, x));
end
