%
% Included:
%   1) Internally-monitored experiment
%      - Subjective Timing Model
%      - Granule-Cell Basis Model
%   2) Externally-cued experiment
%      - Subjective Timing Model with Weber fraction fixed at 0.001
%      - Granule-Cell Basis Model
%
% statistics calculated:
%   - group-level R2 across all five test delays
%   - participant-level RMSE across the three untrained delays,
%     compared with a Wilcoxon signed-rank test
%
% Inputs:
%   data/internallymonitored_model.csv
%   data/externallycued_model.csv
%
% results are printed out in the command window
clear; clc; close all;
rng(1);

BASE_DIR = fileparts(mfilename('fullpath'));

settings.internallymonitored.IN_CSV = ...
    fullfile(BASE_DIR, 'data', 'internallymonitored_model.csv');
settings.externallycued.IN_CSV = ...
    fullfile(BASE_DIR, 'data', 'externallycued_model.csv');



% Inputs are already participant-cleaned.
settings.internallymonitored.exclude_ids = [];
settings.externallycued.exclude_ids = [];



% Preprocessing
settings.rotAmp = 15;
settings.baselineMode = "phase2_only";
settings.desired_rot500  =  1;
settings.desired_rot1500 = -1;
settings.trainDelays = [500 1500];
settings.testDelays  = [500 750 1000 1250 1500];

settings.flipVars = {'handAng','rotDir'};

settings.learning_phase = 21:420;
settings.phase3Mode = "mean_center";
settings.freezeDuringPhase3 = true;

% Optimization
settings.iterations_timing = 20;
settings.iterations_model1 = 100;
settings.iterations_model2 = 100;

% Subjective Timing Model
settings.intervals_sec = [0.5 1.5];
settings.externallycuedWeber = 0.001;

% Granule-Cell Basis Model
settings.granule.centers_ms  = 250:100:1750;
settings.granule.baseSigma   = 100;
settings.granule.widthSlope  = 0.35;
settings.granule.decayTau    = 10000;
settings.granule.timingNoise = 150;

% Test-phase scale factor
settings.testScale = 0.5;

%% PREPROCESS DATA

fprintf('PREPROCESSING DATA\n');

fprintf('\n--- Internally-monitored ---\n');
[internallymonitored, ~, ~] = ...
    preprocessData(settings, settings.internallymonitored.IN_CSV, ...
    settings.internallymonitored.exclude_ids);

assert(isfield(internallymonitored,'rt'), ...
    'Internally-monitored model file must contain rtGoSig.');

fprintf('\n--- Externally-cued ---\n');
[externallycued, ~, ~] = ...
    preprocessData(settings, settings.externallycued.IN_CSV, ...
    settings.externallycued.exclude_ids);

[n_internallymonitored, ~] = size(internallymonitored.ha);
[n_externallycued, ~] = size(externallycued.ha);

fprintf('\nParticipants: internally-monitored = %d; externally-cued = %d\n', ...
    n_internallymonitored, n_externallycued);

basisCache = makeGranuleBasisCache(settings);



    % Internally-monitored: Subjective Timing Model
    fprintf('FITTING INTERNALLY-MONITORED SUBJECTIVE TIMING MODEL\n');

    fits_internallymonitored_timing = initTimingFits(n_internallymonitored, settings);
    fits_internallymonitored_subjectivetiming = initModel1Fits(n_internallymonitored);

    for si = 1:n_internallymonitored
        fprintf('\nSubject %d/%d | id %d\n', ...
            si, n_internallymonitored, internallymonitored.id(si));

        fits_internallymonitored_timing = fitTimingModel( ...
            fits_internallymonitored_timing, si, ...
            internallymonitored.timing(si,:), ...
            internallymonitored.rt(si,:) ./ 1000, settings);

        fits_internallymonitored_subjectivetiming = fitModel1( ...
            fits_internallymonitored_subjectivetiming, si, ...
            internallymonitored.ha(si,:), ...
            internallymonitored.timing(si,:), ...
            internallymonitored.rotation(si,:), settings);
    end

    % Internally-monitored: Granule-Cell Basis Model
    fprintf('FITTING INTERNALLY-MONITORED GRANULE-CELL BASIS MODEL\n');

    fits_internallymonitored_granulecell = initModel2Fits(n_internallymonitored);

    for si = 1:n_internallymonitored
        fprintf('\nSubject %d/%d | id %d\n', ...
            si, n_internallymonitored, internallymonitored.id(si));

        fits_internallymonitored_granulecell = fitModel2( ...
            fits_internallymonitored_granulecell, si, ...
            internallymonitored.ha(si,:), ...
            internallymonitored.timing(si,:), ...
            internallymonitored.rotation(si,:), ...
            basisCache, settings);
    end

    % Externally-cued: Granule-Cell Basis Model

    fprintf('FITTING EXTERNALLY-CUED GRANULE-CELL BASIS MODEL\n');

    fits_externallycued_granulecell = initModel2Fits(n_externallycued);

    for si = 1:n_externallycued
        fprintf('\nSubject %d/%d | id %d\n', ...
            si, n_externallycued, externallycued.id(si));

        fits_externallycued_granulecell = fitModel2( ...
            fits_externallycued_granulecell, si, ...
            externallycued.ha(si,:), ...
            externallycued.timing(si,:), ...
            externallycued.rotation(si,:), ...
            basisCache, settings);
    end

    % Externally-cued: Subjective Timing Model
    % Timing uncertainty is fixed at Weber = 0.001.
    fprintf('FITTING EXTERNALLY-CUED SUBJECTIVE TIMING MODEL\n');
    fprintf('Weber fraction fixed at %.3f\n', settings.externallycuedWeber);

    fits_externallycued_timing = initTimingFits(n_externallycued, settings);
    fits_externallycued_timing.weber(:) = settings.externallycuedWeber;
    fits_externallycued_subjectivetiming = initModel1Fits(n_externallycued);

    for si = 1:n_externallycued
        fprintf('\nSubject %d/%d | id %d\n', ...
            si, n_externallycued, externallycued.id(si));

        fits_externallycued_subjectivetiming = fitModel1( ...
            fits_externallycued_subjectivetiming, si, ...
            externallycued.ha(si,:), ...
            externallycued.timing(si,:), ...
            externallycued.rotation(si,:), settings);
    end


%% FITTED PARAMETERS

fprintf('FITTED PARAMETER SUMMARY\n');

print_parameter_summary( ...
    'INTERNALLY-MONITORED', ...
    fits_internallymonitored_timing, ...
    fits_internallymonitored_subjectivetiming, ...
    fits_internallymonitored_granulecell, true);

print_parameter_summary( ...
    'EXTERNALLY-CUED', ...
    fits_externallycued_timing, ...
    fits_externallycued_subjectivetiming, ...
    fits_externallycued_granulecell, false);



%% TEST-PHASE GENERALIZATION
% Models were fit only to training data.

fprintf('TEST-PHASE GENERALIZATION (fixed s = %.2f)\n', settings.testScale);

t_ms = 0:1:2000;

% Internally-monitored Subjective Timing Model
exp_state_internallymonitored_subjectivetiming = ...
    nan(n_internallymonitored, size(internallymonitored.ha,2), numel(t_ms));

for si = 1:n_internallymonitored
    p = [fits_internallymonitored_subjectivetiming.A(si), ...
         fits_internallymonitored_subjectivetiming.B(si), ...
         fits_internallymonitored_subjectivetiming.Ac(si), ...
         fits_internallymonitored_subjectivetiming.Bc(si)];

    [~, exp_state_internallymonitored_subjectivetiming(si,:,:)] = ...
        simulateModel1Full( ...
            p, internallymonitored.timing(si,:), ...
            internallymonitored.rotation(si,:), ...
            internallymonitored.phase(si,:), ...
            fits_internallymonitored_timing.weber(si), ...
            t_ms, settings);
end

[data_gen_internallymonitored, pred_internallymonitored_subjectivetiming] = ...
    computePhase3FromExpected( ...
        internallymonitored, ...
        exp_state_internallymonitored_subjectivetiming, ...
        t_ms, settings.testDelays, settings.testScale, settings.phase3Mode);

% Internally-monitored Granule-Cell Basis Model
pred_internallymonitored_granulecell = computePhase3GranuleDirect( ...
    internallymonitored, fits_internallymonitored_granulecell, ...
    basisCache, settings, settings.testScale, settings.phase3Mode);

% Externally-cued Subjective Timing Model
exp_state_externallycued_subjectivetiming = ...
    nan(n_externallycued, size(externallycued.ha,2), numel(t_ms));

for si = 1:n_externallycued
    p = [fits_externallycued_subjectivetiming.A(si), ...
         fits_externallycued_subjectivetiming.B(si), ...
         fits_externallycued_subjectivetiming.Ac(si), ...
         fits_externallycued_subjectivetiming.Bc(si)];

    [~, exp_state_externallycued_subjectivetiming(si,:,:)] = ...
        simulateModel1Full( ...
            p, externallycued.timing(si,:), ...
            externallycued.rotation(si,:), ...
            externallycued.phase(si,:), ...
            settings.externallycuedWeber, ...
            t_ms, settings);
end

[data_gen_externallycued, pred_externallycued_subjectivetiming] = ...
    computePhase3FromExpected( ...
        externallycued, ...
        exp_state_externallycued_subjectivetiming, ...
        t_ms, settings.testDelays, settings.testScale, settings.phase3Mode);

% Externally-cued Granule-Cell Basis Model
pred_externallycued_granulecell = computePhase3GranuleDirect( ...
    externallycued, fits_externallycued_granulecell, ...
    basisCache, settings, settings.testScale, settings.phase3Mode);

%% RESULTS

stats_internallymonitored_subjectivetiming = summarizeGenFit( ...
    data_gen_internallymonitored, pred_internallymonitored_subjectivetiming);

stats_internallymonitored_granulecell = summarizeGenFit( ...
    data_gen_internallymonitored, pred_internallymonitored_granulecell);

stats_externallycued_subjectivetiming = summarizeGenFit( ...
    data_gen_externallycued, pred_externallycued_subjectivetiming);

stats_externallycued_granulecell = summarizeGenFit( ...
    data_gen_externallycued, pred_externallycued_granulecell);

untrained = ismember(settings.testDelays, [750 1000 1250]);

[rmse_im_subjectivetiming, rmse_im_granulecell, W_im, p_im] = ...
    compare_untrained_rmse( ...
        data_gen_internallymonitored(:,untrained), ...
        pred_internallymonitored_subjectivetiming(:,untrained), ...
        pred_internallymonitored_granulecell(:,untrained));

[rmse_ec_subjectivetiming, rmse_ec_granulecell, W_ec, p_ec] = ...
    compare_untrained_rmse( ...
        data_gen_externallycued(:,untrained), ...
        pred_externallycued_subjectivetiming(:,untrained), ...
        pred_externallycued_granulecell(:,untrained));

fprintf('\nINTERNALLY-MONITORED\n');
fprintf('  Weber fraction: %.3f +/- %.3f (SD)\n', ...
    mean(fits_internallymonitored_timing.weber,'omitnan'), ...
    std(fits_internallymonitored_timing.weber,0,'omitnan'));
fprintf('  Subjective Timing Model group R2: %.3f\n', ...
    stats_internallymonitored_subjectivetiming.R2_group);
fprintf('  Granule-Cell Basis Model group R2: %.3f\n', ...
    stats_internallymonitored_granulecell.R2_group);
fprintf('  Untrained-delay RMSE, Subjective Timing Model: %.3f\n', ...
    mean(rmse_im_subjectivetiming,'omitnan'));
fprintf('  Untrained-delay RMSE, Granule-Cell Basis Model: %.3f\n', ...
    mean(rmse_im_granulecell,'omitnan'));
fprintf('  Wilcoxon signed-rank: W = %.0f, p = %.6f\n', W_im, p_im);

fprintf('\nEXTERNALLY-CUED\n');
fprintf('  Weber fraction: %.3f (fixed)\n', settings.externallycuedWeber);
fprintf('  Subjective Timing Model group R2: %.3f\n', ...
    stats_externallycued_subjectivetiming.R2_group);
fprintf('  Granule-Cell Basis Model group R2: %.3f\n', ...
    stats_externallycued_granulecell.R2_group);
fprintf('  Untrained-delay RMSE, Subjective Timing Model: %.3f\n', ...
    mean(rmse_ec_subjectivetiming,'omitnan'));
fprintf('  Untrained-delay RMSE, Granule-Cell Basis Model: %.3f\n', ...
    mean(rmse_ec_granulecell,'omitnan'));
fprintf('  Wilcoxon signed-rank: W = %.0f, p = %.6f\n', W_ec, p_ec);





fprintf('\nDone.\n');

%% LOCAL FUNCTIONS

function se = standard_error(x)
    x = x(isfinite(x));
    se = std(x) / sqrt(numel(x));
end

function print_parameter_summary(label, timingFit, subjectiveTimingFit, granuleCellFit, printTiming)

    fprintf('\n--- %s ---\n', label);

    if printTiming
        fprintf('Timing model:\n');
        fprintf('  Weber    = %.3f +/- %.3f\n', ...
            mean(timingFit.weber,'omitnan'), standard_error(timingFit.weber));
        fprintf('  RT floor = %.3f +/- %.3f\n', ...
            mean(timingFit.rtmin,'omitnan'), standard_error(timingFit.rtmin));
        fprintf('  RT scale = %.3f +/- %.3f\n', ...
            mean(timingFit.rtscale,'omitnan'), standard_error(timingFit.rtscale));
    end

    fprintf('Subjective Timing Model:\n');
    fprintf('  A  = %.3f +/- %.3f\n', ...
        mean(subjectiveTimingFit.A,'omitnan'), standard_error(subjectiveTimingFit.A));
    fprintf('  B  = %.3f +/- %.3f\n', ...
        mean(subjectiveTimingFit.B,'omitnan'), standard_error(subjectiveTimingFit.B));
    fprintf('  Ac = %.3f +/- %.3f\n', ...
        mean(subjectiveTimingFit.Ac,'omitnan'), standard_error(subjectiveTimingFit.Ac));
    fprintf('  Bc = %.3f +/- %.3f\n', ...
        mean(subjectiveTimingFit.Bc,'omitnan'), standard_error(subjectiveTimingFit.Bc));

    fprintf('Granule-Cell Basis Model:\n');
    fprintf('  A  = %.3f +/- %.3f\n', ...
        mean(granuleCellFit.A,'omitnan'), standard_error(granuleCellFit.A));
    fprintf('  B  = %.3f +/- %.3f\n', ...
        mean(granuleCellFit.B,'omitnan'), standard_error(granuleCellFit.B));
    fprintf('  Aw = %.3f +/- %.3f\n', ...
        mean(granuleCellFit.Aw,'omitnan'), standard_error(granuleCellFit.Aw));
    fprintf('  Bw = %.3f +/- %.3f\n', ...
        mean(granuleCellFit.Bw,'omitnan'), standard_error(granuleCellFit.Bw));
end

function [rmse_subjectivetiming, rmse_granulecell, W, p] = ...
    compare_untrained_rmse(data, pred_subjectivetiming, pred_granulecell)

    rmse_subjectivetiming = ...
        sqrt(mean((data-pred_subjectivetiming).^2, 2, 'omitnan'));

    rmse_granulecell = ...
        sqrt(mean((data-pred_granulecell).^2, 2, 'omitnan'));

    ok = isfinite(rmse_subjectivetiming) & isfinite(rmse_granulecell);

    [p,~,statsW] = signrank( ...
        rmse_subjectivetiming(ok), rmse_granulecell(ok));

    W = statsW.signedrank;
end



%% PREPROCESSING FUNCTIONS

function [data, T, prepInfo] = preprocessData(settings, in_csv, exclude_ids)

    T = readtable(in_csv);

    mustHave = {'id','subTrial','phase','goDelay','rotDir','handAng'};
    for k = 1:numel(mustHave)
        assert(ismember(mustHave{k}, T.Properties.VariableNames), ...
            'Missing column: %s in %s', mustHave{k}, in_csv);
    end

    if ~isempty(exclude_ids)
        T = T(~ismember(T.id, exclude_ids), :);
    end

    T = sortrows(T, {'id','subTrial'});

    participants = unique(T.id, 'sorted');
    nSubs = numel(participants);

    fprintf('Found %d participants\n', nSubs);

    % Baseline subtraction
    fprintf('Applying baseline mode: %s\n', settings.baselineMode);

    if settings.baselineMode ~= "none"
        for si = 1:nSubs
            pid = participants(si);

            for d = settings.trainDelays
                rows_p1 = T.id==pid & T.phase==1 & T.goDelay==d;
                p1_mean = mean(T.handAng(rows_p1), 'omitnan');

                if ~isfinite(p1_mean)
                    error('Subject %d has no valid Phase 1 baseline for %d ms.', ...
                        pid, d);
                end

                switch settings.baselineMode
                    case "phase2_only"
                        rows_fix = ...
                            T.id==pid & T.phase==2 & T.goDelay==d;

                    otherwise
                        error('Unknown baselineMode: %s', settings.baselineMode);
                end

                T.handAng(rows_fix) = T.handAng(rows_fix) - p1_mean;
            end
        end
    end

    % Detect whether participant-level sign alignment is still required.
    auto_flip_ids = [];
    mapping_summary = table();

    for si = 1:nSubs
        pid = participants(si);

        r500 = T.id==pid & T.phase==2 & ...
            T.goDelay==500 & T.rotDir~=0 & ~isnan(T.rotDir);

        r1500 = T.id==pid & T.phase==2 & ...
            T.goDelay==1500 & T.rotDir~=0 & ~isnan(T.rotDir);

        v500 = cleanUnique(T.rotDir(r500));
        v1500 = cleanUnique(T.rotDir(r1500));

        rot500 = NaN;
        rot1500 = NaN;
        flipFlag = false;
        mappingLabel = "unclear";

        if numel(v500)==1 && numel(v1500)==1
            rot500 = v500(1);
            rot1500 = v1500(1);

            if rot500==settings.desired_rot500 && ...
                    rot1500==settings.desired_rot1500

                mappingLabel = "aligned";

            elseif rot500==-settings.desired_rot500 && ...
                    rot1500==-settings.desired_rot1500

                flipFlag = true;
                mappingLabel = "reversed -- will flip";
                auto_flip_ids(end+1) = pid; %#ok<AGROW>

            else
                error('Unexpected rotation mapping for subject %d.', pid);
            end
        else
            error('Could not determine unique training rotation mapping for subject %d.', pid);
        end

        mapping_summary = [mapping_summary; ...
            table(pid,rot500,rot1500,flipFlag,mappingLabel, ...
            'VariableNames', ...
            {'id','rot500','rot1500','willFlip','mappingLabel'})]; %#ok<AGROW>
    end

    fprintf('Auto-flip IDs: ');
    disp(auto_flip_ids);

    isFlip = ismember(T.id, auto_flip_ids);

    for k = 1:numel(settings.flipVars)
        vn = settings.flipVars{k};
        T.(vn)(isFlip) = -T.(vn)(isFlip);
    end

    % Build model data structure.
    maxTrials = 0;
    for s = 1:nSubs
        maxTrials = max(maxTrials, sum(T.id==participants(s)));
    end

    data.ha       = nan(nSubs, maxTrials);
    data.timing   = nan(nSubs, maxTrials);
    data.rotation = nan(nSubs, maxTrials);
    data.phase    = nan(nSubs, maxTrials);
    data.subTrial = nan(nSubs, maxTrials);
    data.id       = double(participants(:));

    hasRT = ismember('rtGoSig', T.Properties.VariableNames);
    if hasRT
        data.rt = nan(nSubs, maxTrials);
    end

    for s = 1:nSubs
        pid = participants(s);
        S = sortrows(T(T.id==pid,:), 'subTrial');
        nT = height(S);
        idx = 1:nT;

        data.ha(s,idx) = S.handAng';
        data.timing(s,idx) = S.goDelay';
        data.rotation(s,idx) = S.rotDir';
        data.phase(s,idx) = S.phase';
        data.subTrial(s,idx) = S.subTrial';

        if hasRT
            data.rt(s,idx) = S.rtGoSig';
        end
    end

    prepInfo.auto_flip_ids = auto_flip_ids(:);
    prepInfo.mapping_summary = mapping_summary;
end

%% FIT INITIALIZERS

function f = initTimingFits(n, settings)
    f.weber = nan(n,1);
    f.rtmin = nan(n,1);
    f.rtscale = nan(n,1);
    f.error = nan(n,1);
    f.sim_rt = nan(n, numel(settings.intervals_sec));
end

function f = initModel1Fits(n)
    f.A = nan(n,1);
    f.B = nan(n,1);
    f.Ac = nan(n,1);
    f.Bc = nan(n,1);
    f.error = nan(n,1);
end

function f = initModel2Fits(n)
    f.A = nan(n,1);
    f.B = nan(n,1);
    f.Aw = nan(n,1);
    f.Bw = nan(n,1);
    f.error = nan(n,1);
end

function p = getModel2Params(fits, si)
    p = [fits.A(si), fits.B(si), fits.Aw(si), fits.Bw(si)];
end

%% PER-SUBJECT FITTING

function fits = fitTimingModel(fits, si, timing, rt, settings)

    fprintf('  Timing model...\n');

    best_err = inf;
    best_soln = [];

    for iter = 1:settings.iterations_timing
        p0 = [unifrnd(0.05,0.3), ...
              unifrnd(0.1,0.5), ...
              unifrnd(0,1)];

        LB = [0.05, 0.1, 0];
        UB = [0.3, 0.5, 1];

        opts = optimoptions('fmincon','Display','off');

        [soln, err] = fmincon( ...
            @fitTimingObj, p0, [],[],[],[], LB, UB, [], opts, ...
            timing, rt, settings.intervals_sec);

        if err < best_err
            best_err = err;
            best_soln = soln;
        end
    end

    fits.weber(si) = best_soln(1);
    fits.rtmin(si) = best_soln(2);
    fits.rtscale(si) = best_soln(3);
    fits.error(si) = best_err;

    [~, fits.sim_rt(si,:)] = ...
        fitTimingObj(best_soln, timing, rt, settings.intervals_sec);
end

function fits = fitModel1(fits, si, ha, timing, rotation, settings)

    fprintf('  Subjective Timing Model (state-space + discrete contexts)...\n');

    best_err = inf;
    best_soln = [];
    lp = settings.learning_phase;

    LB = [0 0 0 0];
    UB = [1 0.5 1 1];

    opts = optimoptions( ...
        'fmincon','Display','off','Algorithm','sqp', ...
        'MaxFunctionEvaluations',5000);

    for iter = 1:settings.iterations_model1
        p0 = [unifrnd(LB(1),UB(1)), ...
              unifrnd(LB(2),UB(2)), ...
              unifrnd(LB(3),UB(3)), ...
              unifrnd(LB(4),UB(4))];

        [soln, err] = fmincon( ...
            @fitModel1Obj, p0, [],[],[],[], LB, UB, [], opts, ...
            ha, timing, rotation, lp, settings);

        if err < best_err
            best_err = err;
            best_soln = soln;
        end
    end

    fits.A(si) = best_soln(1);
    fits.B(si) = best_soln(2);
    fits.Ac(si) = best_soln(3);
    fits.Bc(si) = best_soln(4);
    fits.error(si) = best_err;
end

function fits = fitModel2(fits, si, ha, timing, rotation, basisCache, settings)

    fprintf('  Granule-Cell Basis Model (state-space + granule basis)...\n');

    best_err = inf;
    best_soln = [];
    lp = settings.learning_phase;

    LB = [0 0 0 0];
    UB = [1 0.5 1 0.10];

    opts = optimoptions( ...
        'fmincon','Display','off','Algorithm','sqp', ...
        'MaxFunctionEvaluations',5000);

    for iter = 1:settings.iterations_model2
        p0 = [unifrnd(LB(1),UB(1)), ...
              unifrnd(LB(2),UB(2)), ...
              unifrnd(LB(3),UB(3)), ...
              unifrnd(LB(4),UB(4))];

        [soln, err] = fmincon( ...
            @fitModel2Obj, p0, [],[],[],[], LB, UB, [], opts, ...
            ha, timing, rotation, lp, basisCache, settings);

        if err < best_err
            best_err = err;
            best_soln = soln;
        end
    end

    fits.A(si) = best_soln(1);
    fits.B(si) = best_soln(2);
    fits.Aw(si) = best_soln(3);
    fits.Bw(si) = best_soln(4);
    fits.error(si) = best_err;
end

%% OBJECTIVE FUNCTIONS

function [err, pred_RT] = fitTimingObj(params, timing, rt, intervals)

    t = linspace(0, 2, 2000);
    w = params(1);

    pdf = zeros(size(t));

    for i = 1:numel(intervals)
        sigma = max(w*intervals(i), 1e-6);
        pdf = pdf + exp(-((t-intervals(i)).^2)/(2*sigma^2));
    end

    pdf = pdf / trapz(t, pdf);
    cdf = cumtrapz(t, pdf);
    haz = pdf ./ max(1-cdf, eps);
    haz(~isfinite(haz)) = 0;

    a = params(2);
    b = params(3);

    pred_RT = nan(1,numel(intervals));

    for k = 1:numel(intervals)
        idx = max(1, min(numel(haz), round(1000*intervals(k))));
        pred_RT(k) = a + b/(haz(idx)+1e-6);
    end

    si1 = find(timing == intervals(1)*1000);
    si2 = find(timing == intervals(2)*1000);

    si1(si1<21|si1>420) = [];
    si2(si2<21|si2>420) = [];

    real_RT = [nanmedian(rt(si1)), nanmedian(rt(si2))];

    err = sqrt(mean((pred_RT - real_RT).^2, 'omitnan'));
end

function err = fitModel1Obj(params, ha, timing, rotation, lp, settings)
    sim_ha = runModel1(params, timing, rotation, lp, settings);
    err = sqrt(mean((sim_ha - ha(lp)).^2, 'omitnan'));
end

function sim_ha = runModel1(params, timing, rotation, lp, settings)

    A = params(1);
    B = params(2);
    Ac = params(3);
    Bc = params(4);

    N = numel(lp);

    Vp = zeros(1,N+1);
    V5 = zeros(1,N+1);
    V15 = zeros(1,N+1);

    sim_ha = nan(1,N);

    for n = 1:N

        ti = lp(n);
        d = timing(ti);
        rot = rotation(ti);

        if d == 500
            y = Vp(n) + V5(n);
        elseif d == 1500
            y = Vp(n) + V15(n);
        else
            y = Vp(n);
        end

        sim_ha(n) = y;

        R = -settings.rotAmp * rot;
        e = R - y;

        Vp(n+1) = A*Vp(n) + B*e;

        V5(n+1) = Ac*V5(n);
        V15(n+1) = Ac*V15(n);

        if d == 500
            V5(n+1) = V5(n+1) + Bc*e;
        elseif d == 1500
            V15(n+1) = V15(n+1) + Bc*e;
        end
    end
end

function err = fitModel2Obj(params, ha, timing, rotation, lp, ...
    basisCache, settings)

    sim_ha = runModel2Learning( ...
        params, timing, rotation, lp, basisCache, settings);

    err = sqrt(mean((sim_ha - ha(lp)).^2, 'omitnan'));
end

function [sim_ha, wFinal, xFinal] = ...
    runModel2Learning(params, timing, rotation, lp, basisCache, settings)

    A = params(1);
    B = params(2);
    Aw = params(3);
    Bw = params(4);

    N = numel(lp);
    nC = numel(basisCache.centers);

    x = 0;
    w = zeros(1,nC);

    sim_ha = nan(1,N);

    for n = 1:N

        ti = lp(n);
        d = timing(ti);
        rot = rotation(ti);

        f = getCachedBasis(d, basisCache);

        y = x + f*w';
        sim_ha(n) = y;

        R = -settings.rotAmp*rot;
        e = R - y;

        x = A*x + B*e;
        w = Aw.*w + Bw.*e.*f;
    end

    wFinal = w;
    xFinal = x;
end

%% FULL SIMULATION

function [SIM, exp_rot] = ...
    simulateModel1Full(params, timing, rotation, phase, ...
    weber, t_ms, settings)

    A = params(1);
    B = params(2);
    Ac = params(3);
    Bc = params(4);

    nT = numel(timing);
    t = t_ms./1000;

    safe_sigma = @(w,c) max(w*c,1e-6);

    k1 = exp(-((t-0.5).^2)/(2*safe_sigma(weber,0.5)^2));
    k2 = exp(-((t-1.5).^2)/(2*safe_sigma(weber,1.5)^2));

    Vp = zeros(1,nT+1);
    V5 = zeros(1,nT+1);
    V15 = zeros(1,nT+1);

    SIM = nan(1,nT);
    exp_rot = nan(nT,numel(t_ms));

    for n = 1:nT

        d = timing(n);
        rot = rotation(n);

        if d==500
            y = Vp(n)+V5(n);
        elseif d==1500
            y = Vp(n)+V15(n);
        else
            y = Vp(n);
        end

        SIM(n) = y;

        exp_rot(n,:) = V5(n)*k1 + V15(n)*k2;

        doUp = ~(settings.freezeDuringPhase3 && phase(n)==3);

        if doUp

            R = -settings.rotAmp * rot;
            e = R - y;

            Vp(n+1) = A*Vp(n)+B*e;

            V5(n+1) = Ac*V5(n);
            V15(n+1) = Ac*V15(n);

            if d==500
                V5(n+1) = V5(n+1)+Bc*e;
            elseif d==1500
                V15(n+1) = V15(n+1)+Bc*e;
            end

        else

            Vp(n+1) = Vp(n);
            V5(n+1) = V5(n);
            V15(n+1) = V15(n);
        end
    end
end

function [data_gen, model_gen] = ...
    computePhase3FromExpected(data, exp_all, t_ms, ...
    testDelays, gain, mode)

    [nS,nT] = size(data.ha);

    data_gen = nan(nS,numel(testDelays));
    model_gen = nan(nS,numel(testDelays));

    for si = 1:nS

        ha = data.ha(si,:);
        timing = data.timing(si,:);
        ixP = data.phase(si,:)==3;

        if mode=="mean_center"
            mu_d = mean(ha(ixP),'omitnan');
        else
            mu_d = 0;
        end

        for k = 1:numel(testDelays)
            data_gen(si,k) = ...
                mean(ha(ixP&timing==testDelays(k))-mu_d,'omitnan');
        end

        mv = nan(1,nT);

        for n = find(ixP)
            cv = squeeze(exp_all(si,n,:)).';
            mv(n) = gain*interp1( ...
                t_ms, cv, timing(n), 'linear', 'extrap');
        end

        if mode=="mean_center"
            mu_m = mean(mv(ixP),'omitnan');
        else
            mu_m = 0;
        end

        for k = 1:numel(testDelays)
            model_gen(si,k) = ...
                mean(mv(ixP&timing==testDelays(k))-mu_m,'omitnan');
        end
    end
end

function model_gen = ...
    computePhase3GranuleDirect(data, fits_granulecell, ...
    basisCache, settings, gain, mode)

    if nargin < 5 || isempty(gain)
        gain = 1;
    end

    if nargin < 6 || isempty(mode)
        mode = settings.phase3Mode;
    end

    nS = size(data.ha,1);
    model_gen = nan(nS,numel(settings.testDelays));

    for si = 1:nS

        timing = data.timing(si,:);
        rotation = data.rotation(si,:);
        phase = data.phase(si,:);

        params = getModel2Params(fits_granulecell, si);

        [~,wFinal,xFinal] = runModel2Learning( ...
            params, timing, rotation, ...
            settings.learning_phase, basisCache, settings);

        ixP = phase==3;
        mv = nan(size(timing));

        for n = find(ixP)

            d = timing(n);

            if ismember(d,basisCache.delays)
                f = getCachedBasis(d,basisCache);
            else
                f = granuleRBF( ...
                    d, basisCache.centers, ...
                    settings.granule.baseSigma, ...
                    settings.granule.widthSlope, ...
                    settings.granule.decayTau, ...
                    settings.granule.timingNoise);
            end

            mv(n) = gain*(xFinal + f*wFinal');
        end

        if mode=="mean_center"
            mu_m = mean(mv(ixP),'omitnan');
        else
            mu_m = 0;
        end

        for k = 1:numel(settings.testDelays)
            model_gen(si,k) = mean( ...
                mv(ixP&timing==settings.testDelays(k))-mu_m, ...
                'omitnan');
        end
    end
end

function stats = summarizeGenFit(data_gen, model_gen)

    md = mean(data_gen,1,'omitnan');
    mm = mean(model_gen,1,'omitnan');

    ok = isfinite(md) & isfinite(mm);

    stats.R2_group = 1 - ...
        sum((md(ok)-mm(ok)).^2) / ...
        sum((md(ok)-mean(md(ok))).^2);
end

%% GRANULE-CELL BASIS

function basisCache = makeGranuleBasisCache(settings)

    allDelays = unique( ...
        [settings.trainDelays settings.testDelays]);

    centers = settings.granule.centers_ms;

    sigma = settings.granule.baseSigma + ...
        settings.granule.widthSlope.*(centers - min(centers));

    sigma_eff = sqrt( ...
        sigma.^2 + settings.granule.timingNoise.^2);

    basisCache.delays = allDelays;
    basisCache.centers = centers;
    basisCache.Phi = ...
        nan(numel(allDelays), numel(centers));

    for i = 1:numel(allDelays)

        d = allDelays(i);

        basisCache.Phi(i,:) = ...
            exp(-centers./settings.granule.decayTau) .* ...
            exp(-0.5.*((d-centers)./sigma_eff).^2);
    end
end

function phi = getCachedBasis(d, basisCache)

    [tf,idx] = ismember(d,basisCache.delays);

    if tf
        phi = basisCache.Phi(idx,:);
    else
        error('Delay %g not in cache',d);
    end
end

function f = granuleRBF(t, centers, ...
    baseSigma, widthSlope, decayTau, timingNoise)

    sigma = baseSigma + ...
        widthSlope.*(centers - min(centers));

    sigma_eff = sqrt(sigma.^2 + timingNoise.^2);

    f = exp(-centers./decayTau) .* ...
        exp(-0.5.*((t-centers)./sigma_eff).^2);
end

function vals = cleanUnique(x)
    x = x(~isnan(x));
    vals = unique(x);
    vals = vals(:)';
end
