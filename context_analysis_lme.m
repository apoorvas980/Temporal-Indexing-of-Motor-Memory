
% mixed-effects models and participant-level
% regressions from the analysis-ready CSV files 

clear; clc;

BASE_DIR = fileparts(mfilename('fullpath'));

explicit = readtable(fullfile(BASE_DIR, 'data', 'behavioral_explicit_lme.csv'));
implicit = readtable(fullfile(BASE_DIR, 'data', 'behavioral_implicit_lme.csv'));

explicit.Subject = categorical(explicit.Subject);
implicit.Subject = categorical(implicit.Subject);

% 1. MIXED-EFFECTS MODELS 

formula = ['HA ~ Time*Rot_t*Rot_t1 + Time*Rot_t2 + ' ...
           '(Time*Rot_t*Rot_t1 + Time*Rot_t2 | Subject)'];

fprintf('EXTERNALLY-CUED / REDUCED-UNCERTAINTY LME\n');
lme_explicit = fitlme(explicit, formula, 'FitMethod', 'REML');
disp(lme_explicit);

fprintf('INTERNALLY-MONITORED / HIGH UNCERTAINTY LME\n');
lme_implicit = fitlme(implicit, formula, 'FitMethod', 'REML');
disp(lme_implicit);

fprintf('COEFFICIENTS\n');

print_key_lme_stats(lme_explicit, 'EXTERNALLY-CUED / reduced uncertainty');
print_key_lme_stats(lme_implicit, 'INTERNALLY-MONITORED / hihg uncertainty');

% 2. PARTICIPANT-LEVEL REGRESSIONS
%
% Per participant:
%   HA ~ context + time + n-1 + n-2
%        + context*time + n-1*time + n-2*time
%
% Predictors are z-scored within participant

betas_explicit = participant_regressions(explicit);
betas_implicit = participant_regressions(implicit);

fprintf('PARTICIPANT-LEVEL REGRESSIONS: EXTERNALLY-CUED\n');
disp(betas_explicit);
print_beta_summary(betas_explicit);

fprintf('PARTICIPANT-LEVEL REGRESSIONS: INTERALLY-MONITORED\n');
disp(betas_implicit);
print_beta_summary(betas_implicit);


% Local functions

function print_key_lme_stats(lme, label)

    C = lme.Coefficients;

    names = {'Rot_t','Rot_t1','Time:Rot_t'};
    labels = {'Context main effect', ...
              'Adaptation main effect (n-1)', ...
              'Context x Time interaction'};

    fprintf('\n%s\n', label);

    for k = 1:numel(names)
        idx = find(strcmp(C.Name, names{k}), 1);

        if isempty(idx) && strcmp(names{k}, 'Time:Rot_t')
            idx = find(strcmp(C.Name, 'Rot_t:Time'), 1);
        end

        if isempty(idx)
            warning('Could not find coefficient: %s', names{k});
            continue
        end

        fprintf('%s:\n', labels{k});
        fprintf('  beta = %.6f\n', C.Estimate(idx));
        fprintf('  SE   = %.6f\n', C.SE(idx));
        fprintf('  t(%g) = %.3f\n', C.DF(idx), C.tStat(idx));
        fprintf('  p = %.6g\n', C.pValue(idx));
        fprintf('  95%% CI = [%.6f, %.6f]\n', C.Lower(idx), C.Upper(idx));
    end
end

function B = participant_regressions(T)

    subjects = categories(T.Subject);

    Subject = NaN(numel(subjects),1);
    context = NaN(numel(subjects),1);
    contextTime = NaN(numel(subjects),1);
    adaptation = NaN(numel(subjects),1);
    adaptationTime = NaN(numel(subjects),1);
    twoBackAdaptation = NaN(numel(subjects),1);
    twoBackAdaptationTime = NaN(numel(subjects),1);

    for i = 1:numel(subjects)

        idx = T.Subject == subjects{i};
        P = T(idx,:);

        y  = P.HA;
        c  = zscore_safe(P.Rot_t);
        tm = zscore_safe(P.Time);
        a1 = zscore_safe(P.Rot_t1);
        a2 = zscore_safe(P.Rot_t2);

        if any(isnan(c)) || any(isnan(tm)) || ...
           any(isnan(a1)) || any(isnan(a2))
            continue
        end

        X = [ones(height(P),1), ...
             c, tm, a1, a2, ...
             c.*tm, a1.*tm, a2.*tm];

        b = X \ y;

        Subject(i) = str2double(subjects{i});
        context(i) = b(2);
        contextTime(i) = b(6);
        adaptation(i) = b(4);
        adaptationTime(i) = b(7);
        twoBackAdaptation(i) = b(5);
        twoBackAdaptationTime(i) = b(8);
    end

    B = table(Subject, context, contextTime, adaptation, adaptationTime, ...
        twoBackAdaptation, twoBackAdaptationTime);
end

function z = zscore_safe(x)
    s = std(x, 0, 1, 'omitnan');
    if s == 0 || isnan(s)
        z = NaN(size(x));
    else
        z = (x - mean(x, 'omitnan')) ./ s;
    end
end

function print_beta_summary(B)

    vars = B.Properties.VariableNames(2:end);

    for i = 1:numel(vars)

        x = B.(vars{i});
        x = x(~isnan(x));

        [~, p_t, ci, stats] = ttest(x);
        p_w = signrank(x);

        fprintf('\n%s\n', vars{i});
        fprintf('  mean = %.6f\n', mean(x));
        fprintf('  SE   = %.6f\n', std(x) / sqrt(numel(x)));
        fprintf('  t(%d) = %.3f, p = %.6g\n', stats.df, stats.tstat, p_t);
        fprintf('  95%% CI = [%.6f, %.6f]\n', ci(1), ci(2));
        fprintf('  Wilcoxon signrank p = %.6g\n', p_w);
    end
end
