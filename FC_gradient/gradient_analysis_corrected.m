
clear; clc; close all;

%% ------------------------------------------------------------------------
% Paths and analysis parameters
% -------------------------------------------------------------------------
addpath('Demography_data');
addpath('data_fmri');
addpath(genpath('BrainSpace'));

if ~exist('results', 'dir')
    mkdir('results');
end

n_components = 10;
sparsity_main = 90;       % BrainSpace sparsity=90 -> retain largest ~10% per row
alpha_fdr = 0.05;
run_sensitivity = true;   % set false if you only want the main analysis
sensitivity_sparsities = [95 90 85 80]; % retain top 5%, 10%, 15%, 20%
random_state = 20260918;

%% ------------------------------------------------------------------------
% Read subject information and FC index table
% -------------------------------------------------------------------------
info_table = readtable('selected_info_rand.xlsx', 'Sheet', 'Sheet1');
fmri_selected = readtable('AllSchaefSubjectsCombined.xlsx', 'Sheet', 'Sheet1');

info_id = string(info_table.Subid);
fmri_id = string(fmri_selected.Subid);

% Preserve the order of info_table
[common_id, idx_info, idx_fmri] = intersect(info_id, fmri_id, 'stable');

% Predefined exclusions
exclude_subids = ["HC_005", "HC_01", "HC_06", "HC_09", "HC_15", "HC_23", ...
                  "PoC_013", "PoC_018", "PoC_020", "PoN_030"];

keep_common = ~ismember(common_id, exclude_subids);
common_id = common_id(keep_common);
idx_info = idx_info(keep_common);
idx_fmri = idx_fmri(keep_common);

% Demographic/clinical table in a single fixed order
info_common = info_table(idx_info, :);

%% ------------------------------------------------------------------------
% Load and align FC matrices to the same subject order
% -------------------------------------------------------------------------
fmri_struct = load('AllSchaefDataCombined.mat');
if ~isfield(fmri_struct, 'allDataCombined')
    error('AllSchaefDataCombined.mat does not contain variable allDataCombined.');
end
fmri_matrix = fmri_struct.allDataCombined;

% Critical check: third dimension must correspond to fmri_selected rows
if size(fmri_matrix, 3) ~= height(fmri_selected)
    error(['FC subject count (%d) does not match AllSchaefSubjectsCombined.xlsx rows (%d). ' ...
           'The matrix-to-Subid correspondence must be verified before analysis.'], ...
           size(fmri_matrix, 3), height(fmri_selected));
end

% idx_fmri already corresponds exactly to info_common because intersect(...,'stable')
fmri_matrix_common = fmri_matrix(:, :, idx_fmri);

nroi = size(fmri_matrix_common, 1);
if size(fmri_matrix_common, 2) ~= nroi
    error('FC matrices are not square.');
end
if nroi ~= 132
    warning('Expected 132 regions, but found %d regions.', nroi);
end

nsubs_initial = size(fmri_matrix_common, 3);
if height(info_common) ~= nsubs_initial
    error('Mismatch between demographic rows and FC matrices after matching.');
end

fprintf('Matched subjects before QC: %d\n', nsubs_initial);
fprintf('FC matrix size: %d x %d x %d\n', size(fmri_matrix_common));

%% ------------------------------------------------------------------------
% Extract covariates and run subject-level QC
% -------------------------------------------------------------------------
age = double(info_common.Age);
education = double(info_common.Education);
group = string(info_common.group);
subid = string(info_common.Subid);

valid_subjects = true(nsubs_initial, 1);
qc_reason = strings(nsubs_initial, 1);
fmri_cell = cell(nsubs_initial, 1);

for subj = 1:nsubs_initial
    A = double(fmri_matrix_common(:, :, subj));

    if ~isequal(size(A), [nroi nroi])
        valid_subjects(subj) = false;
        qc_reason(subj) = "matrix_size_mismatch";
        continue;
    end

    % Exclude non-finite matrices. Do NOT replace missing values with zero.
    if any(~isfinite(A(:)))
        valid_subjects(subj) = false;
        qc_reason(subj) = "NaN_or_Inf";
        continue;
    end

    % Enforce numerical symmetry before gradient construction
    A = (A + A') / 2;
    A(1:nroi+1:end) = 0;

    % Exclude matrices containing an all-zero row
    if any(all(abs(A) < eps, 2))
        valid_subjects(subj) = false;
        qc_reason(subj) = "all_zero_row";
        continue;
    end

    fmri_cell{subj} = A;
end

% Keep only intended groups
valid_group = ismember(group, ["HCs", "PMI"]);
qc_reason(~valid_group & valid_subjects) = "unexpected_group";
valid_subjects = valid_subjects & valid_group;

% Also require finite covariates
valid_cov = isfinite(age) & isfinite(education);
qc_reason(~valid_cov & valid_subjects) = "missing_covariate";
valid_subjects = valid_subjects & valid_cov;

% Save QC table
qc_table = table(subid, group, age, education, valid_subjects, qc_reason, ...
    'VariableNames', {'Subid','Group','Age','Education','Included','QCReason'});
writetable(qc_table, fullfile('results','Subject_QC.csv'));

% Apply the SAME subject mask to every variable
fmri_cell_valid = fmri_cell(valid_subjects);
subid_valid = subid(valid_subjects);
group_valid = group(valid_subjects);
age_valid = age(valid_subjects);
education_valid = education(valid_subjects);

nsubs = numel(fmri_cell_valid);
fprintf('Subjects retained after QC: %d\n', nsubs);
fprintf('  HCs: %d\n', sum(group_valid == "HCs"));
fprintf('  PMI: %d\n', sum(group_valid == "PMI"));

if sum(group_valid == "HCs") < 2 || sum(group_valid == "PMI") < 2
    error('Too few subjects in one or both groups after QC.');
end

% Explicit subject-order audit file
analysis_subjects = table((1:nsubs)', subid_valid, group_valid, age_valid, education_valid, ...
    'VariableNames', {'AnalysisRow','Subid','Group','Age','Education'});
writetable(analysis_subjects, fullfile('results','Analysis_Subject_Order.csv'));

%% ------------------------------------------------------------------------
% Main gradient analysis
% -------------------------------------------------------------------------
[gradient_all, explained_variance, gm_group, gm_all] = ...
    run_gradient_pipeline(fmri_cell_valid, sparsity_main, n_components, random_state);

% gradient_all is nSubjects x nRegions, in EXACTLY the same order as group_valid
assert(size(gradient_all,1) == nsubs, 'Gradient subject count mismatch.');
assert(size(gradient_all,2) == nroi, 'Gradient region count mismatch.');

% Save main gradients and subject order
save(fullfile('results','Gradient_Main.mat'), ...
    'gradient_all','explained_variance','subid_valid','group_valid', ...
    'age_valid','education_valid','sparsity_main','gm_group','gm_all','-v7.3');

%% ------------------------------------------------------------------------
% Explained variance plot
% -------------------------------------------------------------------------
var_HC = explained_variance(group_valid == "HCs", :);
var_PMI = explained_variance(group_valid == "PMI", :);

mean_var_HC = mean(var_HC, 1, 'omitnan');
mean_var_PMI = mean(var_PMI, 1, 'omitnan');
se_var_HC = std(var_HC, 0, 1, 'omitnan') ./ sqrt(size(var_HC,1));
se_var_PMI = std(var_PMI, 0, 1, 'omitnan') ./ sqrt(size(var_PMI,1));

figure('Color','w');
hold on;
b = bar(1:n_components, [mean_var_HC; mean_var_PMI]', 'grouped');
xHC = b(1).XEndPoints;
xPMI = b(2).XEndPoints;
errorbar(xHC, mean_var_HC, se_var_HC, 'k.', 'LineWidth', 1.2);
errorbar(xPMI, mean_var_PMI, se_var_PMI, 'k.', 'LineWidth', 1.2);
plot(xHC, mean_var_HC, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
plot(xPMI, mean_var_PMI, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
legend('HCs','PMI','Location','northeast');
xlabel('Gradient component');
ylabel('Explained variance proportion');
title(sprintf('Explained variance by group (top %d%% retained)', 100-sparsity_main));
set(gca,'XTick',1:n_components);
grid on; box on;
hold off;
exportgraphics(gcf, fullfile('results','ExplainedVariance_byGroup.pdf'), 'ContentType','vector');

%% ------------------------------------------------------------------------
% Covariate-adjusted residuals for visualization only
% Residualize age + education, then z-score EACH REGION across subjects.
% -------------------------------------------------------------------------
X = [ones(nsubs,1), age_valid, education_valid];
resid_region = nan(nsubs, nroi);

for region = 1:nroi
    y = gradient_all(:, region);
    b = X \ y;
    resid_region(:, region) = y - X*b;
end

% Correct direction: standardize across subjects for each region
resid_region_z = zscore(resid_region, 0, 1);

% Distribution plot
x_HC = reshape(resid_region_z(group_valid == "HCs", :), [], 1);
x_PMI = reshape(resid_region_z(group_valid == "PMI", :), [], 1);

figure('Color','w');
h1 = histogram(x_HC, 'Normalization','probability','BinWidth',0.15);
hold on;
h2 = histogram(x_PMI, 'Normalization','probability','BinWidth',0.15);
legend('HCs','PMI','Location','best');
xlabel('FC gradient residuals (z)');
ylabel('Relative frequency');
set(gca,'FontName','Arial');
box on;
hold off;
exportgraphics(gcf, fullfile('results','GradientResidual_Distribution.pdf'), 'ContentType','vector');

%% ------------------------------------------------------------------------
% Optional global mean G1 analysis
% Note: a global mean gradient is not the primary regional inference.
% -------------------------------------------------------------------------
mean_gradient = mean(gradient_all, 2, 'omitnan');
Group = categorical(group_valid);
Group = reordercats(Group, {'HCs','PMI'});

tbl_global = table(mean_gradient, age_valid, education_valid, Group, ...
    'VariableNames', {'MeanGradient','Age','Education','Group'});
lm_global = fitlm(tbl_global, 'MeanGradient ~ Age + Education + Group');

idx_group_global = find(contains(string(lm_global.CoefficientNames), 'Group_PMI'), 1);
if isempty(idx_group_global)
    error('Could not identify the PMI-vs-HCs group coefficient in the global model.');
end

beta_global = lm_global.Coefficients.Estimate(idx_group_global);
t_global = lm_global.Coefficients.tStat(idx_group_global);
p_global = lm_global.Coefficients.pValue(idx_group_global);

fprintf('\nGlobal mean G1: PMI vs HCs adjusted for age + education\n');
fprintf('  beta = %.6f, t = %.4f, p = %.6g\n', beta_global, t_global, p_global);

% Mean of standardized regional residuals per subject
resid_mean = mean(resid_region_z, 2, 'omitnan');
figure('Color','w');
boxplot(resid_mean, cellstr(group_valid), 'Labels', {'HCs','PMI'});
ylabel('Mean regional residual z-score');
title('Mean covariate-adjusted regional gradient residual');
box on;
exportgraphics(gcf, fullfile('results','MeanResidual_byGroup.pdf'), 'ContentType','vector');

%% ------------------------------------------------------------------------
% Regional PMI-vs-HCs group differences
% Model: G1 ~ Age + Education + Group
% Reference group = HCs, therefore positive t/beta means PMI > HCs.
% -------------------------------------------------------------------------
[beta_group, t_group, p_group, q_group] = ...
    regional_group_stats(gradient_all, age_valid, education_valid, group_valid);

sigregs = find(q_group < alpha_fdr);
fprintf('\nRegional analysis: %d / %d regions significant at BH-FDR q < %.3f\n', ...
    numel(sigregs), nroi, alpha_fdr);

% Load region labels if available
if exist('regionLabels.mat','file')
    S = load('regionLabels.mat');
    if isfield(S,'regionLabels') && numel(S.regionLabels) == nroi
        region_labels = string(S.regionLabels(:));
    else
        region_labels = "Region_" + string((1:nroi)');
        warning('regionLabels.mat found, but regionLabels is missing or has the wrong length.');
    end
else
    region_labels = "Region_" + string((1:nroi)');
    warning('regionLabels.mat not found; generic region names will be used.');
end

regional_results = table((1:nroi)', region_labels, beta_group, t_group, p_group, q_group, ...
    q_group < alpha_fdr, ...
    'VariableNames', {'RegionIndex','RegionName','Beta_PMI_minus_HC','T_PMI_minus_HC', ...
                     'P_uncorrected','Q_BHFDR','Significant_FDR05'});
writetable(regional_results, fullfile('results','Regional_Group_Differences.xlsx'));

% Save significant gradient values for downstream correlation analyses
ms_gradient_sigregs = gradient_all(:, sigregs);
region_names_sig = region_labels(sigregs);
save(fullfile('results','ms_gradient_sigregs.mat'),'ms_gradient_sigregs','sigregs');
save(fullfile('results','region_names_sig.mat'),'region_names_sig');

% Map for visualization: non-significant regions set to zero
sigtstat = t_group;
sigtstat(q_group >= alpha_fdr) = 0;

%% ------------------------------------------------------------------------
% Mean regional G1 within each group
% Correct labels: controls = HCs; cases = PMI
% -------------------------------------------------------------------------
meanGradient_HC = mean(gradient_all(group_valid == "HCs", :), 1, 'omitnan')';
meanGradient_PMI = mean(gradient_all(group_valid == "PMI", :), 1, 'omitnan')';

%% ------------------------------------------------------------------------
% Quadrant summary (descriptive only)
% x = mean G1 in controls; y = PMI-vs-HCs t-statistic
% -------------------------------------------------------------------------
xvalues = meanGradient_HC;
yvalues = t_group;

a = mean((xvalues < 0) & (yvalues > 0));
b = mean((xvalues > 0) & (yvalues > 0));
c = mean((xvalues < 0) & (yvalues < 0));
d = mean((xvalues > 0) & (yvalues < 0));

fprintf('\nQuadrant proportions:\n');
fprintf('  x<0, t>0: %.3f\n', a);
fprintf('  x>0, t>0: %.3f\n', b);
fprintf('  x<0, t<0: %.3f\n', c);
fprintf('  x>0, t<0: %.3f\n', d);

%% ------------------------------------------------------------------------
% Spatial association between control mean G1 and regional group t-values
% IMPORTANT: p_parametric below is NOT a spin-test p-value.
% -------------------------------------------------------------------------
[r_spatial, p_parametric] = corr(meanGradient_HC, t_group, ...
    'Rows','complete', 'Type','Pearson');

figure('Color','w');
scatter(meanGradient_HC, t_group, 36, 'x');
hold on;
p_fit = polyfit(meanGradient_HC, t_group, 1);
y_fit = polyval(p_fit, meanGradient_HC);
plot(meanGradient_HC, y_fit, '-', 'LineWidth', 1.5);
yline(0, '--');
xline(0, '--');
xlabel('Mean G1 in HCs');
ylabel('PMI-HCs t-value');
title(sprintf('Spatial association: r = %.3f, parametric p = %.3g', ...
    r_spatial, p_parametric));
grid on; box on;
hold off;
exportgraphics(gcf, fullfile('results','SpatialAssociation_HCGradient_vs_Tstat.pdf'), ...
    'ContentType','vector');

fprintf('\nSpatial correlation across regions:\n');
fprintf('  Pearson r = %.4f\n', r_spatial);
fprintf('  Parametric p = %.6g (NOT P_spin)\n', p_parametric);

%% ------------------------------------------------------------------------
% Save visualization matrix
% -------------------------------------------------------------------------
data2vis = [meanGradient_HC, meanGradient_PMI, t_group, sigtstat];
save(fullfile('results','data2vis.mat'), 'data2vis');

%% ------------------------------------------------------------------------
% Sensitivity analysis for proportional sparsity threshold
% sparsity 95/90/85/80 = retain top 5/10/15/20%
% -------------------------------------------------------------------------
if run_sensitivity
    ns = numel(sensitivity_sparsities);
    retained_percent = 100 - sensitivity_sparsities(:);
    n_sig_fdr = nan(ns,1);
    t_corr_with_main = nan(ns,1);
    t_corr_p = nan(ns,1);
    sig_jaccard_with_main = nan(ns,1);

    sensitivity_t = nan(nroi, ns);
    sensitivity_q = nan(nroi, ns);

    main_sig = q_group < alpha_fdr;

    fprintf('\nRunning sparsity sensitivity analyses...\n');
    for k = 1:ns
        sp = sensitivity_sparsities(k);
        fprintf('  sparsity=%d (retain top %d%%)\n', sp, 100-sp);

        if sp == sparsity_main
            grad_s = gradient_all;
            t_s = t_group;
            q_s = q_group;
        else
            [grad_s, ~] = run_gradient_pipeline( ...
                fmri_cell_valid, sp, n_components, random_state);
            [~, t_s, ~, q_s] = regional_group_stats( ...
                grad_s, age_valid, education_valid, group_valid);
        end

        sensitivity_t(:,k) = t_s;
        sensitivity_q(:,k) = q_s;
        n_sig_fdr(k) = sum(q_s < alpha_fdr);

        [t_corr_with_main(k), t_corr_p(k)] = corr(t_group, t_s, ...
            'Rows','complete', 'Type','Pearson');

        sig_s = q_s < alpha_fdr;
        union_n = sum(main_sig | sig_s);
        if union_n == 0
            sig_jaccard_with_main(k) = NaN;
        else
            sig_jaccard_with_main(k) = sum(main_sig & sig_s) / union_n;
        end
    end

    sensitivity_summary = table(sensitivity_sparsities(:), retained_percent, ...
        n_sig_fdr, t_corr_with_main, t_corr_p, sig_jaccard_with_main, ...
        'VariableNames', {'BrainSpaceSparsity','TopPercentRetained','N_FDR_Significant', ...
                         'TmapCorrelationWithMain','CorrelationP', ...
                         'FDR_Significant_JaccardWithMain'});

    writetable(sensitivity_summary, fullfile('results','Sparsity_Sensitivity_Summary.xlsx'));
    save(fullfile('results','Sparsity_Sensitivity.mat'), ...
        'sensitivity_sparsities','sensitivity_t','sensitivity_q','sensitivity_summary');

    disp(sensitivity_summary);
end

%% ------------------------------------------------------------------------
% Final note on a TRUE spin test
% -------------------------------------------------------------------------
fprintf('\nIMPORTANT: The spatial-correlation p-value above is parametric, not P_spin.\n');
fprintf(['For a true spin test, use only cortical parcels with the matching cortical ' ...
         'parcellation and sphere. Do not spin Tian/subcortical parcels as if they were cortex.\n']);


%% ========================================================================
% Local functions
% ========================================================================
function [gradient_all, explained_variance, gm_group, gm_all] = ...
    run_gradient_pipeline(fmri_cell_valid, sparsity_value, n_components, random_state)
% Compute a common group reference and align ALL subjects together.
%
% BrainSpace sparsity definition (MATLAB): percentage of smallest elements
% zeroed per row. Thus sparsity=90 retains approximately the largest 10%.

    nsubs = numel(fmri_cell_valid);
    nroi = size(fmri_cell_valid{1}, 1);

    % Mean connectivity matrix in the same validated sample
    mean_connectivity = zeros(nroi, nroi);
    for i = 1:nsubs
        mean_connectivity = mean_connectivity + fmri_cell_valid{i};
    end
    mean_connectivity = mean_connectivity / nsubs;
    mean_connectivity = (mean_connectivity + mean_connectivity') / 2;
    mean_connectivity(1:nroi+1:end) = 0;

    % Common reference gradient
    gm_group = GradientMaps('kernel','cs', ...
                            'approach','dm', ...
                            'n_components',n_components, ...
                            'random_state',random_state);
    gm_group = gm_group.fit(mean_connectivity, 'sparsity', sparsity_value);

    % Fit and align ALL individuals jointly to the same initial reference
    gm_all = GradientMaps('kernel','cs', ...
                          'approach','dm', ...
                          'alignment','pa', ...
                          'n_components',n_components, ...
                          'random_state',random_state);

    gm_all = gm_all.fit(fmri_cell_valid, ...
                        'sparsity',sparsity_value, ...
                        'reference',gm_group.gradients{1});

    gradient_all = nan(nsubs, nroi);
    explained_variance = nan(nsubs, n_components);

    for i = 1:nsubs
        % Use aligned first gradient for group comparisons
        gradient_all(i,:) = gm_all.aligned{i}(:,1)';

        % Eigenvalue-based variance proportions
        lambda_i = gm_all.lambda{i};
        lambda_i = lambda_i(:);
        if sum(lambda_i) ~= 0
            prop = lambda_i ./ sum(lambda_i);
        else
            prop = nan(size(lambda_i));
        end
        n_take = min(n_components, numel(prop));
        explained_variance(i,1:n_take) = prop(1:n_take)';
    end
end


function [beta_group, t_group, p_group, q_group] = ...
    regional_group_stats(gradient_all, age, education, group)
% Regional linear model:
%   G1(region) ~ Age + Education + Group
% HCs is the reference category; positive beta/t means PMI > HCs.

    nroi = size(gradient_all, 2);

    Group = categorical(group);
    Group = reordercats(Group, {'HCs','PMI'});

    beta_group = nan(nroi,1);
    t_group = nan(nroi,1);
    p_group = nan(nroi,1);

    for region = 1:nroi
        Y = gradient_all(:,region);
        tbl = table(Y, age, education, Group, ...
            'VariableNames', {'Gradient','Age','Education','Group'});

        lm = fitlm(tbl, 'Gradient ~ Age + Education + Group');
        coef_names = string(lm.CoefficientNames);
        idx_group = find(contains(coef_names, 'Group_PMI'), 1);

        if isempty(idx_group)
            error('Could not identify Group_PMI coefficient at region %d.', region);
        end

        beta_group(region) = lm.Coefficients.Estimate(idx_group);
        t_group(region) = lm.Coefficients.tStat(idx_group);
        p_group(region) = lm.Coefficients.pValue(idx_group);
    end

    % Benjamini-Hochberg FDR correction
    q_group = mafdr(p_group, 'BHFDR', true);
end
