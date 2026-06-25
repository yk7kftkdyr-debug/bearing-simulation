clear; clc;

projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
addpath(fullfile(projectRoot, '球轴承程序'));
addpath(fullfile(projectRoot, '滚子轴承程序'));

outDir = fullfile(projectRoot, 'results_micro_interface_sweep');
figDir = fullfile(outDir, 'figures');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

% Lightweight conference-paper sweep defaults.
% Texture depth uses 2e-6 m as a small deterministic surface texture
% perturbation. It is only a baseline comparison value, not a calibrated
% texture design parameter.
textureDepthDefault = 2e-6;
textureWidthDefault = 2e-4;
textureDensityDefault = 0.35;

cases = build_sweep_cases(textureDepthDefault, textureWidthDefault, textureDensityDefault);
bearingTypes = ["ball", "roller"];
records = repmat(empty_record(), 0, 1);

for iCase = 1:numel(cases)
    for iBearing = 1:numel(bearingTypes)
        fprintf('Running %s / %s = %s\n', bearingTypes(iBearing), cases(iCase).case_type, cases(iCase).case_value);
        record = run_one_case(projectRoot, bearingTypes(iBearing), cases(iCase));
        records(end + 1, 1) = record; %#ok<SAGROW>
        summaryTable = struct2table(records);
        writetable(summaryTable, fullfile(outDir, 'micro_interface_summary.csv'), 'Encoding', 'UTF-8');
    end
end

summaryTable = struct2table(records);
summaryPath = fullfile(outDir, 'micro_interface_summary.csv');
writetable(summaryTable, summaryPath, 'Encoding', 'UTF-8');

make_summary_figures(summaryTable, figDir);
write_markdown_report(summaryTable, outDir, textureDepthDefault, textureWidthDefault, textureDensityDefault);

disp('MICRO_INTERFACE_SWEEP_DONE');
disp(summaryPath);

function cases = build_sweep_cases(textureDepthDefault, textureWidthDefault, textureDensityDefault)
cases = repmat(struct('case_type', "", 'case_value', "", 'micro_config', []), 0, 1);

cfg = make_micro_interface_config();
cases(end + 1) = struct('case_type', "no_micro_effect", 'case_value', "baseline", 'micro_config', cfg);

temperatures = [20, 50, 80, 100];
for t = temperatures
    cfg = make_micro_interface_config();
    cfg.thermal.enabled = true;
    cfg.thermal.temperature_C = t;
    cfg.thermal.viscosity_factor = max(0.55, 1 - 0.006 * (t - 20));
    cfg.thermal.thermal_clearance_factor = 1 + 0.0008 * (t - 20);
    cfg.thermal.Ct1 = cfg.thermal.viscosity_factor;
    cfg.thermal.Ct2 = cfg.thermal.viscosity_factor;
    cfg.lubricant.enabled = true;
    cfg.lubricant.ehl_correction_factor = cfg.thermal.viscosity_factor;
    cases(end + 1) = struct('case_type', "temperature_C", 'case_value', string(t), 'micro_config', cfg); %#ok<AGROW>
end

textureCases = ["none", "longitudinal", "transverse"];
for textureCase = textureCases
    cfg = make_micro_interface_config();
    cfg.texture.texture_type = char(textureCase);
    cfg.texture.texture_depth = textureDepthDefault;
    cfg.texture.texture_width = textureWidthDefault;
    cfg.texture.texture_density = textureDensityDefault;
    if textureCase == "none"
        cfg.texture.enabled = false;
        cfg.texture.Cr = 1;
        cfg.texture.wenli = 6;
    elseif textureCase == "longitudinal"
        cfg.texture.enabled = true;
        cfg.texture.Cr = 1.08;
        cfg.texture.wenli = 1;
    else
        cfg.texture.enabled = true;
        cfg.texture.Cr = 0.94;
        cfg.texture.wenli = 2;
    end
    cases(end + 1) = struct('case_type', "texture_case", 'case_value', string(textureCase), 'micro_config', cfg); %#ok<AGROW>
end

debrisValues = [0, 5e-6, 10e-6, 20e-6];
for ud = debrisValues
    cfg = make_micro_interface_config();
    cfg.debris.enabled = ud > 0;
    cfg.debris.debris_displacement = ud;
    cfg.debris.ud = 0;
    cases(end + 1) = struct('case_type', "debris_displacement", 'case_value', sprintf('%.0e', ud), 'micro_config', cfg); %#ok<AGROW>
end
end

function record = run_one_case(projectRoot, bearingType, caseInfo)
record = empty_record();
record.bearing_type = char(bearingType);
record.case_type = char(caseInfo.case_type);
record.case_value = char(caseInfo.case_value);
record.comment = '';

if bearingType == "ball"
    bearingDir = fullfile(projectRoot, '球轴承程序');
    runner = @qiujieend;
else
    bearingDir = fullfile(projectRoot, '滚子轴承程序');
    runner = @qiujieall;
end

oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(bearingDir);
lastwarn('');
try
    runLog = evalc('[loadj, returndata] = runner([], caseInfo.micro_config);');
    [warningFlag, comment] = classify_warning(runLog);
    metrics = extract_metrics(bearingType, loadj, returndata);
    record = copy_metrics(record, metrics);
    record.warning_flag = warningFlag;
    record.comment = comment;
catch ME
    record.warning_flag = 'error';
    record.comment = sprintf('%s at %s:%d', ME.message, error_file(ME), error_line(ME));
end
end

function [warningFlag, comment] = classify_warning(runLog)
[lastMsg, lastId] = lastwarn;
warningFlag = 'none';
comment = '运行完成';

hasWarningText = contains(runLog, 'Warning', 'IgnoreCase', true) || contains(runLog, '警告');
hasSingular = contains(runLog, 'singular', 'IgnoreCase', true) || contains(lastMsg, 'singular', 'IgnoreCase', true) || ...
    contains(runLog, '奇异') || contains(lastMsg, '奇异');
if hasSingular
    warningFlag = 'singular_matrix';
    comment = '运行完成，但出现奇异矩阵相关警告';
elseif hasWarningText || strlength(string(lastMsg)) > 0 || strlength(string(lastId)) > 0
    warningFlag = 'warning';
    if strlength(string(lastMsg)) > 0
        comment = char(lastMsg);
    else
        comment = '运行完成，但命令输出包含 warning';
    end
end
end

function metrics = extract_metrics(bearingType, loadj, returndata)
metrics.loaded_element_count = loadj;
metrics.max_contact_load = max_from_mats({'Q1.mat', 'Q2.mat'});
metrics.max_contact_stress = max_from_mats({'Ph1.mat', 'Ph2.mat'});
metrics.min_oil_film_thickness = min_from_mats({'oilh1.mat', 'oilh2.mat'});
metrics.PV_value = max_from_mats({'pvzhi1.mat', 'pvzhi2.mat', 'pvzhinonload.mat', 'pvzhi1nonload.mat'});

if bearingType == "ball"
    metrics.cage_slip_ratio = value_or_nan(returndata, 6);
    metrics.equivalent_stiffness_x = value_or_nan(returndata, 1);
    metrics.equivalent_stiffness_y = value_or_nan(returndata, 2);
    metrics.life_or_damage_index = value_or_nan(returndata, 9);
else
    metrics.cage_slip_ratio = value_or_nan(returndata, 5);
    metrics.equivalent_stiffness_x = NaN;
    metrics.equivalent_stiffness_y = value_or_nan(returndata, 1);
    metrics.life_or_damage_index = value_or_nan(returndata, 6);
end
end

function record = copy_metrics(record, metrics)
names = fieldnames(metrics);
for i = 1:numel(names)
    record.(names{i}) = metrics.(names{i});
end
end

function value = max_from_mats(fileNames)
value = NaN;
for i = 1:numel(fileNames)
    if exist(fileNames{i}, 'file') == 2
        data = load(fileNames{i});
        names = fieldnames(data);
        for j = 1:numel(names)
            candidate = data.(names{j});
            if isnumeric(candidate) && ~isempty(candidate)
                value = max([value; candidate(:)], [], 'omitnan');
            end
        end
    end
end
end

function value = min_from_mats(fileNames)
value = NaN;
for i = 1:numel(fileNames)
    if exist(fileNames{i}, 'file') == 2
        data = load(fileNames{i});
        names = fieldnames(data);
        for j = 1:numel(names)
            candidate = data.(names{j});
            if isnumeric(candidate) && ~isempty(candidate)
                value = min([value; candidate(:)], [], 'omitnan');
            end
        end
    end
end
end

function value = value_or_nan(values, index)
if numel(values) >= index
    value = values(index);
else
    value = NaN;
end
end

function record = empty_record()
record = struct( ...
    'bearing_type', '', ...
    'case_type', '', ...
    'case_value', '', ...
    'loaded_element_count', NaN, ...
    'max_contact_load', NaN, ...
    'max_contact_stress', NaN, ...
    'min_oil_film_thickness', NaN, ...
    'PV_value', NaN, ...
    'cage_slip_ratio', NaN, ...
    'equivalent_stiffness_x', NaN, ...
    'equivalent_stiffness_y', NaN, ...
    'life_or_damage_index', NaN, ...
    'warning_flag', '', ...
    'comment', '');
end

function fileName = error_file(ME)
if isempty(ME.stack)
    fileName = 'unknown';
else
    fileName = ME.stack(1).file;
end
end

function lineNo = error_line(ME)
if isempty(ME.stack)
    lineNo = 0;
else
    lineNo = ME.stack(1).line;
end
end

function make_summary_figures(summaryTable, figDir)
set(groot, 'defaultFigureVisible', 'off');
set(groot, 'defaultAxesFontName', 'Arial Unicode MS');
set(groot, 'defaultTextFontName', 'Arial Unicode MS');

plot_case_metric(summaryTable, 'temperature_C', 'max_contact_stress', ...
    '温度效应对最大接触应力的影响', '温度 (degC)', '最大接触应力 (Pa)', ...
    fullfile(figDir, 'temperature_effect_contact_stress.png'));
plot_case_metric(summaryTable, 'temperature_C', 'min_oil_film_thickness', ...
    '温度效应对最小油膜厚度的影响', '温度 (degC)', '最小油膜厚度 (m)', ...
    fullfile(figDir, 'temperature_effect_oil_film.png'));
plot_case_metric(summaryTable, 'texture_case', 'min_oil_film_thickness', ...
    '表面纹理对最小油膜厚度的影响', '纹理类型', '最小油膜厚度 (m)', ...
    fullfile(figDir, 'texture_effect_oil_film.png'));
plot_case_metric(summaryTable, 'texture_case', 'max_contact_stress', ...
    '表面纹理对最大接触应力的影响', '纹理类型', '最大接触应力 (Pa)', ...
    fullfile(figDir, 'texture_effect_contact_stress.png'));
plot_case_metric(summaryTable, 'debris_displacement', 'max_contact_stress', ...
    '杂质扰动对最大接触应力的影响', '杂质位移 (m)', '最大接触应力 (Pa)', ...
    fullfile(figDir, 'debris_effect_contact_stress.png'));
plot_case_metric(summaryTable, 'debris_displacement', 'cage_slip_ratio', ...
    '杂质扰动对保持架打滑率的影响', '杂质位移 (m)', '保持架打滑率', ...
    fullfile(figDir, 'debris_effect_slip_ratio.png'));
plot_stiffness(summaryTable, fullfile(figDir, 'ball_vs_roller_stiffness_comparison.png'));
end

function plot_case_metric(summaryTable, caseType, metricName, titleText, xText, yText, outPath)
rows = strcmp(summaryTable.case_type, caseType);
data = summaryTable(rows, :);
figure('Color', 'w', 'Position', [100, 100, 760, 480]);
hold on; grid on; box on;
bearings = {'ball', 'roller'};
markers = {'-o', '-s'};
for i = 1:numel(bearings)
    bRows = strcmp(data.bearing_type, bearings{i});
    bData = data(bRows, :);
    [xValues, xLabels] = case_x_values(bData.case_value, caseType);
    yValues = bData.(metricName);
    plot(xValues, yValues, markers{i}, 'LineWidth', 1.6, 'MarkerSize', 6);
end
title(titleText);
xlabel(xText);
ylabel(yText);
legend({'球轴承', '滚子轴承'}, 'Location', 'best');
if ~isempty(xLabels)
    xticks(1:numel(xLabels));
    xticklabels(xLabels);
end
exportgraphics(gcf, outPath, 'Resolution', 220);
close(gcf);
end

function [xValues, xLabels] = case_x_values(caseValues, caseType)
xLabels = {};
if strcmp(caseType, 'texture_case')
    xValues = 1:numel(caseValues);
    xLabels = cellstr(caseValues);
elseif strcmp(caseType, 'debris_displacement')
    xValues = zeros(height(table(caseValues)), 1);
    for i = 1:numel(caseValues)
        xValues(i) = sscanf(caseValues{i}, '%e');
    end
else
    xValues = str2double(caseValues);
end
end

function plot_stiffness(summaryTable, outPath)
rows = strcmp(summaryTable.case_type, 'no_micro_effect');
data = summaryTable(rows, :);
figure('Color', 'w', 'Position', [100, 100, 760, 480]);
labels = categorical({'球轴承 k_x', '球轴承 k_y', '滚子轴承 k_y'});
values = [
    data.equivalent_stiffness_x(strcmp(data.bearing_type, 'ball'));
    data.equivalent_stiffness_y(strcmp(data.bearing_type, 'ball'));
    data.equivalent_stiffness_y(strcmp(data.bearing_type, 'roller'))];
bar(labels, values);
grid on; box on;
title('球轴承与滚子轴承基准等效刚度对比');
ylabel('等效刚度 (N/m)');
exportgraphics(gcf, outPath, 'Resolution', 220);
close(gcf);
end

function write_markdown_report(summaryTable, outDir, textureDepthDefault, textureWidthDefault, textureDensityDefault)
reportPath = fullfile(outDir, 'micro_interface_sweep_report.md');
fid = fopen(reportPath, 'w', 'n', 'UTF-8');
cleanupObj = onCleanup(@() fclose(fid));

ballBase = summaryTable(strcmp(summaryTable.bearing_type, 'ball') & strcmp(summaryTable.case_type, 'no_micro_effect'), :);
rollerBase = summaryTable(strcmp(summaryTable.bearing_type, 'roller') & strcmp(summaryTable.case_type, 'no_micro_effect'), :);
warningRows = summaryTable(~strcmp(summaryTable.warning_flag, 'none'), :);

fprintf(fid, '# 微观界面因素轻量批处理扫描报告\n\n');
fprintf(fid, '## 1. 批处理工况说明\n\n');
fprintf(fid, '本次扫描覆盖球轴承和滚子轴承两个默认算例。工况包括基准无扰动、温度 20/50/80/100 degC、纹理 none/longitudinal/transverse、杂质位移 0/5e-6/10e-6/20e-6 m。\n\n');
fprintf(fid, '纹理默认深度为 %.3g m，宽度为 %.3g m，密度为 %.2f；该值用于会议论文阶段的轻量对比，不代表最终优化纹理参数。\n\n', textureDepthDefault, textureWidthDefault, textureDensityDefault);

fprintf(fid, '## 2. 球轴承与滚子轴承基准结果\n\n');
write_base_line(fid, '球轴承', ballBase);
write_base_line(fid, '滚子轴承', rollerBase);

fprintf(fid, '\n## 3. 热效应影响规律\n\n');
write_trend_line(fid, summaryTable, 'temperature_C', 'max_contact_stress', '最大接触应力');
write_trend_line(fid, summaryTable, 'temperature_C', 'min_oil_film_thickness', '最小油膜厚度');

fprintf(fid, '\n## 4. 纹理影响规律\n\n');
write_trend_line(fid, summaryTable, 'texture_case', 'min_oil_film_thickness', '最小油膜厚度');
write_trend_line(fid, summaryTable, 'texture_case', 'max_contact_stress', '最大接触应力');

fprintf(fid, '\n## 5. 杂质扰动影响规律\n\n');
write_trend_line(fid, summaryTable, 'debris_displacement', 'max_contact_stress', '最大接触应力');
write_trend_line(fid, summaryTable, 'debris_displacement', 'cage_slip_ratio', '保持架打滑率');

fprintf(fid, '\n## 6. 对核心动力学指标的影响\n\n');
fprintf(fid, '本次汇总表给出了接触载荷、接触应力、油膜厚度、PV 值、打滑率、等效刚度与寿命/损伤指数。温度和纹理主要通过油膜修正进入轻量对比，杂质扰动通过接触变形附加位移进入，因此异常工况需要优先检查收敛性和接触变形合理性。\n\n');

fprintf(fid, '## 7. 可用于会议文章的结果\n\n');
fprintf(fid, '可直接用于会议文章的内容包括：基准球/滚子轴承对比、温度对接触应力和油膜厚度的影响图、纹理方向对油膜和接触应力的对比图、杂质位移对接触应力和打滑率的影响图、以及球/滚子轴承基准等效刚度对比。\n\n');

fprintf(fid, '## 8. 仍需验证或改进的问题\n\n');
if isempty(warningRows)
    fprintf(fid, '本次扫描未记录 warning。后续仍建议补充实验或高保真模型校核热修正系数、纹理等效系数和杂质位移等效方式。\n');
else
    fprintf(fid, '本次扫描记录到以下 warning/异常，会议论文使用前需要复核：\n\n');
    for i = 1:height(warningRows)
        fprintf(fid, '- %s / %s=%s：%s，%s\n', warningRows.bearing_type{i}, warningRows.case_type{i}, warningRows.case_value{i}, warningRows.warning_flag{i}, warningRows.comment{i});
    end
end
end

function write_base_line(fid, name, row)
if isempty(row)
    fprintf(fid, '- %s：未得到基准结果。\n', name);
else
    fprintf(fid, '- %s：承载滚动体数 %.0f，最大接触载荷 %.4g N，最大接触应力 %.4g Pa，最小油膜厚度 %.4g m，打滑率 %.4g，等效刚度 y %.4g N/m。\n', ...
        name, row.loaded_element_count(1), row.max_contact_load(1), row.max_contact_stress(1), row.min_oil_film_thickness(1), row.cage_slip_ratio(1), row.equivalent_stiffness_y(1));
end
end

function write_trend_line(fid, summaryTable, caseType, metricName, metricLabel)
rows = strcmp(summaryTable.case_type, caseType);
data = summaryTable(rows, :);
bearings = {'ball', 'roller'};
for i = 1:numel(bearings)
    bData = data(strcmp(data.bearing_type, bearings{i}), :);
    bearingLabel = "球轴承";
    if strcmp(bearings{i}, 'roller')
        bearingLabel = "滚子轴承";
    end
    values = bData.(metricName);
    if isempty(values) || all(isnan(values))
        fprintf(fid, '- %s %s：无有效数据。\n', bearingLabel, metricLabel);
    else
        [minValue, minIndex] = min(values, [], 'omitnan');
        [maxValue, maxIndex] = max(values, [], 'omitnan');
        fprintf(fid, '- %s %s：范围 %.4g 至 %.4g，对应工况 %s 至 %s。\n', bearingLabel, metricLabel, minValue, maxValue, bData.case_value{minIndex}, bData.case_value{maxIndex});
    end
end
end
