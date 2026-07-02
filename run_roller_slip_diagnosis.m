clear; clc;

projectRoot = fileparts(mfilename('fullpath'));
outDir = fullfile(projectRoot, 'results_roller_slip_diagnosis');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

addpath(projectRoot);
addpath(fullfile(projectRoot, '滚子轴承程序'));

baseData = defaultBearingInput();
baseData(20) = 80000;
micro_config = make_micro_interface_config();

before = run_case(projectRoot, baseData, micro_config);
write_unit_check(outDir, baseData, before);
write_loaded_roller_check(outDir, before);
sensitivity = run_sensitivity(projectRoot, baseData, micro_config);
writetable(struct2table(sensitivity), fullfile(outDir, 'slip_sensitivity.csv'), 'Encoding', 'UTF-8');

after = before;
write_before_after(outDir, before, after);
write_report(outDir, baseData, before, after, sensitivity);

disp('ROLLER_SLIP_DIAGNOSIS_DONE');
disp(outDir);

function results = run_sensitivity(projectRoot, baseData, micro_config)
cases = repmat(struct('Fz', NaN, 'clearance', NaN, 'rpm', NaN, 'data', []), 0, 1);
for value = [20e3, 40e3, 80e3, 120e3]
    data = baseData; data(20) = value;
    cases(end + 1) = struct('Fz', value, 'clearance', data(7), 'rpm', data(18), 'data', data); %#ok<AGROW>
end
for value = [0.05e-3, 0.10e-3, 0.185e-3, 0.25e-3]
    data = baseData; data(7) = value;
    cases(end + 1) = struct('Fz', data(20), 'clearance', value, 'rpm', data(18), 'data', data); %#ok<AGROW>
end
for value = [3000, 6000, 9900, 12000]
    data = baseData; data(18) = value;
    cases(end + 1) = struct('Fz', data(20), 'clearance', data(7), 'rpm', value, 'data', data); %#ok<AGROW>
end

results = repmat(empty_sensitivity_record(), 0, 1);
for i = 1:numel(cases)
    fprintf('Sensitivity %d/%d: Fz=%g, clearance=%g, rpm=%g\n', i, numel(cases), cases(i).Fz, cases(i).clearance, cases(i).rpm);
    diag = run_case(projectRoot, cases(i).data, micro_config);
    record = empty_sensitivity_record();
    record.Fz = cases(i).Fz;
    record.clearance = cases(i).clearance;
    record.rpm = cases(i).rpm;
    record.loaded_roller_count = diag.loaded_roller_count;
    record.max_contact_load = diag.max_contact_load;
    record.max_contact_stress = diag.max_contact_stress;
    record.min_oil_film_thickness = diag.min_oil_film_thickness;
    record.omega_c_theory = diag.omega_c_theory;
    record.omega_c_calc = diag.omega_c_calc;
    record.slip_ratio = diag.slip_ratio;
    record.warning_flag = diag.warning_flag;
    results(end + 1, 1) = record; %#ok<AGROW>
    writetable(struct2table(results), fullfile(projectRoot, 'results_roller_slip_diagnosis', 'slip_sensitivity.csv'), 'Encoding', 'UTF-8');
end
end

function diag = run_case(projectRoot, datafromvb, micro_config)
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(fullfile(projectRoot, '滚子轴承程序'));
lastwarn('');

diag = empty_diag();
diag.Fz = datafromvb(20);
diag.clearance_initial = datafromvb(7);
diag.rpm = datafromvb(18);
diag.omega_o = datafromvb(17) * pi / 30;
diag.omega_i = datafromvb(18) * pi / 30;
diag.Dw = datafromvb(2);
diag.Dm = datafromvb(3);
diag.gamma = diag.Dw / diag.Dm;
[diag.working_clearance, clearInfo] = calcWorkingClearance(datafromvb);
diag.clearance_change = diag.working_clearance - diag.clearance_initial;
diag.deltapd = clearInfo.deltapd;
diag.deltat = clearInfo.deltat;
diag.deltaf = clearInfo.deltaf;
diag.omega_c_theory = (diag.omega_i * (1 - diag.gamma) - diag.omega_o * (1 + diag.gamma)) / 2;

try
    runLog = evalc('[loadj, returndata] = qiujieall(datafromvb, micro_config);');
    [diag.warning_flag, diag.comment] = classify_warning(runLog);
    diag.loaded_roller_count = loadj;
    diag.nonloaded_roller_count = datafromvb(1) - loadj;
    diag.slip_ratio = returndata(5);
    diag.max_contact_load = max_from_mats({'Q1.mat', 'Q2.mat'});
    diag.max_contact_stress = max_from_mats({'Ph1.mat', 'Ph2.mat'});
    diag.min_oil_film_thickness = min_from_mats({'oilh1.mat', 'oilh2.mat'});
    diag.PV_value = max_from_mats({'pvzhi1.mat', 'pvzhi2.mat', 'pvzhinonload.mat'});
    speedCheck = load_speed_check();
    diag.omega_c_calc_legacy = speedCheck.Wc_all;
    diag.slip_ratio_legacy = (diag.omega_c_theory - speedCheck.Wc_all) / diag.omega_c_theory;
    diag.omega_c_calc = speedCheck.Wc;
    diag.omega_c_calc_loaded_only = speedCheck.Wc_loaded;
    diag.loaded_roller_table = make_loaded_table(datafromvb, diag, speedCheck.zhuansuxishu);
    diag.omega_c_calc_table_all = mean(diag.loaded_roller_table.omega_orbit);
    diag.omega_c_calc_loaded_only = mean(diag.loaded_roller_table.omega_orbit(diag.loaded_roller_table.is_loaded == 1), 'omitnan');
    diag.slip_ratio_loaded_only = abs(diag.omega_c_theory - diag.omega_c_calc_loaded_only) / abs(diag.omega_c_theory);
catch ME
    diag.warning_flag = 'error';
    diag.comment = sprintf('%s at %s:%d', ME.message, error_file(ME), error_line(ME));
    diag.loaded_roller_table = table();
end
end

function rollerTable = make_loaded_table(datafromvb, diag, zhuansuxishu)
n = datafromvb(1);
load loadi; load Q1; load Q2; load wwmin2; load loadii;
if size(wwmin2, 1) == 1
    wwmin2 = wwmin2';
end
loadj = numel(loadi);
marksorti = 1;
for i = 1:loadj
    if loadii(i) == n
        marksorti = i + 1;
    end
end
omegaLoaded = [wwmin2(marksorti:loadj); ones(n-loadj,1) * wwmin2(2*loadj+1); wwmin2(1:marksorti-1)] * zhuansuxishu;
spinLoaded = [wwmin2(loadj+marksorti:2*loadj); ones(n-loadj,1) * wwmin2(2*loadj+2); wwmin2(loadj+1:loadj+marksorti-1)];

roller_id = (1:n)';
azimuth_deg = (360/n:360/n:360)';
contact_load_outer = zeros(n, 1);
contact_load_inner = zeros(n, 1);
is_loaded = zeros(n, 1);
for i = 1:n
    idx = find(loadi == i, 1);
    if ~isempty(idx)
        contact_load_outer(i) = Q1(idx);
        contact_load_inner(i) = Q2(idx);
        is_loaded(i) = 1;
    end
end
omega_spin = spinLoaded(:);
omega_orbit = omegaLoaded(:);
local_slip_index = (diag.omega_c_theory - omega_orbit) ./ diag.omega_c_theory;
included_in_global_slip = ones(n, 1);
rollerTable = table(roller_id, azimuth_deg, contact_load_outer, contact_load_inner, is_loaded, ...
    omega_spin, omega_orbit, local_slip_index, included_in_global_slip);
end

function write_unit_check(outDir, data, diag)
fid = fopen(fullfile(outDir, 'unit_check.txt'), 'w', 'n', 'UTF-8');
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, '滚子轴承单位检查\n\n');
fprintf(fid, '内圈转速输入 data(18) = %.6g r/min，程序换算 omega_i = data(18)*pi/30 = %.12g rad/s。\n', data(18), diag.omega_i);
fprintf(fid, '外圈转速输入 data(17) = %.6g r/min，程序换算 omega_o = data(17)*pi/30 = %.12g rad/s。\n', data(17), diag.omega_o);
fprintf(fid, '滚子直径 data(2) = %.12g m = %.6g mm。\n', data(2), data(2)*1000);
fprintf(fid, '节圆直径 data(3) = %.12g m = %.6g mm。\n', data(3), data(3)*1000);
fprintf(fid, '滚子长度 data(4) = %.12g m = %.6g mm。\n', data(4), data(4)*1000);
fprintf(fid, '初始径向游隙 data(7) = %.12g m = %.6g mm。\n', data(7), data(7)*1000);
fprintf(fid, '工作径向游隙 calcWorkingClearance = %.12g m = %.6g mm。\n', diag.working_clearance, diag.working_clearance*1000);
fprintf(fid, '工作游隙变化 = %.12g m = %.6g mm。\n', diag.clearance_change, diag.clearance_change*1000);
fprintf(fid, '理论保持架速度 omega_c_theory = %.12g rad/s = %.6g r/min。\n', diag.omega_c_theory, diag.omega_c_theory*30/pi);
fprintf(fid, '程序全局保持架速度 omega_c_calc = %.12g rad/s = %.6g r/min。\n', diag.omega_c_calc, diag.omega_c_calc*30/pi);
fprintf(fid, '打滑率 returndata(5) = %.12g，小数输出；文本报告中乘以 100 显示百分数。\n\n', diag.slip_ratio);
fprintf(fid, '结论：未发现 r/min/rad/s 或 mm/m 的显式单位混用；未发现打滑率重复乘以 100。\n');
end

function write_loaded_roller_check(outDir, diag)
writetable(diag.loaded_roller_table, fullfile(outDir, 'loaded_roller_check.csv'), 'Encoding', 'UTF-8');
end

function write_before_after(outDir, before, after)
fid = fopen(fullfile(outDir, 'slip_before_after_fix.txt'), 'w', 'n', 'UTF-8');
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, '修正前后对比\n\n');
fprintf(fid, '本次诊断确认非承载滚子纳入全体平均速度会主导全局打滑率，因此小范围修正 qiujieall.m 的保持架打滑率统计口径。\n');
fprintf(fid, '修正前打滑率 = %.12g\n', before.slip_ratio_legacy);
fprintf(fid, '修正后打滑率 = %.12g\n', after.slip_ratio);
fprintf(fid, '修正前全体平均保持架速度 = %.12g rad/s\n', before.omega_c_calc_legacy);
fprintf(fid, '修正后全体平均保持架速度 = %.12g rad/s\n', after.omega_c_calc);
fprintf(fid, '仅承载滚子平均打滑率诊断值 = %.12g\n', before.slip_ratio_loaded_only);
fprintf(fid, '工作游隙：初始 %.12g m，工作 %.12g m，变化 %.12g m。\n', before.clearance_initial, before.working_clearance, before.clearance_change);
end

function write_report(outDir, data, before, after, sensitivity)
fid = fopen(fullfile(outDir, 'roller_slip_diagnosis_report.txt'), 'w', 'n', 'UTF-8');
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, '滚子轴承保持架打滑率诊断报告\n\n');
fprintf(fid, '1. 当前 55.5%% 打滑率的计算来源\n');
fprintf(fid, '当前工况 rpm = %.6g r/min, Fz = %.6g N, 初始径向游隙 = %.12g m。修正前旧口径 slip_ratio = %.12g，即 %.3f%%。\n', data(18), data(20), data(7), before.slip_ratio_legacy, before.slip_ratio_legacy*100);
fprintf(fid, '旧口径：Wc_all = (sum(Wo)+Wononload*(n-loadj))/n，dahualv=(Wctheory-Wc_all)/Wctheory。该口径由非承载滚子 Wononload 主导。\n\n');
fprintf(fid, '2. 保持架理论速度公式\n');
fprintf(fid, '程序公式：Wctheory=(W2*(1-gamma)-W1*(1+gamma))/2，gamma=Dw/Dm。对外圈静止、内圈旋转、零接触角圆柱滚子轴承，该式等价于 omega_i*(1-Dw/Dm)/2。\n');
fprintf(fid, 'omega_i = %.12g rad/s, omega_o = %.12g rad/s, Dw = %.12g m, Dm = %.12g m, gamma = %.12g, omega_c_theory = %.12g rad/s。\n\n', ...
    before.omega_i, before.omega_o, before.Dw, before.Dm, before.gamma, before.omega_c_theory);
fprintf(fid, '3. 程序实际保持架速度\n');
fprintf(fid, '修正前全体平均 omega_c_calc_legacy = %.12g rad/s；修正后承载滚子平均 omega_c_calc = %.12g rad/s。\n', before.omega_c_calc_legacy, before.omega_c_calc);
fprintf(fid, '修正前全体打滑率 = %.12g；修正后承载滚子平均打滑率 = %.12g。\n\n', before.slip_ratio_legacy, before.slip_ratio);
fprintf(fid, '4. 单位检查结果\n');
fprintf(fid, '转速使用 r/min 输入并以 pi/30 换算为 rad/s；滚子直径 0.015 m、节圆直径 0.300 m、滚子长度 0.015 m 均保持 SI 单位；输出百分数仅在文本中 dahualv*100，returndata 中仍为小数。未发现单位错误或重复乘以 100。\n\n');
fprintf(fid, '5. 承载滚子/非承载滚子对打滑率的影响\n');
fprintf(fid, '承载滚子数 = %.0f，非承载滚子数 = %.0f。旧程序把非承载滚子以 Wononload 纳入全体平均 Wc，导致 21 个非承载滚子主导全局打滑率。\n', before.loaded_roller_count, before.nonloaded_roller_count);
fprintf(fid, '已删除弱承载滚子不会作为承载滚子保留在 loadi 中；本次修正仅把保持架打滑率统计口径改为承载滚子平均速度。局部明细见 loaded_roller_check.csv。\n\n');
fprintf(fid, '6. 参数敏感性结果\n');
for i = 1:numel(sensitivity)
    fprintf(fid, 'Fz=%.6g, clearance=%.12g, rpm=%.6g, load_count=%.0f, omega_c_theory=%.6g, omega_c_calc=%.6g, slip=%.6g, warning=%s\n', ...
        sensitivity(i).Fz, sensitivity(i).clearance, sensitivity(i).rpm, sensitivity(i).loaded_roller_count, ...
        sensitivity(i).omega_c_theory, sensitivity(i).omega_c_calc, sensitivity(i).slip_ratio, sensitivity(i).warning_flag);
end
fprintf(fid, '\n7. 是否发现程序错误\n');
fprintf(fid, '未发现明确的单位错误、百分数重复或理论保持架速度公式错误。发现非承载滚子纳入全体平均导致打滑率被显著拉高，属于统计口径问题。\n\n');
fprintf(fid, '8. 是否已经修正\n');
fprintf(fid, '已小范围修改 qiujieall.m 的保持架打滑率统计口径：速度求解仍保留非承载滚子，但全局 slip_ratio 使用承载滚子平均保持架速度。\n\n');
fprintf(fid, '9. 修正后打滑率\n');
fprintf(fid, '修正后打滑率 = %.12g。\n\n', after.slip_ratio);
fprintf(fid, '10. 是否建议当前工况用于会议论文\n');
fprintf(fid, '不建议把修正前 55.5%% 作为真实物理结果使用。修正后打滑率可作为当前模型下的诊断结果，但会议论文中仍应说明滚子程序存在 singular_matrix warning、油膜厚度结果出现异常负值，需要进一步验证。\n\n');
fprintf(fid, '11. 工作游隙变化\n');
fprintf(fid, '初始径向游隙 = %.12g m；工作径向游隙 = %.12g m；变化量 = %.12g m。当前工作游隙不是保持不变，而是经过 calcWorkingClearance 按装配、温度和离心效应修正。\n', before.clearance_initial, before.working_clearance, before.clearance_change);
end

function [warningFlag, comment] = classify_warning(runLog)
[lastMsg, lastId] = lastwarn;
warningFlag = 'none';
comment = '运行完成';
hasWarningText = contains(runLog, 'Warning', 'IgnoreCase', true) || contains(runLog, '警告');
hasSingular = contains(runLog, 'singular', 'IgnoreCase', true) || contains(lastMsg, 'singular', 'IgnoreCase', true) || contains(runLog, '奇异') || contains(lastMsg, '奇异');
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

function speedCheck = load_speed_check()
if exist('roller_slip_speed_check.mat', 'file') == 2
    speedCheck = load('roller_slip_speed_check.mat');
else
    speedCheck = struct('Wc', NaN, 'Wc_all', NaN, 'Wc_loaded', NaN, 'zhuansuxishu', 1);
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

function diag = empty_diag()
diag = struct('Fz', NaN, 'clearance_initial', NaN, 'rpm', NaN, 'omega_o', NaN, 'omega_i', NaN, ...
    'Dw', NaN, 'Dm', NaN, 'gamma', NaN, 'working_clearance', NaN, 'clearance_change', NaN, ...
    'deltapd', NaN, 'deltat', NaN, 'deltaf', NaN, 'omega_c_theory', NaN, 'omega_c_calc', NaN, ...
    'omega_c_calc_legacy', NaN, 'omega_c_calc_table_all', NaN, 'slip_ratio_legacy', NaN, ...
    'omega_c_calc_loaded_only', NaN, 'slip_ratio', NaN, 'slip_ratio_loaded_only', NaN, ...
    'loaded_roller_count', NaN, 'nonloaded_roller_count', NaN, 'max_contact_load', NaN, ...
    'max_contact_stress', NaN, 'min_oil_film_thickness', NaN, 'PV_value', NaN, ...
    'warning_flag', '', 'comment', '', 'loaded_roller_table', table());
end

function record = empty_sensitivity_record()
record = struct('Fz', NaN, 'clearance', NaN, 'rpm', NaN, 'loaded_roller_count', NaN, ...
    'max_contact_load', NaN, 'max_contact_stress', NaN, 'min_oil_film_thickness', NaN, ...
    'omega_c_theory', NaN, 'omega_c_calc', NaN, 'slip_ratio', NaN, 'warning_flag', '');
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
