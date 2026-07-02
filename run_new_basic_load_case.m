clear; clc;

projectRoot = fileparts(mfilename('fullpath'));
desktopOut = fullfile('/Users/bailingguan/Desktop', '轴承新基本工况_球径向20000_滚子径向80000');
ballOut = fullfile(desktopOut, 'ball');
rollerOut = fullfile(desktopOut, 'roller');
if exist(desktopOut, 'dir') ~= 7
    mkdir(desktopOut);
end
if exist(ballOut, 'dir') ~= 7
    mkdir(ballOut);
end
if exist(rollerOut, 'dir') ~= 7
    mkdir(rollerOut);
end

addpath(projectRoot);
addpath(fullfile(projectRoot, '球轴承程序'));
addpath(fullfile(projectRoot, '滚子轴承程序'));

micro_config = make_micro_interface_config();
records = repmat(empty_record(), 0, 1);

ballData = default_ball_input();
ballData(16) = 20000;
ballData(17) = 0;
ballData(18) = 20000;
records(end + 1, 1) = run_ball_case(projectRoot, ballData, micro_config, ballOut);

rollerData = defaultBearingInput();
rollerData(19) = 0;
rollerData(20) = 80000;
records(end + 1, 1) = run_roller_case(projectRoot, rollerData, micro_config, rollerOut);

summaryTable = struct2table(records);
writetable(summaryTable, fullfile(desktopOut, 'new_basic_load_case_summary.csv'), 'Encoding', 'UTF-8');

fid = fopen(fullfile(desktopOut, 'new_basic_load_case_report.txt'), 'w', 'n', 'UTF-8');
fprintf(fid, '新基本工况运行结果\n\n');
fprintf(fid, '球轴承：Fx = 20000 N, Fy = 0 N, Fz = 20000 N；其余参数沿用 qiujieend.m 默认值；micro_config 为无扰动基准。\n');
fprintf(fid, '滚子轴承：Fy = 0 N, Fz = 80000 N；其余参数沿用 defaultBearingInput.m 默认值；micro_config 为无扰动基准。\n\n');
for i = 1:height(summaryTable)
    fprintf(fid, '%s：loaded_element_count = %.0f, max_contact_load = %.6g N, max_contact_stress = %.6g Pa, min_oil_film_thickness = %.6g m, cage_slip_ratio = %.6g, warning_flag = %s, comment = %s\n', ...
        summaryTable.bearing_type{i}, summaryTable.loaded_element_count(i), summaryTable.max_contact_load(i), ...
        summaryTable.max_contact_stress(i), summaryTable.min_oil_film_thickness(i), summaryTable.cage_slip_ratio(i), ...
        summaryTable.warning_flag{i}, summaryTable.comment{i});
end
fclose(fid);

detailPath = fullfile(desktopOut, '当前基本工况详细说明报告.txt');
fid = fopen(detailPath, 'w', 'n', 'UTF-8');
fprintf(fid, '当前轴承基本工况详细说明报告\n');
fprintf(fid, '生成时间：%s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, '一、工况来源\n');
fprintf(fid, '本次只运行基本工况，不进行参数扫描。球轴承参数沿用 qiujieend.m 默认输入，仅将径向载荷 Fz 改为 20000 N；滚子轴承参数沿用 defaultBearingInput.m 默认输入，仅将径向载荷 Fz 改为 80000 N。\n');
fprintf(fid, 'micro_config 使用无扰动基准：thermal.enabled=false, texture.enabled=false, debris.enabled=false。\n\n');
fprintf(fid, '二、球轴承基本工况\n');
fprintf(fid, '滚动体数：15\n');
fprintf(fid, '球直径：22.23 mm\n');
fprintf(fid, '节圆直径：125.3 mm\n');
fprintf(fid, '内/外沟曲率系数：0.5232 / 0.5232\n');
fprintf(fid, '初始接触角：40 deg\n');
fprintf(fid, '内圈转速：10000 r/min，外圈转速：0 r/min\n');
fprintf(fid, '载荷：Fx=20000 N，Fy=0 N，Fz=20000 N\n');
fprintf(fid, '力矩：My=0 N*m，Mz=0 N*m\n');
fprintf(fid, '粗糙度：内圈/外圈/滚动体均为 1e-7 m\n');
fprintf(fid, '润滑油密度：970 kg/m3，黏度：0.0318 Pa*s，黏压系数：1.28e-8\n');
fprintf(fid, '温度输入：油温 20 degC；其他热装配参数保持默认。\n\n');
fprintf(fid, '三、滚子轴承基本工况\n');
fprintf(fid, '滚子数：30\n');
fprintf(fid, '滚子直径：15 mm\n');
fprintf(fid, '节圆直径：300 mm\n');
fprintf(fid, '滚子长度：15 mm\n');
fprintf(fid, '有效接触线长：15 mm\n');
fprintf(fid, '初始径向游隙：0.15 mm\n');
fprintf(fid, '内圈转速：9900 r/min，外圈转速：0 r/min\n');
fprintf(fid, '载荷：Fy=0 N，Fz=80000 N\n');
fprintf(fid, '力矩：My=0 N*m，Mz=0 N*m\n');
fprintf(fid, '粗糙度：内圈/外圈/滚子均为 1e-7 m\n');
fprintf(fid, '润滑油密度：885 kg/m3，黏度：0.0318 Pa*s，黏压系数：1.25e-8\n');
fprintf(fid, '温度输入：各部件默认 20 degC。\n\n');
fprintf(fid, '四、滚子轴承打滑率说明\n');
fprintf(fid, '当前滚子轴承程序已采用修正后的保持架打滑率统计口径：速度求解仍保留承载与非承载滚子，但全局打滑率使用承载滚子平均公转速度与理论保持架速度比较。\n');
fprintf(fid, '该修正用于避免非承载滚子统一低公转速度主导整体平均，从而把打滑率异常放大。\n\n');
fprintf(fid, '五、本次运行结果摘要\n');
for i = 1:height(summaryTable)
    fprintf(fid, '%s：载荷 Fx=%.6g N, Fy=%.6g N, Fz=%.6g N；承载滚动体数=%.0f；最大接触载荷=%.6g N；最大接触应力=%.6g Pa；最小油膜厚度=%.6g m；打滑率=%.6g；等效刚度x=%.6g N/m；等效刚度y=%.6g N/m；寿命/损伤指标=%.6g；warning=%s；说明=%s。\n', ...
        summaryTable.bearing_type{i}, summaryTable.Fx_N(i), summaryTable.Fy_N(i), summaryTable.Fz_N(i), ...
        summaryTable.loaded_element_count(i), summaryTable.max_contact_load(i), summaryTable.max_contact_stress(i), ...
        summaryTable.min_oil_film_thickness(i), summaryTable.cage_slip_ratio(i), summaryTable.equivalent_stiffness_x(i), ...
        summaryTable.equivalent_stiffness_y(i), summaryTable.life_or_damage_index(i), summaryTable.warning_flag{i}, summaryTable.comment{i});
end
fprintf(fid, '\n六、注意事项\n');
fprintf(fid, '本次结果是当前程序、当前输入和当前统计口径下的基本工况输出。滚子轴承仍可能出现 singular_matrix 数值警告；若用于论文，应继续复核刚度矩阵条件数、油膜厚度异常值和接触载荷量级。\n');
fclose(fid);

disp('NEW_BASIC_LOAD_CASE_DONE');
disp(desktopOut);

function record = run_ball_case(projectRoot, datafromvb, micro_config, outDir)
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(fullfile(projectRoot, '球轴承程序'));
record = empty_record();
record.bearing_type = 'ball';
record.Fx_N = datafromvb(16);
record.Fy_N = datafromvb(17);
record.Fz_N = datafromvb(18);
lastwarn('');
try
    runLog = evalc('[loadj, returndata] = qiujieend(datafromvb, micro_config);');
    [record.warning_flag, record.comment] = classify_warning(runLog);
    record = fill_metrics(record, "ball", loadj, returndata);
    copyfile(fullfile(projectRoot, 'results_bearing_base', 'ball', '*'), outDir);
catch ME
    record.warning_flag = 'error';
    record.comment = sprintf('%s at %s:%d', ME.message, error_file(ME), error_line(ME));
end
end

function record = run_roller_case(projectRoot, datafromvb, micro_config, outDir)
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(fullfile(projectRoot, '滚子轴承程序'));
record = empty_record();
record.bearing_type = 'roller';
record.Fx_N = NaN;
record.Fy_N = datafromvb(19);
record.Fz_N = datafromvb(20);
lastwarn('');
try
    runLog = evalc('[loadj, returndata] = qiujieall(datafromvb, micro_config);');
    [record.warning_flag, record.comment] = classify_warning(runLog);
    record = fill_metrics(record, "roller", loadj, returndata);
    copyfile(fullfile(projectRoot, 'results_bearing_base', 'roller', '*'), outDir);
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

function record = fill_metrics(record, bearingType, loadj, returndata)
record.loaded_element_count = loadj;
record.max_contact_load = max_from_mats({'Q1.mat', 'Q2.mat'});
record.max_contact_stress = max_from_mats({'Ph1.mat', 'Ph2.mat'});
record.min_oil_film_thickness = min_from_mats({'oilh1.mat', 'oilh2.mat'});
record.PV_value = max_from_mats({'pvzhi1.mat', 'pvzhi2.mat', 'pvzhinonload.mat', 'pvzhi1nonload.mat'});
if bearingType == "ball"
    record.cage_slip_ratio = value_or_nan(returndata, 6);
    record.equivalent_stiffness_x = value_or_nan(returndata, 1);
    record.equivalent_stiffness_y = value_or_nan(returndata, 2);
    record.life_or_damage_index = value_or_nan(returndata, 9);
else
    record.cage_slip_ratio = value_or_nan(returndata, 5);
    record.equivalent_stiffness_x = NaN;
    record.equivalent_stiffness_y = value_or_nan(returndata, 1);
    record.life_or_damage_index = value_or_nan(returndata, 6);
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
    'Fx_N', NaN, ...
    'Fy_N', NaN, ...
    'Fz_N', NaN, ...
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

function datafromvb = default_ball_input()
datafromvb = zeros(75, 1);
datafromvb(1)=15; datafromvb(2)=0.02223; datafromvb(3)=0.1253; datafromvb(4)=0.5232; datafromvb(5)=0.5232; datafromvb(6)=40; datafromvb(7)=0.0115; datafromvb(35)=5;
datafromvb(61)=7800; datafromvb(68)=7800; datafromvb(51)=7800; datafromvb(50)=7800; datafromvb(9)=2.06e+11; datafromvb(10)=2.06e+11; datafromvb(11)=2.06e+11;
datafromvb(65)=2.06e+11; datafromvb(12)=0.3; datafromvb(13)=0.3; datafromvb(14)=0.3; datafromvb(66)=0.3; datafromvb(69)=0; datafromvb(15)=10000; datafromvb(16)=20000; datafromvb(17)=0;
datafromvb(18)=3000; datafromvb(19)=0; datafromvb(20)=0; datafromvb(21)=0.275; datafromvb(22)=0.225; datafromvb(23)=1e-07; datafromvb(24)=1e-07; datafromvb(25)=1e-07; datafromvb(26)=970;
datafromvb(27)=20; datafromvb(28)=0.0966; datafromvb(29)=0.0318; datafromvb(30)=1.28e-008; datafromvb(31)=3.2e-002; datafromvb(32)=1; datafromvb(33)=0.4e-03; datafromvb(34)=3.8;
datafromvb(74)=2; datafromvb(36)=370e-03; datafromvb(37)=220e-03; datafromvb(38)=120.0e-03; datafromvb(39)=225e-03; datafromvb(40)=1.96e11; datafromvb(41)=2.18e11; datafromvb(42)=0.3; datafromvb(43)=0.3;
datafromvb(44)=0.5e-03; datafromvb(45)=-0.5e-03; datafromvb(54)=11.6e-06; datafromvb(55)=11.8e-06; datafromvb(56)=11.8e-06; datafromvb(57)=11.8e-06; datafromvb(58)=11.8e-06; datafromvb(46)=180;
datafromvb(47)=170; datafromvb(48)=190; datafromvb(49)=27; datafromvb(59)=195; datafromvb(60)=160; datafromvb(50)=7870; datafromvb(51)=7870; datafromvb(52)=7860; datafromvb(53)=8360;
datafromvb(70)=0; datafromvb(71)=0.2; datafromvb(72)=0.2; datafromvb(73)=0.2;
end
