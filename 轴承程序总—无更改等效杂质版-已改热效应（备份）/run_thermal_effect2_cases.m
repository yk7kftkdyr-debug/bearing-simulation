clear; clc;

projectRoot = fileparts(mfilename('fullpath'));
desktopOut = fullfile('/Users/bailingguan/Desktop', '热效应2');
if exist(desktopOut, 'dir') ~= 7
    mkdir(desktopOut);
end

addpath(projectRoot);
addpath(fullfile(projectRoot, '球轴承程序'));
addpath(fullfile(projectRoot, '滚子轴承程序'));

% Table supplied by user:
% ambient, outer ring, inner ring, rolling element, shaft, housing.
thermalTable = [
    20,  40,  60,  55,  65,  35;
    50,  75, 105,  95, 110,  70;
    80, 115, 150, 138, 155, 105;
   100, 140, 180, 165, 185, 130];

rho_oil = 970; % kg/m^3, from uploaded thermal-effect document
eta20 = oil_dynamic_viscosity(20, rho_oil);
records = repmat(empty_record(), 0, 1);

for row = 1:size(thermalTable, 1)
    ambient_T = thermalTable(row, 1);
    outer_T = thermalTable(row, 2);
    inner_T = thermalTable(row, 3);
    rolling_T = thermalTable(row, 4);
    shaft_T = thermalTable(row, 5);
    housing_T = thermalTable(row, 6);

    nu_T = oil_kinematic_viscosity(ambient_T);
    eta_T = oil_dynamic_viscosity(ambient_T, rho_oil);
    alpha_p_T = pressure_viscosity_coefficient(ambient_T);
    viscosity_factor = eta_T / eta20;

    cfg = make_micro_interface_config();
    cfg.thermal.enabled = true;
    cfg.thermal.temperature_C = ambient_T;
    cfg.thermal.viscosity_factor = viscosity_factor;
    cfg.thermal.thermal_clearance_factor = 1;
    cfg.thermal.Ct1 = 1;
    cfg.thermal.Ct2 = 1;
    cfg.lubricant.enabled = true;
    cfg.lubricant.eta0 = eta_T;
    cfg.lubricant.alpha_p = alpha_p_T;
    cfg.lubricant.ehl_correction_factor = 1;

    tempDir = fullfile(desktopOut, sprintf('ambient_%03dC', ambient_T));
    ballOut = fullfile(tempDir, 'ball');
    rollerOut = fullfile(tempDir, 'roller');
    ensure_dir(ballOut);
    ensure_dir(rollerOut);

    ballData = default_ball_input();
    ballData(16) = 20000;
    ballData(17) = 0;
    ballData(18) = 20000;
    ballData(26) = rho_oil;
    ballData(27) = ambient_T;
    ballData(29) = eta_T;
    ballData(30) = alpha_p_T;
    % Ball program mapping: 46 rolling element, 47 outer, 48 inner,
    % 49 ambient, 59 shaft, 60 housing.
    ballData(46) = rolling_T;
    ballData(47) = outer_T;
    ballData(48) = inner_T;
    ballData(49) = ambient_T;
    ballData(59) = shaft_T;
    ballData(60) = housing_T;
    records(end + 1, 1) = run_ball_case(projectRoot, ballData, cfg, ballOut, ...
        ambient_T, outer_T, inner_T, rolling_T, shaft_T, housing_T, ...
        nu_T, eta_T, alpha_p_T, viscosity_factor); %#ok<SAGROW>

    rollerData = defaultBearingInput();
    rollerData(19) = 0;
    rollerData(20) = 80000;
    rollerData(26) = rho_oil;
    rollerData(27) = ambient_T;
    rollerData(29) = eta_T;
    rollerData(30) = alpha_p_T;
    % Roller program mapping: 44 roller, 45 outer, 46 inner,
    % 47 ambient, 57 shaft, 58 housing.
    rollerData(44) = rolling_T;
    rollerData(45) = outer_T;
    rollerData(46) = inner_T;
    rollerData(47) = ambient_T;
    rollerData(57) = shaft_T;
    rollerData(58) = housing_T;
    records(end + 1, 1) = run_roller_case(projectRoot, rollerData, cfg, rollerOut, ...
        ambient_T, outer_T, inner_T, rolling_T, shaft_T, housing_T, ...
        nu_T, eta_T, alpha_p_T, viscosity_factor); %#ok<SAGROW>

    writetable(struct2table(records), fullfile(desktopOut, 'thermal_effect2_summary.csv'), 'Encoding', 'UTF-8');
end

writetable(struct2table(records), fullfile(desktopOut, 'thermal_effect2_summary.csv'), 'Encoding', 'UTF-8');
disp('THERMAL_EFFECT2_CASES_DONE');
disp(desktopOut);

function nu = oil_kinematic_viscosity(T)
nu = 27.6 * exp(-0.02815 * (T - 40)); % cSt
end

function eta = oil_dynamic_viscosity(T, rho)
eta = rho * oil_kinematic_viscosity(T) * 1e-6; % Pa*s
end

function alpha_p = pressure_viscosity_coefficient(T)
alpha_p = 1.28e-8 * exp(-0.010 * (T - 20)); % Pa^-1
end

function ensure_dir(pathName)
if exist(pathName, 'dir') ~= 7
    mkdir(pathName);
end
end

function record = run_ball_case(projectRoot, datafromvb, cfg, outDir, ambient_T, outer_T, inner_T, rolling_T, shaft_T, housing_T, nu_T, eta_T, alpha_p_T, viscosity_factor)
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(fullfile(projectRoot, '球轴承程序'));
record = base_record('ball', ambient_T, outer_T, inner_T, rolling_T, shaft_T, housing_T, nu_T, eta_T, alpha_p_T, viscosity_factor);
record.Fx_N = datafromvb(16);
record.Fy_N = datafromvb(17);
record.Fz_N = datafromvb(18);
lastwarn('');
try
    runLog = evalc('[loadj, returndata] = qiujieend(datafromvb, cfg);');
    [record.warning_flag, record.comment] = classify_warning(runLog);
    record = fill_metrics(record, "ball", loadj, returndata);
    copyfile(fullfile(projectRoot, 'results_bearing_base', 'ball', '*'), outDir);
catch ME
    record.warning_flag = 'error';
    record.comment = sprintf('%s at %s:%d', ME.message, error_file(ME), error_line(ME));
end
end

function record = run_roller_case(projectRoot, datafromvb, cfg, outDir, ambient_T, outer_T, inner_T, rolling_T, shaft_T, housing_T, nu_T, eta_T, alpha_p_T, viscosity_factor)
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(fullfile(projectRoot, '滚子轴承程序'));
record = base_record('roller', ambient_T, outer_T, inner_T, rolling_T, shaft_T, housing_T, nu_T, eta_T, alpha_p_T, viscosity_factor);
record.Fx_N = NaN;
record.Fy_N = datafromvb(19);
record.Fz_N = datafromvb(20);
lastwarn('');
try
    runLog = evalc('[loadj, returndata] = qiujieall(datafromvb, cfg);');
    [record.warning_flag, record.comment] = classify_warning(runLog);
    record = fill_metrics(record, "roller", loadj, returndata);
    copyfile(fullfile(projectRoot, 'results_bearing_base', 'roller', '*'), outDir);
catch ME
    record.warning_flag = 'error';
    record.comment = sprintf('%s at %s:%d', ME.message, error_file(ME), error_line(ME));
end
end

function record = base_record(bearingType, ambient_T, outer_T, inner_T, rolling_T, shaft_T, housing_T, nu_T, eta_T, alpha_p_T, viscosity_factor)
record = empty_record();
record.bearing_type = bearingType;
record.ambient_temperature_C = ambient_T;
record.outer_ring_temperature_C = outer_T;
record.inner_ring_temperature_C = inner_T;
record.rolling_element_temperature_C = rolling_T;
record.shaft_temperature_C = shaft_T;
record.housing_temperature_C = housing_T;
record.kinematic_viscosity_cSt = nu_T;
record.dynamic_viscosity_Pa_s = eta_T;
record.pressure_viscosity_Pa_inv = alpha_p_T;
record.viscosity_factor_vs_20C = viscosity_factor;
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
    'ambient_temperature_C', NaN, ...
    'outer_ring_temperature_C', NaN, ...
    'inner_ring_temperature_C', NaN, ...
    'rolling_element_temperature_C', NaN, ...
    'shaft_temperature_C', NaN, ...
    'housing_temperature_C', NaN, ...
    'kinematic_viscosity_cSt', NaN, ...
    'dynamic_viscosity_Pa_s', NaN, ...
    'pressure_viscosity_Pa_inv', NaN, ...
    'viscosity_factor_vs_20C', NaN, ...
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
