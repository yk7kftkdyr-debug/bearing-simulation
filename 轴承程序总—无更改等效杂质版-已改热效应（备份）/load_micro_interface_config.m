function micro_config = load_micro_interface_config()
%LOAD_MICRO_INTERFACE_CONFIG Load runtime micro-interface configuration.

if exist('micro_config_runtime.mat', 'file') == 2
    data = load('micro_config_runtime.mat', 'micro_config');
    if isfield(data, 'micro_config')
        micro_config = data.micro_config;
        return;
    end
end

micro_config = make_micro_interface_config();
end
