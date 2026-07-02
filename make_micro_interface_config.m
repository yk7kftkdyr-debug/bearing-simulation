function micro_config = make_micro_interface_config(varargin)
%MAKE_MICRO_INTERFACE_CONFIG Unified micro-interface parameter defaults.
% Defaults describe the baseline undisturbed case used by the base examples.

micro_config = struct();

micro_config.thermal = struct();
micro_config.thermal.enabled = false;
micro_config.thermal.Ct1 = 1;
micro_config.thermal.Ct2 = 1;

micro_config.texture = struct();
micro_config.texture.enabled = false;
micro_config.texture.texture_type = 'none';
micro_config.texture.texture_depth = 0;
micro_config.texture.texture_width = 0;
micro_config.texture.texture_density = 0;
micro_config.texture.Cr = 1;
micro_config.texture.wenli = 6;

micro_config.debris = struct();
micro_config.debris.enabled = false;
micro_config.debris.debris_displacement = 0;
micro_config.debris.ud = 0;

if nargin == 0
    return;
end

if nargin == 1 && isstruct(varargin{1})
    micro_config = merge_micro_config(micro_config, varargin{1});
else
    if mod(nargin, 2) ~= 0
        error('make_micro_interface_config:InvalidInput', ...
            'Overrides must be provided as name-value pairs or one struct.');
    end
    overrides = struct();
    for k = 1:2:nargin
        overrides.(varargin{k}) = varargin{k + 1};
    end
    micro_config = merge_micro_config(micro_config, overrides);
end
end

function base = merge_micro_config(base, overrides)
names = fieldnames(overrides);
for i = 1:numel(names)
    name = names{i};
    if isfield(base, name) && isstruct(base.(name)) && isstruct(overrides.(name))
        base.(name) = merge_micro_config(base.(name), overrides.(name));
    else
        base.(name) = overrides.(name);
    end
end
end
