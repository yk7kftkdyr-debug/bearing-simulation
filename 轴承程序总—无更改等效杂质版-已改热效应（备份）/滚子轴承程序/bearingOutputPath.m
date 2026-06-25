function outpath = bearingOutputPath(varargin)
%BEARINGOUTPUTPATH Return and create the shared project output path.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
parts = varargin;

if numel(parts) >= 1 && strcmp(parts{1}, 'bearing')
    parts = parts(2:end);
end
if numel(parts) >= 1 && strcmp(parts{1}, 'qiu')
    parts{1} = 'ball';
elseif numel(parts) >= 1 && strcmp(parts{1}, 'gunzi')
    parts{1} = 'roller';
end
if isempty(parts) || ~any(strcmp(parts{1}, {'ball','roller'}))
    parts = [{'roller'}, parts];
end

outpath = fullfile(projectRoot, 'results_bearing_base', parts{:});
outdir = fileparts(outpath);
if ~isempty(outdir) && exist(outdir, 'dir') ~= 7
    mkdir(outdir);
end
end
