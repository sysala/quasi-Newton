function paths = mp2d_paths()
script_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(script_dir);

paths = struct();
paths.root = root_dir;
paths.scripts = script_dir;
paths.images = fullfile(root_dir, 'images');
paths.results = fullfile(root_dir, 'results');

if ~exist(paths.images, 'dir')
    mkdir(paths.images);
end
if ~exist(paths.results, 'dir')
    mkdir(paths.results);
end
end
