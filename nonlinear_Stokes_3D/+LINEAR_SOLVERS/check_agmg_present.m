function is_agmg_present = check_agmg_present(agmg_folder)
%CHECK_AGMG_PRESENT Adds AGMG to path and checks wrapper plus MEX backend.

if nargin < 1 || isempty(agmg_folder) || ~isfolder(agmg_folder)
    is_agmg_present = 0;
    return;
end

addpath(agmg_folder);

has_wrapper = exist(fullfile(agmg_folder, 'agmg.m'), 'file') == 2 && exist('agmg', 'file') == 2;
has_backend = exist('dmtlagmg', 'file') == 3 || exist('zmtlagmg', 'file') == 3 || ...
    exist('dmtlagmg_sm', 'file') == 3 || exist('zmtlagmg_sm', 'file') == 3;

is_agmg_present = has_wrapper && has_backend;
end
