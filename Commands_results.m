
resultsFolder = fullfile(pwd(), '/results');

% Find ONLY the final STRIKE-GOLDD result files
files = dir(fullfile(resultsFolder, 'id_results_*.mat'));

% Preallocate
nFiles = numel(files);

Model = strings(nFiles,1);
N_states = zeros(nFiles,1);
N_parameters = zeros(nFiles,1);
N_unknowns = zeros(nFiles,1);
Rank = zeros(nFiles,1);
N_observable_states = zeros(nFiles,1);
N_identifiable_parameters = zeros(nFiles,1);
N_unobservable_states = zeros(nFiles,1);
N_unidentifiable_parameters = zeros(nFiles,1);
FISPO = strings(nFiles,1);
Lie_derivatives = zeros(nFiles,1);

Identifiable_parameters = strings(nFiles,1);
Unidentifiable_parameters = strings(nFiles,1);
Observable_states = strings(nFiles,1);
Unobservable_states = strings(nFiles,1);

%% Read each result file

for i = 1:nFiles

    filePath = fullfile(files(i).folder, files(i).name);

    % Load only the variables we need
    S = load(filePath, ...
        'rango', ...
        'p', 'p_id', 'p_un', ...
        'x', 'obs_states', 'unobs_states', ...
        'isFISPO', 'nd');

    % Model name
    name = erase(files(i).name, 'id_results_');
    
    % Remove date at the end, e.g. "_19-Aug-2026"
    name = regexprep(name, '_\d{2}-[A-Za-z]{3}-\d{4}\.mat$', '');
    
    Model(i) = string(name);

    % Number of states and parameters
    N_states(i) = numel(S.x);
    N_parameters(i) = numel(S.p);

    % Total quantities tested by the OI matrix
    N_unknowns(i) = N_states(i) + N_parameters(i);

    % Rank
    Rank(i) = double(S.rango);

    % Identifiability
    N_identifiable_parameters(i) = numel(S.p_id);
    N_unidentifiable_parameters(i) = numel(S.p_un);

    % Observability
    N_observable_states(i) = numel(S.obs_states);
    N_unobservable_states(i) = numel(S.unobs_states);

    % FISPO
    if isempty(S.unobs_states) && isempty(S.p_un)
        FISPO(i) = "Yes";
    else
        FISPO(i) = "No";
    end

    % Number of Lie derivatives
    Lie_derivatives(i) = S.nd;

    % Store variable names as text
    if isempty(S.p_id)
        Identifiable_parameters(i) = "";
    else
        Identifiable_parameters(i) = strjoin( ...
            string(cellstr(char(S.p_id(:)))), ', ');
    end

    if isempty(S.p_un)
        Unidentifiable_parameters(i) = "";
    else
        Unidentifiable_parameters(i) = strjoin( ...
            string(cellstr(char(S.p_un(:)))), ', ');
    end

    if isempty(S.obs_states)
        Observable_states(i) = "";
    else
        Observable_states(i) = strjoin( ...
            string(cellstr(char(S.obs_states(:)))), ', ');
    end

    if isempty(S.unobs_states)
        Unobservable_states(i) = "";
    else
        Unobservable_states(i) = strjoin( ...
            string(cellstr(char(S.unobs_states(:)))), ', ');
    end

end

%% Create summary table

Summary = table( ...
    Model, ...
    N_states, ...
    N_parameters, ...
    N_unknowns, ...
    Rank, ...
    N_observable_states, ...
    N_identifiable_parameters, ...
    N_unobservable_states, ...
    N_unidentifiable_parameters, ...
    FISPO, ...
    Lie_derivatives, ...
    Identifiable_parameters, ...
    Unidentifiable_parameters, ...
    Observable_states, ...
    Unobservable_states);

%% Display

disp(Summary)

%% Export to Excel

writetable(Summary, ...
    fullfile(pwd(), 'SI_result_summary.xlsx'));

%%

S = load('id_results_CombinedHSC_19-Aug-2026.mat');

whos('-file', 'id_results_CombinedHSC_19-Aug-2026.mat')

S.isFISPO
class(S.isFISPO)

% for i = 1:numel(files)
%     fprintf('\n\n%s\n', files(i).name);
%     disp(whos('-file', fullfile(files(i).folder, files(i).name)));
% end
% 
% Summary = table( ...
%     strings(0,1), ...
%     strings(0,1), ...
%     zeros(0,1), ...
%     'VariableNames', {'File','Variable','Size'});
% 
% for i = 1:numel(files)
% 
%     filePath = fullfile(files(i).folder, files(i).name);
%     info = whos('-file', filePath);
% 
%     for j = 1:numel(info)
% 
%         newRow = table( ...
%             string(files(i).name), ...
%             string(info(j).name), ...
%             string(mat2str(info(j).size)), ...
%             'VariableNames', {'File','Variable','Size'});
% 
%         Summary = [Summary; newRow];
% 
%     end
% end
% 
% disp(Summary)

% 
% for i = 1:numel(files)
% 
%     filePath = fullfile(files(i).folder, files(i).name);
% 
%     fprintf('\n============================================\n');
%     fprintf('%s\n', files(i).name);
%     fprintf('============================================\n');
% 
%     info = whos('-file', filePath);
% 
%     for j = 1:numel(info)
%         fprintf('  %-35s  %s\n', ...
%             info(j).name, ...
%             mat2str(info(j).size));
%     end
% end