function displayParameters(params)
%% DISPLAYPARAMETERS
% The function displays the alphabetically organized parameters in params
% if the params.verbose variable is larger than 0.
%
% INPUT
%   params: Struct containing subfields with the parameter name and value,
%       for example params.smoothing = 2. 
%
% OUTPUT
%   prints to workspace values for all parameters.
%
% Written by Andreas S Lande 2019

% Display output if verbose is positive
global verbose
if verbose 
    fprintf('-----------------\nParameters used:\n');

    % Find all parameters and sort them alphabetically
    [~, neworder] = sort(lower(fieldnames(params)));
    params = orderfields(params, neworder);
    param_fields = fields(params);

    % Print name and value of all parameters
    for x = 1:length(param_fields)

        % For all fields except "display_output" print fieldname and value
        if ~strcmp(param_fields{x},"verbose")
            if isstring(params.(param_fields{x}))
                fprintf('+ %s: %s\n',param_fields{x},params.(param_fields{x}));
            else
                fprintf('+ %s: %s\n',param_fields{x},num2str(params.(param_fields{x})));
            end
        end

    end
    
    fprintf('-----------------\n');
    
end

end