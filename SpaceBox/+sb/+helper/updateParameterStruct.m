function params = updateParameterStruct(params,user_input)
%% UPDATE PARAMETER STRUCT
% This is a helper function that takes a params struct as input and a cell
% array of alternating keyword and value given as user input to some
% functions. This function then updates necessary values in the params
% struct.
%
% INPUT
%   params: Struct containing key:value pairs. 
%   user_input: Cell array with alternating ("key1",value,"key2", value ...)
%
% OUTPUT
%   params: Struct containing key:value pairs.
%
% Written by Andreas S Lande 2019

% Grab global verbose variable
try
    global verbose
    
    if isempty(verbose)
        params.verbose = 1;
    else
        
        params.verbose = verbose;
    end
    
catch
   params.verbose = 1; 
end

% Update fields in params struct that are present in user_input
if ~isempty(user_input)
    update_fields = user_input(1:2:end);
    update_values = user_input(2:2:end);

    for f = 1:length(update_fields)
        field_name = update_fields{f};
        
        % Update parameter value
        params.(field_name) = update_values{f};  
    end
end


end

