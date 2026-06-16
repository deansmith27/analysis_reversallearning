function [output 1, output2] = choiceAnalyzer_DS(input1, input2)
% Created by DS on 06/16/26
% Python adaptation of Stephanie Prince's choice_analyzer python code 

% need to add python package and module equivalents up here
    % emphasis on the modules update_project.general.results_io,
    % update_project.general.timeseries, and
    % update_project.base_analysis_class (these are file paths)

% 0a.Creating BaseAnalysisClass that is the input into ChoiceAnalyzer
classdef BaseAnalysisClass
    properties
    end
end



% 1.creating the ChoiceAnalyzer class
classdef ChoiceAnalyzer(BaseAnalysisClass) % fix formatting
    properties
        nwbfile
        session_id
        target_var = 'choice'; 
        velocity_only = False;
    end
        
    methods
        function obj = ChoiceAnalyzer(nwbfile, session_id, target_var, velocity_only)
            if nargin >= 1
                obj.nwbfile = nwbfile;
            end
            if nargin >= 2
                obj.session_id = session_id;
            end
            if nargin >= 3
                obj.target_var = target_var;
            end
            if nargin >= 4
                obj.velocity_only = velocity_only;
            end
        end
    end

        
        % velocity_only
        % mask_value 
        % params(% dict then used: MATLAB equivalent needed)
        % grid_search_params(% dict then used: MATLAB equivalent needed)