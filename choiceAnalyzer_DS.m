function [outputPlaceholder] = choiceAnalyzer_DS(inputPlaceholder) 
% input should maybe be helper functions created, input should be whatever otiginal .py returns

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

        % setup paramas
        velocity_only = False;
        mask_value = -9999
        params = dictionary([batch_size, epochs, regularizer, learning_rate, predict_update] , [32, 20, None, 0.1, True])
        grid_search_params = dictionary([batch_size, epochs, regularizer, learning_rate, predict_update], [{20 50 100, 10 20 30}, (regularizer = {[], 0.01, 0.1}),
                                                    ("L2Regularization", 0.1), {0.01 0.1}, False])
        % setup data
        trials_df = trialsTable = nwbfile.intervals_trials.toTable(); % check format of this and the below
        session_ts =

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
        % 2. run_analysis function 
        % 3. _grid_search function
        % 4. get_dynamic_choice function 
        % 5. _get_classifier function 
        % 6. _setup_data function
        % 7._pad_data function
        % 8. _preprocess_data function
        % 9. _get_trial_inds function
        % 10. _aggregate_data function
        % 11. _get_decoder_data function
        % 12._log2_likelihood function
        % 13. _get_repeated_fold_average function



    end

        
       