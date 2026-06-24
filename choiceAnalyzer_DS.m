
function [outputPlaceholder] = choiceAnalyzer_DS(inputPlaceholder) 
% input should maybe be helper functions created, input should be whatever otiginal .py returns

% Created by DS on 06/16/26
% Python adaptation of Stephanie Prince's choice_analyzer python code 

% need to add python package and module equivalents up here
    % Block one: built in python libraries
    % Block two: importation of built in pyton modules AND final line
    % imports from a built in file
    % Block Three: imports from built in files from elsewhere 

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
        % several params throughout the entire funtion needs to be inserted 

        % setup paramas
        velocity_only = False;
        mask_value = -9999
        params = dictionary([batch_size, epochs, regularizer, learning_rate, predict_update] , [32, 20, None, 0.1, True])
        grid_search_params = dictionary([batch_size, epochs, regularizer, learning_rate, predict_update], [{20 50 100, 10 20 30}, (regularizer = {[], 0.01, 0.1}),
                                                    ("L2Regularization", 0.1), {0.01 0.1}, False])
        
        % setup data
        trials_df = util.table.fromNWB(nwbfile.intervals_trials); % this uses the MATLAB version of NWB
        session_ts = nwbfile.processing.get('behavior') ...
        .nwbdatainterface.get('view_angle').view_angle;
    
        [input_data, target_data, timestamp_data] = ...
            setup_data(nwbfile, target_var);
        
        [~, non_update_index] = ...
            get_trial_inds('with_update', false, 'ret_index', true);

        
        if isfield(params, 'predict_update') && params.predict_update

            [input, target, timestamp] = setup_data(nwbfile, target_var, true);
            update_input_data = input;
            update_target_data = target;
            update_timestamp_data = timestamp;
        
            [~, update_index] = get_trial_inds(true, true);
        
            max_pad_length = max( ...
                [ max(cellfun(@length, input_data)), ...
                  max(cellfun(@length, update_input_data)) ] );
        
        else
            update_input_data = {};  %might beed to change later depending on what data type is being stored 
            update_target_data = {};
            update_index = {};
            max_pad_length = max(cellfun(@length, input_data)); % returns length of longest cell in array
        end

        % Get results to save
        if velocity_only
            add_tags = 'velocity_only';
        else 
           add_tags = '';
        end    
        
       % still in Python
       results_io = ResultsIO(creator_file=__file__, session_id=session_id, folder_name=Path(__file__).parent.stem,
                                    tags=f'{target_var}{add_tags}') % uses separate file (ResultsIO, also calls another file within repo)
       % back to MATLAB
       data_files = struct( ...
        'dynamic_choice_output', struct( ...
        'vars', { {'output_data', 'agg_data', 'decoder_data', 'params'} }, ...
        'format', 'pkl' ...
            ) ...
        );


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
            % etc... for the rest of the inserted properties 
        end

        % 2. run_analysis function 
        function obj = run_analysis(overwrite, grid_search)
            if nargin < 1
                overwrite = false;
            end
            if nargin < 2
               grid_search = false; 
            end

            if grid_search
                data_files = struct( ...
                    "grid_search", struct( ...
                    'vars', {{'grid_search_data', 'grid_search_params'}}, ...
                    ' format', 'pkl' ...
                        )...
                    );
            end
            if overwrite
                if grid_search 
                    obj.grid_search();
                else
                    obj.get_dynamic_choice();
                    obj.aggregate_data();
                    obj.get_decoder_data();
                    obj.export_data();
                end
                obj.export_data();

            else 
                if obj.results_io.data_exists(obj.data_files)
                    obj.load_data();
                else
                    warning('Data with those input parameters does not exist, setting overwrite to true')
                    obj.run_analysis(true)
                end
            end
        end 
        

        % 3. _grid_search function
        function obj = grid_search()
            grid_search_data = {};
            for batch_size, epochs, regularizer, learning_rate in 1:length(itertools.product(*list(self.grid_search_params.values()))) % fix syntax
            for i = 1:length(itertools, {grid_search_params}, {values})
                batch_size =
                epochs =
                regularizer =
                learning_rate =
        end

        % 4. get_dynamic_choice function 
        function obj = get_dynamic_choice
        end

        % 5. _get_classifier function 
        function obj = get_classifier
        end

        % 6. setup_data function
        function obj = setup_data
        end

        % 7._pad_data function
        function obj = pad_data
        end

        % 8. _preprocess_data function
        function obj = preprocess_data
        end

        % 9. get_trial_inds function
        function obj = get_trial_inds
        end

        % 10. _aggregate_data function
        function obj = aggregate_data
        end

        % 11. _get_decoder_data function
        function obj = get_decoder_data
        end

        % 12._log2_likelihood function
        function obj = log2_likelihood
        end

        % 13. _get_repeated_fold_average function
        function  obj = get_repeated_fold_average
        end
    end

end

        
       