
function [outputPlaceholder] = choiceAnalyzer_DS(inputPlaceholder) 
% input should maybe be helper functions created, input should be whatever otiginal .py returns

% Created by DS on 06/16/26
% Python adaptation of Stephanie Prince's choice_analyzer python code 

% need to add python package and module equivalents up here
    % Block one: built in python libraries
    % Block two: importation of built in pyton modules AND final line
    % imports from a built in file
    % Block Three: imports from built in files from elsewhere 

% 0a.Creating BaseAnalysisClass: loads, saves and exports data. Can be
% implemented elsewhere 

% 1.creating the ChoiceAnalyzer class
    classdef ChoiceAnalyzer < handle 
        properties
            nwbfile
            session_id
            target_var = 'choice'; 
            velocity_only = false;


            % setup paramas
            mask_value -9999;
            params = struct(...
               'batch_size', 32, ...
               'epochs', 20, ...
               'regularizer', [], ...
               'learning_rate', 0.1, ...
               'predict_update' = true ...
               );            
            grid_search_params = struct('batch_size', [20, 50, 100], ...
                                       'epochs', [10, 20, 30], ...
                                       'regularizer', {{[], 'l2(0.01)', 'l2(0.1)'}}, ...
                                       'learning_rate', [0.01, 0.1], ...
                                       'predict_update', false ...
                                       ); % has placeholders for tensorflow


            % setup data
            trials_df = table()
            session_ts = []
            
            input_data = {}
            target_data = {}
            timestamp_data = {}

            non_update_index = []

            update_input_data = {}
            update_target_data = {}
            update_timestamp_data = {}
            update_index = []

            max_pad_length = []

            % Get results to save
            results_io = []
            data_files = struct()

        end
            
        methods
            function obj = ChoiceAnalyzer(nwbfile, session_id, target_var, velocity_only)
                if nargin >= 1 && ~empty(nwbfile)
                    obj.nwbfile = nwbfile;
                end
                if nargin >= 2  && ~empty(session_id)
                    obj.session_id = session_id;
                end
                if nargin >= 3  && ~empty(target_var)
                    obj.target_var = target_var;
                end
                if nargin >= 4  && ~empty(velocity_only)
                    obj.velocity_only = velocity_only;
                end
                
               % setup data methods
               obj.trials_df = util.table.fromNWB(nwbfile.intervals_trials);
               obj.session_ts = nwbfile.processing.get('behavior') ...
                .nwbdatainterface.get('view_angle').view_angle;
            
                [input_data, obj.target_data, obj.timestamp_data] = ...
                    obj.setup_data(nwbfile, obj.target_var);
                
                [~, obj.non_update_index] = ...
                    obj.get_trial_inds('with_update', false, 'ret_index', true);


                
                if obj.params.predict_update
                    [input, target, timestamp] = ...
                        obj.setup_data(nwbfile, obj.target_var, true);

                    obj.update_input_data = input;
                    obj.update_target_data = target;
                    obj.update_timestamp_data = timestamp;
                
                    [~, obj.update_index] = get_trial_inds(true, true);
                
                    obj.max_pad_length = max( ...
                        [ max(cellfun(@length, obj.input_data)), ...
                          max(cellfun(@length, obj.update_input_data)) ] );
                
                else
                    obj.update_input_data = {};  %might beed to change later depending on what data type is being stored 
                    obj.update_target_data = {};
                    obj.update_index = {};
                    obj.max_pad_length = max(cellfun(@length, input_data)); % returns length of longest cell in array
                end
                % Get results to save
                if obj.velocity_only
                    add_tags = 'velocity_only';
                else 
                   add_tags = '';
                end    
                
               creator_file = [mfilename('fullpath'), '.m']
               class_file_folder = fileparts(mfilename('fullpath'));
               [~, folder_name] = fileparts(class_file_folder)

               obj.results_io = resultsIO( ...
                  creator_file, ...
                  obj.session_id, ...
                  folder_name, ...
                  [obj.target_var, add_tags]...
                );

               obj.data_files = struct();


               obj.data_files.dynamic_choice_output = struct( ...
                'vars', { {'output_data', 'agg_data', 'decoder_data', 'params'} }, ...
                'format', 'pkl' ...
                );
    
            end
    
            % 2. run_analysis function 
            function run_analysis(obj, overwrite, grid_search)
                if nargin < 2
                    overwrite = false;
                end
                if nargin < 3
                   grid_search = false; 
                end
    
                if grid_search
                    obj.data_files = struct( ...
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
                function grid_search(obj)
                    paramLists = struct2cell(obj.grid_search_params);
                    T = combinations(paramLists{:});
                    obj.grid_search_data = T;
                    for i = 1:height(T)
                        batch_size    = T{i, 1};
                        epochs        = T{i, 2};
                        regularizer   = T{i, 3};
                        learning_rate = T{i, 4};
                    end
                
                    cv = cvpartition(target_data(:,1), 'KFold', 5);
                    for i = 1:cv.NumTestSets
                        train_index = training(cv, i);   
                        test_index  = test(cv, i);       
                    end
                   
                    
                end  
    
            % 4. get_dynamic_choice function 
            function get_dynamic_choice(obj)

            end
    
            % 5. _get_classifier function 
            function get_classifier(obj)

            end
    
            % 6. setup_data function
            function setup_data(obj)

            end
    
            % 7._pad_data function
            function pad_data(obj)

            end
    
            % 8. _preprocess_data function
            function preprocess_data(obj)

            end
    
            % 9. get_trial_inds function
            function get_trial_inds(obj)

            end
    
            % 10. _aggregate_data function
            function aggregate_data(obj)

            end
    
            % 11. _get_decoder_data function
            function get_decoder_data(obj)

            end
    
            % 12._log2_likelihood function
            function log2_likelihood(obj)

            end
    
            % 13. _get_repeated_fold_average function
            function get_repeated_fold_average(obj)
                
            end
        end    
    end
        
        
       