
        %% From getBehaviorROC, behaves like Python align_by_time_intervals module
        function [data, groupData] = extractZoneBinsByLapGroup_DS(data, groupData, lapData, statsByLap, params, currrEnv, bp_lapGroups, sessionInfoRow) 
            
            % makes sure current enviornment exists in the data
            if ~isfield(data, currEnv)
                data.(currEnv).az = [];
                data.(currEnv).cz = [];
                data.(currEnv).nevrz = [];
                data.(currEnv).sessionInfo = [];
            end


            % checks if each field exists in currEnv
            if ~isfield(data.(currEnv), 'az')
                data.(currEnv).az = [];
            end

            if ~isfield(data.(currEnv), 'cz')
                data.(currEnv).cz = [];
            end

            if ~isfield(data.(currEnv), 'nevrz')
                data.(currEnv).nevrz = [];
            end

            if ~isfield(data.(currEnv), 'sessionInfo')
                data.(currEnv).sessionInfo = [];
            end

            
            if ~isfield(data, groupData)
                groupData.(currEnv) = struct('az', {}, 'cz', {}, 'nevrz', {});
            end

            numGroups = max(bp_lapGroups)

             for g = 1:numGroups

                 if length(groupData.(currEnv)) < g
                     groupData.(currEnv)(g).az = [];
                     groupData.(currEnv)(g).cz = [];
                     groupData.(currEnv)(g).nevrz = [];
                 end
            
                 if ~isfield(groupData.(currEnv)(g), 'az')
                     groupData.(currEnv)(g).az = [];
                 end

                 if ~isfield(groupData.(currEnv)(g), 'cz')
                     dagroupDatata.(currEnv)(g).cz = [];
                 end

                 if ~isfield(groupData.(currEnv)(g), 'nevrz')
                     groupData.(currEnv)(g).nevrz = [];
                 end


         

                 for lp = 1:size(lapData,1)
                     currGroup = bp_lapGroups(lp); %loops through each lap number using lookup vector and stores it in new vector
                     if currGroup ==0
                         continue
                     end


                % Deandra's code

                for zn = 1:length(params.Azones)
                    tmpBins = [];
                    tmpBins = params.Azones(zn):params.binsize_deg:params.Azones(zn)+statsByLap.fileInfo.cueSize-params.binsize_deg;
                    tmpBins = round(tmpBins/params.binsize_deg);%DC adding round to handle offset RZ with new projector
                    data.(currEnv).az = [data.(currEnv).az; lapData(lp,tmpBins)];

                    % Deandra:  fills az box with az values
                    groupData.(currEnv)(currGroup).az = [groupData.(currEnv)(currGroup).az; lapData(lp,tmpBins)];

                end
                for zn = 1:length(params.NRzones)
                    tmpBins = [];
                    tmpBins = params.NRzones(zn):params.binsize_deg:params.NRzones(zn)+statsByLap.fileInfo.cueSize-params.binsize_deg;
                    tmpBins = round(tmpBins/params.binsize_deg);
                    data.(currEnv).cz = [data.(currEnv).cz; lapData(lp,tmpBins)];
                    % Deandra:  fills cz box with cz values
                    groupData.(currEnv)(currGroup).cz = [groupData.(currEnv)(currGroup).cz; lapData(lp,tmpBins)];

                end
                for zn = 1:length(params.NevRzones)
                    if ~isnan(params.NevRzones(zn))
                        tmpBins = [];
                        tmpBins = params.NevRzones(zn):params.binsize_deg:params.NevRzones(zn)+statsByLap.fileInfo.cueSize-params.binsize_deg;
                        tmpBins = round(tmpBins/params.binsize_deg);
                        data.(currEnv).nevrz = [data.(currEnv).nevrz; lapData(lp,tmpBins)];
                        % Deandra: fills nevrz box with nevrz values
                        groupData.(currEnv)(currGroup).nevrz = [groupData.(currEnv)(currGroup).nevrz; lapData(lp,tmpBins)];

                    else
                        data.(currEnv).nevrz = [];
                    end
                end

                     % Save zone/bin metadata into data.
                     data.(currEnv).azBins_deg = params.Azones;
                     data.(currEnv).czBins_deg = params.NRzones;
                     data.(currEnv).nevrzBins_deg = params.NevRzones;
    
                     % Save session info into data if it was provided.
                     if nargin >= 8 && ~isempty(sessionInfoRow)a
                         data.(currEnv).sessionInfo = [data.(currEnv).sessionInfo; sessionInfoRow];
    
                     end
                 
                 end
             end
        end


