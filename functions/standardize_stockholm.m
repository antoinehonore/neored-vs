function [vs,studyid,name] = standardize_stockholm(vs,studyid,files,f)
        % Here, we may have to add some extra lines of code to make the
        % Berlin and Stockholm data similar to the Oxford data. 
        vs.birthweight = 1000*vs.birthweight;


        if isfield(vs,'evt_start_days')
            vs.evt_start = vs.evt_start_days*24;
            vs.evt_end = vs.evt_end_days*24;
        end
        vs.time_events_vital_signs = [vs.evt_start; vs.evt_end];
        % if ~isfield(vs,'time')
        vs.time = vs.pna_days* 24;
        % end
        vs.events_vital_signs = {'rbctransfusion_start','rbctransfusion_stop'};
        vs.HR = vs.btb;
        vs.pulse = vs.btb;   %%%%   /!\ not the pulse derived from pulse oxymetry  
        vs.sats = vs.spo2;
        vs.RR = vs.rf;
        if vs.sex == 'U'
            fprintf('%s,%s\n',vs.socsecurity,vs.sex )
            vs.sex = NaN;
        end
        if isnan(vs.sex)
            fprintf('SEX,%s\n',vs.socsecurity)
            vs.sex='-';
        end

        rates_data = eval(strrep(vs.evt_rates_dose,' ',' '));
        if size(rates_data,2) >1
            %fprintf('RATE %s: %s\n',vs.socsecurity, strrep(vs.evt_rates_dose,' ',' '))
            rates_data=rates_data(1);
        end
        
        vs.rate = mean(rates_data);

        if isfloat(vs.respirator_pre)
            if isnan(vs.respirator_pre)
                vs.respirator_pre='StandBy';
            end
        end

        if isfloat(vs.respirator_post)
            if isnan(vs.respirator_post)
                vs.respirator_post='StandBy';
            end
        end
        
        vs.respirator = unique({vs.respirator_pre, vs.respirator_post});

        % Correct potential problems in data
        if vs.weight_kg_pre>5
            vs.weight_kg_pre=NaN;
        end
        if isnan(vs.weight_kg_pre)
            fprintf('WEIGHT,%s,%f\n',vs.socsecurity, vs.evt_start_days)
        end
        if isnan(vs.hb_pre)
            %fprintf('HB %s\n',vs.socsecurity);
        end

        if isfield(vs,'fio2_pre')
            vs.fiO2_pre = vs.fio2_pre;
            vs.fiO2_post = vs.fio2_post;
        end
        
        vs = rmfield(vs, {'btb','spo2','rf',});
        
        % the filename contains the patid: /home/.../[patid]/RBC/[number].mat
        [filepath, name, ext] = fileparts(files{f,1});
        [filepath, name, ext] = fileparts(filepath);
        [filepath, name, ext] = fileparts(filepath);
        studyid = [studyid;f];

end