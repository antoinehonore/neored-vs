
function print_demographics(BIRTH_weight,BIRTH_Gender,PMA)
    n = size(BIRTH_weight,2);
    if n==0
        fprintf('No Data to print\n');
    else
        assert(size(BIRTH_Gender,2) == n);
        assert(size(PMA,2) == n);
        %B = regexp(BIRTH_Gender,'\S+','match');
        %T = cell2table([B{:}].');

        counts = [0,0,0]; % Always 3 categories

if all(cellfun(@isnumeric, BIRTH_Gender))  % Case: {0,1,0,1,...}
    % Convert cell array to numeric vector
    bg = cell2mat(BIRTH_Gender);
    counts(1) = sum(bg == 0);              % female
    counts(2) = sum(bg == 1);              % male
    counts(3) = sum(~ismember(bg,[0 1]));  % other / invalid

        elseif iscellstr(BIRTH_Gender) || isstring(BIRTH_Gender)
            % for datasets: '-', 'M', 'F'
            ic = 1;
            for c = {'-','M','F'}
                counts(ic) = sum(strcmp(BIRTH_Gender, c));
                ic = ic + 1;
            end
        else
            error('BIRTH_Gender has unexpected type');
        end
    
        %BIRTH_Gender_table = groupsummary(T,'Var1');
        %if size(BIRTH_Gender_table,1)==2
        %    BIRTH_Gender_table = [{{'-'},0};BIRTH_Gender_table];
        %end
        
        BIRTH_Gender_str = sprintf('Gender:\t\t\tF: %d (%.2f), M: %d (%.2f)\t(Missing: %d/%d)',counts(3),counts(3)/(n-counts(1)),counts(2),counts(2)/(n-counts(1)), counts(1),n);
        PMA_str = sprintf('Birth PMA (weeks): \t%.2f (%.2f-%.2f)\t(Missing: %d/%d)', median(PMA,'omitnan'), quantile(PMA,0.25),quantile(PMA,0.75),sum(isnan(PMA)),n);
        BIRTH_weight_str = sprintf('Birth weight (g):\t%d (%d-%d)\t(Missing: %d/%d)', round(median(BIRTH_weight,'omitnan')), round(quantile(BIRTH_weight,0.25)),round(quantile(BIRTH_weight,0.75)),sum(isnan(BIRTH_weight)),n);
        n_b750 = size(BIRTH_weight(BIRTH_weight<750),2);
        n_b1000 = size(BIRTH_weight(BIRTH_weight<1000),2);
        n_b1500 = size(BIRTH_weight(BIRTH_weight<1500),2);
        BIRTH_weight_str_counts = sprintf('%d (%d%%)\n%d (%d%%)\n%d (%d%%)',n_b750,round(100*n_b750/n), n_b1000, round(100*n_b1000/n), n_b1500, round(100*n_b1500/n));
        fprintf('%s\n%s\n%s\n%s\n', PMA_str, BIRTH_weight_str, BIRTH_weight_str_counts,BIRTH_Gender_str);
    end
end
