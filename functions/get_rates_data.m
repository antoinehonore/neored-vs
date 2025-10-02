function rates_data = get_rates_data(evt_rates_dose)
    rates_data = eval(strrep(evt_rates_dose,' ',' '));
    if size(rates_data,2) >1
        %fprintf('RATE %s: %s\n',vs.socsecurity, strrep(vs.evt_rates_dose,' ',' '))
        rates_data=rates_data(1);
    end
    rates_data = mean(rates_data);
end

