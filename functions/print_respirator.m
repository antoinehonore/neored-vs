
function print_respirator(respirators)
    all_respirators = unique(respirators);
    n = size(all_respirators,2);
    ntot = size(respirators,2);
    %B = regexp(BIRTH_Gender,'\S+','match');
    %T = cell2table([B{:}].');scsdsd
    counts=zeros(n,1);
    ic=1;
    respirator_str='Respirators:\t\t';
    for c = all_respirators
        counts(ic) = sum(strcmp(c,respirators));
        
        respirator_str=sprintf('%s%s: %d (%d%%), ',respirator_str,c{1,1},counts(ic),round(100*counts(ic)/ntot));
        ic = ic +1;
    end
    fprintf('%s (Missing: 0/%d)\n',respirator_str,ntot);
end
