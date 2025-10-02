function time_events=identify_desatbradytachy(tld)
%
%Function to identify desats, brady and tachy according to the definitions
%used in Poppi.
%note code assumes sampling rate of 1Hz.

time_secs=tld.time_ref*3600; 
no_transfusions=size(tld.hr,2);

%preallocate
clear time_events
time_events = struct;
time_events.desat={};
time_events.brady={};
time_events.tachy={};
counter = 1;
%%

for nt=1:no_transfusions
    
    HR=tld.hr(:,nt);
    sats=tld.sats(:,nt);

    % an episode of bradycardia will be defined as a pulse rate < 100 beats/min for at least 15 s
    
    count=1;
    
    all_st_brady_heartrate=[]; all_en_brady_heartrate=[];
    
    while count<length(HR)
        
        brady_f=find(HR(count:end)<100,1,'first');
        
        if ~isempty(brady_f)
            
            en_brady=find(HR((count+brady_f-1):end)>=100,1,'first');
            
            if (time_secs(count+en_brady+brady_f-2)-time_secs(count+brady_f-1))>=15
                
                all_st_brady_heartrate=[all_st_brady_heartrate,(count+brady_f-1)];
                all_en_brady_heartrate=[all_en_brady_heartrate,(count+brady_f+en_brady-2)];
                
                count=count+brady_f+en_brady;
                
                k=1;
    
                while k==1
                    
                    brady_g=find(HR(count:end)<100,1,'first');
                    
                    if (time_secs(count+brady_g-2)-time_secs(count-1))<=60
                        count=count+brady_g;
                    else
                        k=2;
                    end
                end
                
            else 
              
                
                    count=count+brady_f+en_brady; 
                
            end
            
            
            
        else
            count=length(HR);
            
        end
    end
    
    % an episode of tachycardia will be defined as a pulse rate > 200 beats/min for at least 15 s
    
    count=1;
    
    all_st_tachy_heartrate=[]; all_en_tachy_heartrate=[]; 
    
        
    while count<length(HR)
        
        tachy_f=find(HR(count:end)>200,1,'first');
        
        if ~isempty(tachy_f)
            
            en_tachy=find(HR((count+tachy_f-1):end)<=200,1,'first');
            
            if (time_secs(count+en_tachy+tachy_f-2)-time_secs(count+tachy_f-1))>=15
                
                all_st_tachy_heartrate=[all_st_tachy_heartrate,(count+tachy_f-1)];
                all_en_tachy_heartrate=[all_en_tachy_heartrate,(count+tachy_f+en_tachy-2)];
                
                count=count+tachy_f+en_tachy;
                
                k=1;
    
                while k==1
                    
                    tachy_g=find(HR(count:end)>200,1,'first');
                    
                    if (time_secs(count+tachy_g-2)-time_secs(count-1))<=60
                        count=count+tachy_g;
                    else
                        k=2;
                    end
                end
                
            else 
              
                
                 count=count+tachy_f+en_tachy; 
            
                   
            end
            
            
            
        else
            count=length(HR);
            
        end
    end
        
    

    
    % an episode of desaturation will be defined as oxygen saturation < 80% for at least 10 s
    
    count=1;
    
    all_st_desat=[]; all_en_desat=[]; 
    
    while count<length(sats)
        
        desat_f=find(sats(count:end)<80,1,'first');
        
        if ~isempty(desat_f)
            
            en_desat=find(sats((count+desat_f-1):end)>=80,1,'first');
            
            if (time_secs(count+en_desat+desat_f-2)-time_secs(count+desat_f-1))>=10
                
                all_st_desat=[all_st_desat,(count+desat_f-1)];
                all_en_desat=[all_en_desat,(count+desat_f+en_desat-2)];
                
                
                count=count+desat_f+en_desat;
            
                k=1;
    
                while k==1
    
                    desat_g=find(sats(count:end)<80,1,'first');
    
                    if (time_secs(count+desat_g-2)-time_secs(count-1))<=60 %if less than this apart count as one desat
                        count=count+desat_g;
                    else
                        k=2;
                    end
                end
                
            else 
              
                
                 count=count+desat_f+en_desat; 
            
                    
            end
                
            
            
        else
            count=length(sats);
            
        end
    end
    
    %% check signal quality during episodes. remove those where there are nans in the signal
    % find desats too near start of signal to check signal quality in baseline
    if ~isempty(all_st_desat)
    if all_st_desat(1)<21
       all_st_desat(1)=[]; all_en_desat(1)=[];
    end
    if all_st_desat(end)+21>length(sats)
       all_st_desat(end)=[]; all_en_desat(end)=[];
    end
    end
    
    %loop through all desats to check those with nan 20 seconds before and
    %after start
    if ~isempty(sats)
    i=1;
    while i<=length(all_st_desat)
        f=length(find(isnan(sats(all_st_desat(i)-20:all_st_desat(i)+20))));
        if f>0
            all_st_desat(i)=[]; all_en_desat(i)=[];
        else
            i=i+1;
        end
    end
    end
    
    %% go through brady
    all_st_brady=all_st_brady_heartrate;
    all_en_brady=all_en_brady_heartrate;
    
    if ~isempty(all_st_brady)
    if all_st_brady(1)<21
       all_st_brady(1)=[]; all_en_brady(1)=[];
    end
    end
    
    i=1;
    while i<=length(all_st_brady)
        f=length(find(isnan(HR(all_st_brady(i)-20:all_st_brady(i)+25))));
        if f>0
            all_st_brady(i)=[]; all_en_brady(i)=[];
        else
            i=i+1;
        end
    end
    
    
    
    %check tachy
    all_st_tachy=all_st_tachy_heartrate;
    all_en_tachy=all_en_tachy_heartrate;
    
    if ~isempty(all_st_tachy)
    if all_st_tachy(1)<21
       all_st_tachy(1)=[]; all_en_tachy(1)=[];
    end
    end
    
    
    i=1;
    while i<=length(all_st_tachy)
        f=length(find(isnan(HR(all_st_tachy(i)-20:all_st_tachy(i)+25))));
        if f>0
            all_st_tachy(i)=[]; all_en_tachy(i)=[];
        else
            i=i+1;
        end
    end
    
    
    start_time_desat=time_secs(all_st_desat)/3600; %time in hours
    start_time_brady=time_secs(all_st_brady)/3600;
    start_time_tachy=time_secs(all_st_tachy)/3600;
    
    if ~isempty(start_time_desat)
        time_events.desat{counter}=start_time_desat;
    else
        time_events.desat{counter}=NaN;
    end
    if ~isempty(start_time_brady)
        time_events.brady{counter}=start_time_brady;
    else
        time_events.brady{counter}=NaN;
    end
    if ~isempty(start_time_tachy)
        time_events.tachy{counter}=start_time_tachy;
    else
        time_events.tachy{counter}=NaN;
    end

    counter=counter+1;
end
end
