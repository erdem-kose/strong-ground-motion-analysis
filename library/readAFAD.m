function eqData=readAFAD(filename)
    %1 gal ~ 1 cm/s2
    g2cms2=980.6;%cm/s2

    output_folder_name=['outputs\' filename];
    if exist(output_folder_name,'dir')==0
        mkdir(output_folder_name);
    else
        delete([output_folder_name '\*']);
    end
    
    mat_file_name=['data\' filename '\' filename '.mat'];
    if exist(mat_file_name,'file')~=0
        load(mat_file_name);
        return;
    end
    
    dinfo = dir(['data\' filename '\*.txt']);
    minRecTime=inf; maxN=-inf;
    for k = 1 : length(dinfo)
        active_filename = dinfo(k).name;  %just the name

        fid = fopen(active_filename);
        tline = fgetl(fid);
        i=1;

        disp(['Record ID:' num2str(k)]);
        while ischar(tline)
            tline = fgetl(fid);
            i=i+1;

            if i==2
                strp=find(tline==':')+1;
                endp=strp+15;
                eqData{k}.stationName=(tline(strp:endp));
            elseif i==4
                coordRaw=(tline((find(tline==':')+1):(end)));
                strPos=find(isletter(coordRaw));
                
                eqData{k}.epi_latit=str2num(coordRaw(1:(strPos(1)-1)));
                if coordRaw(strPos(1))=='S'
                    eqData{k}.epi_latit=-eqData{k}.epi_latit;
                end
                eqData{k}.epi_longi=str2num(coordRaw((strPos(1)+2):(strPos(2)-1)));
                if coordRaw(strPos(2))=='W'
                    eqData{k}.epi_longi=-eqData{k}.epi_longi;
                end
            elseif i==6
                eqData{k}.ML=str2num((tline((find(tline==':')+1):(end-2))));
            elseif i==8
                coordRaw=(tline((find(tline==':')+1):(end)));
                strPos=find(isletter(coordRaw));
                
                eqData{k}.latit=str2num(coordRaw(1:(strPos(1)-1)));
                if coordRaw(strPos(1))=='S'
                    eqData{k}.latit=-eqData{k}.latit;
                end
                eqData{k}.longi=str2num(coordRaw((strPos(1)+2):(strPos(2)-1)));
                if coordRaw(strPos(2))=='W'
                    eqData{k}.longi=-eqData{k}.longi;
                end
            elseif i==12
                end_ind=strfind(tline,'(GMT)')-1;
                str_ind=strfind(tline,' : ')+14;
                
                recTime=tline(str_ind:end_ind);
                eqData{k}.recStartTime=(str2num(recTime(1:2))*60+str2num(recTime(4:5)))*60+str2num(recTime(7:end));
            elseif i==13
                eqData{k}.N=str2num((tline((find(tline==':')+1):end)));
            elseif i==14
                eqData{k}.fs=1/(eps+str2num((tline((find(tline==':')+1):end))));
            elseif i<=17
                disp(tline);
            elseif i>=19
                if ischar(tline)==1
                    arr=str2num(tline);
                    eqData{k}.EQ_NS(i-18)=arr(1);
                    eqData{k}.EQ_EW(i-18)=arr(2);
                    eqData{k}.EQ_UD(i-18)=arr(3);
                end
            end
        end
        disp(' ');
        fclose(fid);

        eqData{k}.EQ_NS=eqData{k}.EQ_NS;
        eqData{k}.EQ_EW=eqData{k}.EQ_EW;
        eqData{k}.EQ_UD=eqData{k}.EQ_UD;
        ds_h_dinfo = dir(['data\' filename '\*_ds_h.csv']);
        ds_h = csvread(ds_h_dinfo(1).name,1,0);
        eqData{k}.T_h=ds_h(:,1);
        eqData{k}.ds_h=ds_h(:,2).*g2cms2;
        ds_v_dinfo = dir(['data\' filename '\*_ds_v.csv']);
        ds_v = csvread(ds_v_dinfo(1).name,1,0);
        eqData{k}.T_v=ds_v(:,1);
        eqData{k}.ds_v=ds_v(:,2).*g2cms2;
        
        if eqData{k}.recStartTime<minRecTime
            minRecTime=eqData{k}.recStartTime;
        end
        recN=eqData{k}.N;
        if recN>maxN
            maxN=recN;
        end
    end
    
    for k = 1 : length(dinfo)
        eqData{k}.recStartTime=eqData{k}.recStartTime-minRecTime;
        eqData{k}.t=0:(1/eqData{k}.fs):((eqData{k}.N-1)/eqData{k}.fs);
        eqData{k}.t_full=0:(1/eqData{k}.fs):((maxN-1)/eqData{k}.fs);
        eqData{k}.t=eqData{k}.t+eqData{k}.recStartTime;
    end
    
    save(mat_file_name);
end