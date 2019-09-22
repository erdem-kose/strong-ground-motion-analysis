function arlData = extractARLdata(arlSettings,avxData)  
    g2cms2=980.6;%cm/s2
    
    N=size(avxData,2);
    
    arlData.coord=zeros(3,N);
    arlData.epi_coord=zeros(3,N);
    acc=cell(1,N);
    
    arlData.N=-inf;
    for i=1:N       
        arlData.coord(:,i)=ll2xy([avxData{i}.longi;  avxData{i}.latit]);
        arlData.epi_coord(:,i)=ll2xy([avxData{i}.epi_longi; avxData{i}.epi_latit]);
        
        if avxData{i}.N>arlData.N
            arlData.N=avxData{i}.N;
        end
        
        arlData.stationName{i}=avxData{i}.stationName;
    end

    arlData.type=arlSettings.type;
    arlData.class='NA';
    
    arlData.zone=zone_create(min(xy2ll(arlData.coord-arlSettings.range),[],2),max(xy2ll(arlData.coord+arlSettings.range),[],2),2,arlData.N);
    arlData.zone=ll2xy(arlData.zone);

    ML=avxData{1}.ML;

    switch arlSettings.type
        case 'campbell'
            arlData.PHA=campbell(arlData,arlData.epi_coord(:,1),ML).*g2cms2;
        case 'boore'
            arlData.class=arlSettings.boore_site_class;
            arlData.PHA=boore(arlData,arlData.epi_coord(:,1),ML).*g2cms2;
        otherwise
            arlData.PHA=campbell(arlData,arlData.epi_coord(:,1),ML).*g2cms2;
    end
    
    arlData.zone=xy2ll(arlData.zone);
    arlData.coord=xy2ll(arlData.coord);
end

function [PHAvals,R]=campbell(arlData,coord,ML) %numerical indefinite integral as filter
    R=sqrt(sum((arlData.zone-coord).^2,1));
    PHAvals=abs(exp(-4.141+0.868.*ML-1.09.*log(R+0.0606*exp(0.7*ML))+0.37));
end

function [PHAvals,R]=boore(arlData,coord,ML) %numerical indefinite integral as filter
    Gb=0; Gc=0;
    switch arlData.class
        case 'A'
            Gb=0; Gc=0;
        case 'B'
            Gb=1; Gc=0;
        case 'C'
            Gb=0; Gc=1;
    end
    
    R=sqrt(sum((arlData.zone-coord).^2,1));
    PHAvals=abs(10.^((-0.038+0.216.*(ML-6)-0.777.*log10(sqrt(R.^2+30.03))+0.158*Gb+0.254*Gc)));
end

function x=zone_create(range_min,range_max,d,N)
    nV=ceil(N^(1/d));

    v=linspace(0,1,nV);
    [Y{d:-1:1}] = ndgrid(1:(nV));

    i = reshape(cat(d+1,Y{:}),[],d);
    x = v(i)';

    for i=1:d
        x_min=x(i,1); x_max=x(i,end);
        if (x_min~=range_min(i)) || (x_max~=range_max(i))
            x(i,:)=(x(i,:)).*(range_max(i)-range_min(i))+range_min(i);
        end
    end
end

function xy=ll2xy(latlon)
    deg2rad=pi/180;

%     a=6378137;
%     b=6356752;
%     e2=((a.^2)-(b.^2))./(eps+(a.^2));
%     N=a./(eps+sqrt(1-e2.*sin(latlon(1).*deg2rad)));
%     N=N/1000;
    N=6391;
    
    xy=zeros(3,size(latlon,2));
    xy(1,:) = N .* cos(latlon(1,:)*deg2rad) .* cos(latlon(2,:)*deg2rad);
    xy(2,:) = N .* cos(latlon(1,:)*deg2rad) .* sin(latlon(2,:)*deg2rad);
    xy(3,:) = N .* sin(latlon(1,:)*deg2rad);
end

function latlon=xy2ll(xy)
    rad2deg=180/pi;

    N=6391;
    latlon(2,:)=atan2(xy(2,:),xy(1,:)).*rad2deg;
    latlon(1,:)=asin(xy(3,:)/N)*rad2deg;
end