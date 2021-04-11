function dist=great_circle_distance(coor1, coor2)
%[xlat_sub(ni) ylon_sub(ni)],[xlat_sub(nj) ylon_sub(nj)]);

%coor1=[xlat_sub(ni) ylon_sub(ni)];
%coor2=[xlat_sub(nj) ylon_sub(nj)];

xlat1=coor1(:,1);
ylon1=coor1(:,2);

xlat2=coor2(:,1);
ylon2=coor2(:,2);

%{
Earth_radius=6371.0087714; %WGS84 mean radius

dist_lon=abs(ylon1 - ylon2).*(abs(ylon1 - ylon2)<=180)+...
            (360-abs(ylon1 - ylon2)).*((abs(ylon1 - ylon2))>180);


dist = Earth_radius*acos(sin(xlat1*pi/180).*sin(xlat2*pi/180) ...
    + cos(xlat1*pi/180).*cos(xlat2*pi/180).*cos(dist_lon*pi/180));

dist = real(dist);
%}

dist = sqrt((xlat1-xlat2).^2+(ylon1-ylon2).^2);