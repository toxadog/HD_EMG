function Point = get_points(MAX,MIN,SIZE,ang,step)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
MID=abs(round((MAX-MIN)/2));
XD=MID(1);
YD=MID(2);
Point=[];
if (cos(ang)>=0)&&(sin(ang)>=0)
    Xmax=SIZE(1);
    Ymax=SIZE(2);
    Xmin=0;
    Ymin=0;
else if (cos(ang)<=0)&&(sin(ang)>=0)
        Xmax=0;
        Ymax=SIZE(2);
        Xmin=SIZE(1);
        Ymin=0;
    else if (cos(ang)<=0)&&(sin(ang)<=0)
            Xmax=0;
            Ymax=0;
            Xmin=SIZE(1);
            Ymin=SIZE(2);
        else
            Xmax=SIZE(1);
            Ymax=0;
            Xmin=0;
            Ymin=SIZE(2);
        end
    end
end
            
while (sign(cos(ang))*(Xmax-XD)>=0)&&(sign(sin(ang))*(Ymax-YD)>=0)
    Point=[Point;XD YD];
    XD=round(XD+cos(ang)*step);
    YD=round(YD+sin(ang)*step);    
end
XD=MID(1)-cos(ang)*step;
YD=MID(2)-sin(ang)*step;
while (sign(cos(ang))*(XD-Xmin)>=0)&&(sign(sin(ang))*(YD-Ymin)>=0)
    Point=[XD YD;Point];
    XD=round(XD-cos(ang)*step);
    YD=round(YD-sin(ang)*step);    
end

end

