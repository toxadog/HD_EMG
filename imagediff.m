function Idiff = imagediff(I,ang,step)
%Diff?renciation
ly=size(I,1);
lx=size(I,2);
Idiff=zeros(size(I));
dy=round(step*sin(ang));
dx=round(step*cos(ang));
if (dx>=0)&&(dy>=0)
    Idiff(dy+1:ly,dx+1:lx)=I(1:ly-dy,1:lx-dx)-I(dy+1:ly,dx+1:lx);
else if (dx<=0)&&(dy>=0)
        Idiff(dy+1:ly,1:lx+dx)=I(1:ly-dy,1-dx:lx)-I(dy+1:ly,1:lx+dx);
    else if (dx>=0)&&(dy<=0)
            Idiff(1:ly+dy,dx+1:lx)=I(1-dy:ly,1:lx-dx)-I(1:ly+dy,dx+1:lx);
            else if (dx<=0)&&(dy<=0)
                    Idiff(1:ly+dy,1:lx+dx)=I(1-dy:ly,1-dx:lx)-I(1:ly+dy,1:lx+dx);
                end
        end
    end
end
end

