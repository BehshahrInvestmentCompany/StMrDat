function d=Cal_conv(A,M)
% convert Farsi Cal to Georgian
% A is string like 69.12: after point is tne month
% p is 4 or 12 is not need
% M is Boolean: M==0 means from Persian To Georgian

if ~exist('M','var')
    M=0;
end
if M==0
    fyds=fix(A);    %--- First Year of Data Set
    if fyds<100
        fyds=fyds+1300;
    end
    fqds=round(100*(A-fix(A)));
    
    DayCount=fyds*365+fix(fyds/4)+fqds*31-max(0,fqds-6)-max(0,fqds-11);
    
    % The Kabise Year
    DayCount((fqds==12) .* (mod(fyds+1,4)==0)==1 )=DayCount((fqds==12) .* (mod(fyds+1,4)==0)==1 )+1;
    
    % March 21, 1921 Or 621-03-10 == 001/1/0
    d=datenum(621,03,10)+DayCount;
    %datestr(d)
else % M==1
    A=A-datenum(621,03,10);
    Y=fix(A/365);
    TY=fix(abs(Y/4)/365);
    Y=Y-TY;
    A=A-Y*365-fix(Y/4);
    Y(A<=0)=Y(A<=0)-1;
    A(A<=0)=A(A<=0)+366;
    M=nan(size(A));
    for i=1:length(A)
        if A(i)<=(31*6)
           M(i)=fix(A(i)/31);
        elseif A(i)<=(31*6+30*5)
            M(i)=fix((A(i)-31*6)/30)+6;
        else
            M(i)=12;
            if A(i)>(31*6+30*6)
                warning('Moth Problem');
            end
        end
    end
    d=Y+M/100;
    
    
    
end

end