function d=Cal_conv(A,M)
% convert Farsi Cal to Georgian
% A is string like 69.12: after point is tne month
% p is 4 or 12 is not need
% M is Boolean: M==0 means from Persian To Georgian
fyds=fix(A);    %--- First Year of Data Set
if ~exist('M','var')
    M=0;
end
if M==0
    if fyds<1300
        fyds=fyds+1300;
    end
    fqds=round(100*(A-fix(A)));
    fqds=fqds+3;
    
    if fqds>12
        fyds=fyds+1;
        fqds=fqds-12;
    end
    fyds=fyds+621;
    d=datenum(fyds,fqds,28);
else % M==1
    fyds=str2double(datestr(A,'yyyy'));
    fyds=fyds+621;
    %--- First Month of Data Set
    fqds=str2double(datestr(A,'mm'));
    fqds=fqds-3;
    if fqds<0
        fyds=fyds-1;
        fqds=fqds+12;
    end
    
    
    
end

end