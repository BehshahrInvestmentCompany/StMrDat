function Cal_conv(A,p)
% convert Farsi Cal to Georgian
% A is string like 69.1
% p is 4 or 12
fyds=fix(Date(1));    %--- First Year of Data Set
if fyds>1300
    fyds=fyds-1300;
end
fqds=round(10*(Date(1)-fix(Date(1))));

if p==4
    fqds=fqds-1;
    if fqds<1
        fyds=fyds+1;
    end
    fyds=fyds+621;
    datenum(fyds,fqds*3)
elseif p==12
    
end


end