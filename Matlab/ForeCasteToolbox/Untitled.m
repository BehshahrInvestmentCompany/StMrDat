for i=1:horizon
     Actual_date(nd+i)=Actual_date(nd+i-1)+0.01;
     if Actual_date(nd+i)-fix(Actual_date(nd+i))>0.115
         Actual_date(nd+i)=fix(Actual_date(nd+i))+1;
     end
end
% mode