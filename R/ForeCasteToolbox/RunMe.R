
#graphics.off()
cat("\014") # Clear Screen
rm(list=ls()) #remove all variable
setwd(dirname(sys.frame(1)$ofile))
{
Target_var_names=c('CPI' , 'Gold' , 'USD'  , 'SI','GDP')
Dum_var_name=c('D1','D2')
FRT<-1
nd<-20
Train.Start<-c(1390,2)
horizon<-4
} # Use Inputs
#options(echo=FALSE)
# require("data.table") '''
if (Sys.getenv("JAVA_HOME")!="")
{
  Sys.setenv(JAVA_HOME='C:\\Program Files\\MATLAB\\R2016a\\sys\\java\\jre\\win64\\jre')
}# Load Java Path
{
require("rJava")
require("tseries")
require("xlsx")
require("forecast")
require("graphics")
require("stats")
require("tsDyn")

} #Load Libraries
#Adr=getwd ();
Data0<- read.xlsx("Input/Data.xlsx",'SHEET','Data')
#script.dir <- dirname(sys.frame(1)$ofile)
{
  Transformation <- function(Target, Exp_Var,T,Target_var_name)
    {
    # This Function Do Transformation Base On Below
    # 'no_change', 'Ln', 'Diff', 'Double_Diff', 'Diff_Ln', 'Double_Diff_Ln'
    #     0         1       2         3             4           5
    
   # Target_=T[]
    Tr=T[,Target_var_name]
    if (Tr==1){Target_=log(Target)} # Logarith
    else if (Tr==2){Target_=diff(Target)}# Diff
    else if (Tr==3){Target_=diff(Target,2)}# Double Diff
    else if (Tr==4){Target_=diff(log(Target),1)} #Log Diff
    else if (Tr==5){Target_=diff(log(Target),2)} # Log Double Diff
  #  else if (Tr=6){Target_=log(data_adj(5:T,i))-log(data_adj(1:T-4,i))}
    else {Target_=Target}
    
    Exp_Var_=Exp_Var
    for (Vvv in Exp_Var_names)
    { # Deseasonal Data
      Tr=T[,Vvv]
      if (Tr==1){Exp_Var_[,Vvv]=log(Exp_Var[,Vvv])} # Logarith
      else if (Tr==2){Exp_Var_[,Vvv]=diff(Exp_Var[,Vvv])}# Diff
      else if (Tr==3){Exp_Var_[,Vvv]=diff(Exp_Var[,Vvv],2)}# Double Diff
      else if (Tr==4){Exp_Var_[,Vvv]=diff(log(Exp_Var[,Vvv]),1)} #Log Diff
      else if (Tr==5){Exp_Var_[,Vvv]=diff(log(Exp_Var[,Vvv]),2)} # Log Double Diff
      #  else if (Tr=6){Target_=log(data_adj(5:T,i))-log(data_adj(1:T-4,i))}
      else {Exp_Var_=Exp_Var} 
    }
    
    return(list(T=Target_,E=Exp_Var_))
  }  
} # User Functions

if( FRT==1)
{ # Check Wheter To load First row as Transformation
  T=Data0[1,]
  Data0 = Data0[-1,]
  Frq<-T$Date[1] # frequency
  T=T[,-which(names(Data0) %in% c("Date",Dum_var_name))]
}

Data<-ts(data=Data0,frequency=Frq,start =  c(round(Data0$Date[1]),Frq*100/12*(Data0$Date[1]-round(Data0$Date[1]))))#Data0$Date[1]
Tvarnameindx<-1
#for (Tvarnameindx in 1:length(Target_var_names))
{ # The Loop Of Target Variable 
  Model_Count<-1
  model=list();
  Target=Data[,Target_var_names[Tvarnameindx]];
  ddd<-decompose(Target,type = "additive")
  Target<-Target-ddd$seasonal
  Dummy=Data[,Dum_var_name]
  Exp_Var=Data[,-which(colnames(Data) %in% c("Date",Target_var_names[Tvarnameindx],Dum_var_name))]
  Exp_Var_names=colnames(Exp_Var)
  for (Vvv in Exp_Var_names)
  { # Deseasonal Data
    ddd<-decompose(Exp_Var[,Vvv],type = "additive")
    Exp_Var[,Vvv]<-Exp_Var[,Vvv]-ddd$seasonal 
  }
 
  # Biul Banker of Exogenius Combination
  C=combn(1:length(Exp_Var_names),3)
  # Biuld Banker of Transformation Name 
  Trans_mod=cbind('no_change', 'Ln', 'Diff', 'Double_Diff', 'Diff_Ln', 'Double_Diff_Ln')
  
  #xx <- window(x, start = 1988)
  Tt2=floor(log(T[1,1])/log(10))+1;
  if  (!is.finite (Tt2))
  {
  Tt2=1
  }
  
  for (tt in 1:Tt2)
  {
    T[tt+1,]<- (T[1,]%%10^tt-T[1,]%%10^(tt-1))/10^(tt-1)
  }
  T=T[-1,]# Remove First Row
  tt<-1
 # for (tt in 1:Tt2)
  {
  # Loop For Transformation
  TT<-Transformation(Target = Target,Exp_Var = Exp_Var ,T=T[tt,],Target_var_name = Target_var_names[Tvarnameindx] )
  Target_=TT$T
  Exp_Var_=TT$E
  for (i in 1:nd){# Training Loop
    Ynd=nd%/%Frq
    Qnd=nd%%Frq
    Target__<-window(Target_,end=Train.Start+c(Ynd,Qnd))
    Exp_Var__<-window(Exp_Var_,end=Train.Start+c(Ynd,Qnd))
    Dummy__<-window(Dummy,end=Train.Start+c(Ynd,Qnd))
    
    model<-list(auto.arima(Target__),
                aar(Target__, m=1, d=1, steps=horizon),
                aar(Target__, m=1, d=2, steps=horizon),
                aar(Target__, m=1, d=3, steps=horizon),
                star(Target__, m=2, noRegimes=3, d = 1, steps = horizon, rob = FALSE,thDelay=1, sig=0.05),
                setar(Target__,),
                lstar(Target__,),
                nnetTs(),
                TVAR(),
                lineVar(),
                
                )
    Model_Count<-Model_Count+1
  }# Training Loop
  }# Loop For Transformation
  Model_Count=length(model);
  
} # The Loop Of Target Variable
