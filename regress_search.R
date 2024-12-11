
if(!exists("only.fixef"))
  library(lme4)
if(exists("only.fixef"))
 if(!only.fixef)
  library(lme4)


regress.ic.dredge=function(data,response,covs,family="normal",
 offset="",strata="",use.glmmADMB=F, use.glmmTMB=F,
 IC="BIC",threshold="flexible",
 check.est=F,check.se=F, do.return.status=F, talkative=T,
 top=1, max.complex=NULL)
{
  numcov=length(covs)
  
  if(IC=="AICc")
    library(MuMIn)

  if(use.glmmADMB)
    library("glmmADMB")

  if(use.glmmTMB)
    library("glmmTMB")

  if(family=="clogit")
    library("survival")

  if(family=="ordinal" | family=="ordinal_cauchit")
    library("ordinal")
  
  # Assume random factors start with "(":
  isrand=substr(covs,1,1)=="("  

  numtop=0
  top.res=list()
  top.ic=rep(0,0)
  top.statuses=rep(0,0)
  if(top>1)
    do.return.status=TRUE
  
  # Create regression string:
  getmodel=function(status)
  {
   if(sum(status)==0)
    {
     ret=sprintf("%s~1",response)  
    }
   else
    {
       # only random effects?
       if((family=="ordinal" | family=="ordinal_cauchit") & 
          sum(status==1)==sum(isrand[status==1])) 
       {
         ret=sprintf("%s~1+",response)
       }
       else
       {
         ret=sprintf("%s~",response)
       }
       index=which(status==1)
       ret=sprintf("%s%s",ret,covs[index[1]])
       if(length(index)>1)
        for(i in index[2:length(index)])
         ret=sprintf("%s+%s", ret, covs[i]) 
    }
   if(offset!="")
    ret=sprintf("%s+%s",ret,offset)

   if(strata!="")
    ret=sprintf("%s+strata(%s)",ret,strata)
  
   return(ret)
  }
  
  
  # Make zero model:
  if(family=="normal")
  {
    if(offset!="")
      res0=lm(as.formula(sprintf("%s~1+%s",response,offset)),data=data)
    if(offset=="")
      res0=lm(as.formula(sprintf("%s~1",response)),data=data)
  }
  if(family=="clogit")
  {
    if(strata!="")
      res0=clogit(as.formula(sprintf("%s~%s+strata(%s)",response,covs[1],strata)),data=data)
    if(strata=="")
      res0=clogit(as.formula(sprintf("%s~%s",response,covs[1])),data=data)
  }
  if(family=="ordinal")
  {
    res0=clm(as.formula(sprintf("%s~1",response)),data=data,threshold=threshold)
  }
  if(family=="ordinal_cauchit")
  {
    res0=clm(as.formula(sprintf("%s~1",response)),data=data,threshold=threshold,link="cauchit")
  }
  if(family!="normal" & family!="clogit" & family!="ordinal" & family!="ordinal_cauchit")
  {
    if(offset!="")
      res0=glm(as.formula(sprintf("%s~1+%s",response,offset)),data=data,family=family)
    if(offset=="")
      res0=glm(as.formula(sprintf("%s~1",response)),data=data,family=family)
  }

  # Run regression:
  getres=function(status)
  {
    str=getmodel(status)
    if(family=="normal")
    {
      if(use.glmmADMB)
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
               family="gaussian") ), 
               abort=function() {return(res0)}))
	       
      if(use.glmmTMB)
        return(withRestarts(tryCatch(glmmTMB(as.formula(str), data=data, 
               family="gaussian") ), 
               abort=function() {return(res0)}))
	       
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(lmer(as.formula(str), data=data, REML=F) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(lm(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }
    
    if(family=="clogit")
    {
      return(withRestarts(tryCatch(clogit(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }

    if(family=="ordinal")
    {
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(clmm(as.formula(str), data=data,threshold=threshold) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(clm(as.formula(str), data=data,threshold=threshold) ), 
             abort=function() {return(res0)}))
    }
    
    if(family=="ordinal_cauchit")
    {
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(clmm(as.formula(str), link="cauchit", data=data,threshold=threshold) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(clm(as.formula(str), link="cauchit", data=data,threshold=threshold) ), 
             abort=function() {return(res0)}))
    }
    
    if(use.glmmADMB)
      {
        #print.srcref(sprintf("glmmadmb(%s, data=data, family=%s)", str,family))
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
             family=family) ), 
             abort=function() {return(res0)}))
      }
 	      
    if(use.glmmTMB)
      {
        #print.srcref(sprintf("glmmadmb(%s, data=data, family=%s)", str,family))
        return(withRestarts(tryCatch(glmmTMB(as.formula(str), data=data, 
             family=family) ), 
             abort=function() {return(res0)}))
      }
 	      
    if(sum(isrand[status!=0])>0)
    {
      print.srcref(sprintf("glmer(as.formula(%s), data=data, family=%s)",str,family))
      return(withRestarts(tryCatch(glmer(as.formula(str), data=data,family=family) ), 
             abort=function() {return(res0)}))
    }
    
    return(withRestarts(tryCatch(glm(as.formula(str), data=data,family=family) ), 
           abort=function() {return(res0)}))
  }
  

  update.top=function(new.ic,new.status,new.res,top,numtop,top.ic,top.statuses,top.res)
  {
    if(is.na(new.ic) | is.infinite(new.ic))
      return(list(numtop=numtop,top.ic=top.ic,top.statuses=top.statuses,top.res=top.res))
   
    if(numtop>0)
    {
      #show(top.ic)
      if(new.ic<top.ic[numtop])
      {
        equalstatus=function(x) 
           sum(new.status==x)==numcov      
        
        if(sum(apply(top.statuses,1,equalstatus))==0) # not in already?
        {
          #show(sprintf("%f<%f",new.ic,top.ic[numtop]))
          top.ic=c(top.ic,new.ic)
          top.statuses=rbind(top.statuses,new.status)
          numtop=length(top.ic)
          #show(numtop)
          top.res[[numtop]]=new.res
          
	  # Order:
          o=order(top.ic)
          if(numtop>top)
          {
            o=o[1:top]
            numtop=top
            #show("numtop2")
            #show(numtop)
          }
          top.ic=top.ic[o]
          top.statuses=top.statuses[o,]
          top.res2=list()
          for(i in 1:top)
            top.res2[[i]]=top.res[[o[i]]]
          top.res=top.res2

        }  
      }    
    }
    if(numtop==0)
    {
      top.res[[numtop+1]]=new.res
      top.ic=c(new.ic)
      top.statuses=rbind(new.status)
      numtop=1
      #show(numtop)
    } 

   return(list(numtop=numtop,top.ic=top.ic,top.statuses=top.statuses,top.res=top.res))
  }

  increment.status=function(status.old)
  {
    status.new=status.old
    numcov=length(status.old)
    i=numcov
    status.new[i]=1-status.new[i]
    while(status.new[i]==0 & i>1)
    {
      i=i-1
      status.new[i]=1-status.new[i]
    }

   status.new
  }

  status1=rep(0,numcov)
  
  for(i in 1:(2^numcov))
  {
   dobreak=FALSE
   if(!is.null(max.complex))
    if(sum(status1)>max.complex)
     dobreak=TRUE  

   if(!dobreak)
   {
    res1=getres(status1)
    ic1=BIC(res1)
    if(IC=="AIC")
     ic1=AIC(res1)
    if(IC=="AICc")
     ic1=AICc(res1)
    if(IC!="BIC" & IC!="AIC" & IC!="AICc")
     ic1=AIC(res1,k=as.numeric(IC))
    if(talkative)
      show(c(i,ic1))
    if(check.est)
     if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
    if(check.est)
     if(sum(is.na(summary(res1)$coef[,2]))>0)
       ic1=NA

    if(sum(status1)==0)
    {
      ic0=ic1
      res0=res1

      if(!is.na(ic0))
      {
       best.ic=ic0
       best.status=status1
      }
    }

    if(!is.na(ic1))
    {
     if(is.na(ic0))
     {
      best.ic=ic1
      best.status=status1
     }

     if(!is.na(ic0))
      if(ic1<best.ic)
     {
      best.ic=ic1
      best.status=status1
     }
    }

    if(top>1)
    {
     utop=update.top(ic1,status1,res1,top,numtop,top.ic,top.statuses,top.res)
     numtop=utop$numtop
     top.ic=utop$top.ic
     top.statuses=utop$top.statuses
     top.res=utop$top.res
    }
   }
   
   status1=increment.status(status1)
  }
  status=best.status

  if(top>1)
    return(list(status=status, res=getres(status),top.res=top.res,numtop=numtop,top.statuses=top.statuses,top.ic=top.ic,top.weights=exp(-(top.ic-min(top.ic))/2)/sum(exp(-(top.ic-min(top.ic))/2))  ))
  if(!do.return.status)
    return(getres(status))
  return(list(status=status, res=getres(status)))
}



regress.functioncriterion.dredge=function(data,response,covs,family="normal",fun=BIC,
 do.maximize=F,offset="",strata="",use.glmmADMB=F, threshold="flexible",
 check.est=F,check.se=F, do.return.status=F, start.status=NA, talkative=T,
 top=1, max.complex=NULL)
{
  numcov=length(covs)
  
  if(use.glmmADMB)
    library("glmmADMB")

  if(family=="clogit")
    library("survival")

  if(family=="ordinal" | family=="ordinal_cauchit")
    library("ordinal")
  
  # Assume random factors start with "(":
  isrand=substr(covs,1,1)=="("  

  numtop=0
  top.res=list()
  top.fun=rep(0,0)
  top.statuses=rep(0,0)
  if(top>1)
    do.return.status=TRUE
  
  # Create regression string:
  getmodel=function(status)
  {
   if(sum(status)==0)
    {
     ret=sprintf("%s~1",response)  
    }
   else
    {
       # only random effects?
       if((family=="ordinal" | family=="ordinal_cauchit") & 
          sum(status==1)==sum(isrand[status==1])) 
       {
         ret=sprintf("%s~1+",response)
       }
       else
       {
         ret=sprintf("%s~",response)
       }
       index=which(status==1)
       ret=sprintf("%s%s",ret,covs[index[1]])
       if(length(index)>1)
        for(i in index[2:length(index)])
         ret=sprintf("%s+%s", ret, covs[i]) 
    }
   if(offset!="")
    ret=sprintf("%s+%s",ret,offset)

   if(strata!="")
    ret=sprintf("%s+strata(%s)",ret,strata)
  
   return(ret)
  }
  
  
  # Make zero model:
  if(family=="normal")
  {
    if(offset!="")
      res0=lm(as.formula(sprintf("%s~1+%s",response,offset)),data=data)
    if(offset=="")
      res0=lm(as.formula(sprintf("%s~1",response)),data=data)
  }
  if(family=="clogit")
  {
    if(strata!="")
      res0=clogit(as.formula(sprintf("%s~%s+strata(%s)",response,covs[1],strata)),data=data)
    if(strata=="")
      res0=clogit(as.formula(sprintf("%s~%s",response,covs[1])),data=data)
  }
  if(family=="ordinal")
  {
    res0=clm(as.formula(sprintf("%s~1",response)),data=data,threshold=threshold)
  }
  if(family=="ordinal_cauchit")
  {
    res0=clm(as.formula(sprintf("%s~1",response)),data=data,threshold=threshold,link="cauchit")
  }
  if(family!="normal" & family!="clogit" & family!="ordinal" & family!="ordinal_cauchit")
  {
    if(offset!="")
      res0=glm(as.formula(sprintf("%s~1+%s",response,offset)),data=data,family=family)
    if(offset=="")
      res0=glm(as.formula(sprintf("%s~1",response)),data=data,family=family)
  }

  # Run regression:
  getres=function(status)
  {
    str=getmodel(status)
    if(family=="normal")
    {
      if(use.glmmADMB)
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
               family="gaussian") ), 
               abort=function() {return(res0)}))
	       
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(lmer(as.formula(str), data=data, REML=F) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(lm(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }
    
    if(family=="clogit")
    {
      return(withRestarts(tryCatch(clogit(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }

    if(family=="ordinal")
    {
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(clmm(as.formula(str), data=data,threshold=threshold) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(clm(as.formula(str), data=data,threshold=threshold) ), 
             abort=function() {return(res0)}))
    }
    
    if(family=="ordinal_cauchit")
    {
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(clmm(as.formula(str), link="cauchit", data=data,threshold=threshold) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(clm(as.formula(str), link="cauchit", data=data,threshold=threshold) ), 
             abort=function() {return(res0)}))
    }
    
    if(use.glmmADMB)
      {
        #print.srcref(sprintf("glmmadmb(%s, data=data, family=%s)", str,family))
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
             family=family) ), 
             abort=function() {return(res0)}))
      }
 	      
    if(sum(isrand[status!=0])>0)
    {
      print.srcref(sprintf("glmer(as.formula(%s), data=data, family=%s)",str,family))
      return(withRestarts(tryCatch(glmer(as.formula(str), data=data,family=family) ), 
             abort=function() {return(res0)}))
    }
    
    return(withRestarts(tryCatch(glm(as.formula(str), data=data,family=family) ), 
           abort=function() {return(res0)}))
  }
  

  update.top=function(new.fun,new.status,new.res,top,numtop,top.fun,top.statuses,top.res)
  {
    if(is.na(new.fun) | is.infinite(new.fun))
      return(list(numtop=numtop,top.fun=top.fun,top.statuses=top.statuses,top.res=top.res))
   
    if(numtop>0)
    {
      #show(top.fun)
      if(new.fun<top.fun[numtop])
      {
        equalstatus=function(x) 
           sum(new.status==x)==numcov      
        
        if(sum(apply(top.statuses,1,equalstatus))==0) # not in already?
        {
          #show(sprintf("%f<%f",new.fun,top.fun[numtop]))
          top.fun=c(top.fun,new.fun)
          top.statuses=rbind(top.statuses,new.status)
          numtop=length(top.fun)
          #show(numtop)
          top.res[[numtop]]=new.res
          
	  # Order:
          o=order(top.fun)
          if(numtop>top)
          {
            o=o[1:top]
            numtop=top
            #show("numtop2")
            #show(numtop)
          }
          top.fun=top.fun[o]
          top.statuses=top.statuses[o,]
          top.res2=list()
          for(i in 1:top)
            top.res2[[i]]=top.res[[o[i]]]
          top.res=top.res2

        }  
      }    
    }
    if(numtop==0)
    {
      top.res[[numtop+1]]=new.res
      top.fun=c(new.fun)
      top.statuses=rbind(new.status)
      numtop=1
      #show(numtop)
    } 

   return(list(numtop=numtop,top.fun=top.fun,top.statuses=top.statuses,top.res=top.res))
  }

  increment.status=function(status.old)
  {
    status.new=status.old
    numcov=length(status.old)
    i=numcov
    status.new[i]=1-status.new[i]
    while(status.new[i]==0 & i>1)
    {
      i=i-1
      status.new[i]=1-status.new[i]
    }

   status.new
  }

  status1=rep(0,numcov)
  
  for(i in 1:(2^numcov))
  {
   dobreak=FALSE
   if(!is.null(max.complex))
    if(sum(status1)>max.complex)
     dobreak=TRUE

   if(!dobreak)
   {
    res1=getres(status1)
    fun1=fun(res1)
    if(talkative)
      show(c(i,fun1))
    if(check.est)
     if(sum(is.na(summary(res1)$coef[,1]))>0)
      fun1=NA
    if(check.est)
     if(sum(is.na(summary(res1)$coef[,2]))>0)
      fun1=NA

    if(sum(status1)==0)
    {
     fun0=fun1
     res0=res1

     if(!is.na(fun0))
     {
      best.fun=fun0
      best.status=status1
     }
    }

    if(!is.na(fun1))
    {
     if(is.na(fun0))
     {
      best.fun=fun1
      best.status=status1
     }

     if(!is.na(fun0))
      if(fun1<best.fun)
     {
      best.fun=fun1
      best.status=status1
     }
    }

    if(top>1)
    {
     utop=update.top(fun1,status1,res1,top,numtop,top.fun,top.statuses,top.res)
     numtop=utop$numtop
     top.fun=utop$top.fun
     top.statuses=utop$top.statuses
     top.res=utop$top.res
    }
   }

   status1=increment.status(status1)
  }
  status=best.status

  if(top>1)
    return(list(status=status, res=getres(status),top.res=top.res,numtop=numtop,top.statuses=top.statuses,top.fun=top.fun  ))
  if(!do.return.status)
    return(getres(status))
  return(list(status=status, res=getres(status)))
}




regress.ic.search=function(data,response,covs,family="normal",
 offset="",strata="",use.glmmADMB=F,  use.glmmTMB=F,
 IC="BIC",threshold="flexible",
 check.est=F,check.se=F, do.return.status=F, start.status=NA, talkative=T,
 top=1)
{
  #print.srcref(sprintf("offset=%s",offset))

  numcov=length(covs)
  
  if(IC=="AICc")
    library(MuMIn)

  if(use.glmmADMB)
    library("glmmADMB")

  if(use.glmmTMB)
    library("glmmTMB")

  if(family=="clogit")
    library("survival")

  if(family=="ordinal" | family=="ordinal_cauchit")
    library("ordinal")
  
  # Assume random factors start with "(":
  isrand=substr(covs,1,1)=="("  

  numtop=0
  top.res=list()
  top.ic=rep(0,0)
  top.statuses=rep(0,0)
  if(top>1)
    do.return.status=TRUE
  
  # Create regression string:
  getmodel=function(status)
  {
   if(sum(status)==0)
    {
     ret=sprintf("%s~1",response)  
    }
   else
    {
       # only random effects?
       if((family=="ordinal" | family=="ordinal_cauchit") & 
          sum(status==1)==sum(isrand[status==1])) 
       {
         ret=sprintf("%s~1+",response)
       }
       else
       {
         ret=sprintf("%s~",response)
       }
       index=which(status==1)
       ret=sprintf("%s%s",ret,covs[index[1]])
       if(length(index)>1)
        for(i in index[2:length(index)])
         ret=sprintf("%s+%s", ret, covs[i]) 
    }
   if(offset!="")
    ret=sprintf("%s+%s",ret,offset)

   if(strata!="")
    ret=sprintf("%s+strata(%s)",ret,strata)

   #show(strata)
   #show(ret)
   return(ret)
  }
  
  
  # Make zero model:
  if(family=="normal")
  {
    if(offset!="")
      res0=lm(as.formula(sprintf("%s~1+%s",response,offset)),data=data)
    if(offset=="")
      res0=lm(as.formula(sprintf("%s~1",response)),data=data)
  }
  if(family=="clogit")
  {
    if(strata!="")
      res0=clogit(as.formula(sprintf("%s~%s+strata(%s)",response,covs[1],strata)),data=data)
    if(strata=="")
      res0=clogit(as.formula(sprintf("%s~%s",response,covs[1])),data=data)
  }
  if(family=="ordinal")
  {
    res0=clm(as.formula(sprintf("%s~1",response)),data=data,threshold=threshold)
  }
  if(family=="ordinal_cauchit")
  {
    res0=clm(as.formula(sprintf("%s~1",response)),data=data,threshold=threshold,link="cauchit")
  }
  if(family!="normal" & family!="clogit" & family!="ordinal" & family!="ordinal_cauchit")
  {
    if(offset!="")
      res0=glm(as.formula(sprintf("%s~1+%s",response,offset)),data=data,family=family)
    if(offset=="")
      res0=glm(as.formula(sprintf("%s~1",response)),data=data,family=family)
  }

  # Run regression:
  getres=function(status)
  {
    str=getmodel(status)
    if(family=="normal")
    {
      if(use.glmmADMB)
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
               family="gaussian") ), 
               abort=function() {return(res0)}))
	       
      if(use.glmmTMB)
        return(withRestarts(tryCatch(glmmTMB(as.formula(str), data=data, 
               family="gaussian") ), 
               abort=function() {return(res0)}))
	       
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(lmer(as.formula(str), data=data, REML=F) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(lm(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }
    
    if(family=="clogit")
    {
      return(withRestarts(tryCatch(clogit(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }

    if(family=="ordinal")
    {
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(clmm(as.formula(str), data=data,threshold=threshold) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(clm(as.formula(str), data=data,threshold=threshold) ), 
             abort=function() {return(res0)}))
    }
    
    if(family=="ordinal_cauchit")
    {
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(clmm(as.formula(str), link="cauchit", data=data,threshold=threshold) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(clm(as.formula(str), link="cauchit", data=data,threshold=threshold) ), 
             abort=function() {return(res0)}))
    }
    
    if(use.glmmADMB)
      {
        #print.srcref(sprintf("glmmadmb(%s, data=data, family=%s)", str,family))
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
             family=family) ), 
             abort=function() {return(res0)}))
      }
 	        
    if(use.glmmTMB)
      {
        #print.srcref(sprintf("glmmTMB(%s, data=data, family=%s)", str,family))
        return(withRestarts(tryCatch(glmmTMB(as.formula(str), data=data, 
             family=family) ), 
             abort=function() {return(res0)}))
      }
 	        
    if(sum(isrand[status!=0])>0)
    {
      print.srcref(sprintf("glmer(as.formula(%s), data=data, family=%s)",str,family))
      return(withRestarts(tryCatch(glmer(as.formula(str), data=data,family=family) ), 
             abort=function() {return(res0)}))
    }
    
    return(withRestarts(tryCatch(glm(as.formula(str), data=data,family=family) ), 
           abort=function() {return(res0)}))
  }
  
  # IC instead of stepwise hypothesis testing
  # allows for going up and down and sideways in the model search
  
  update.top=function(new.ic,new.status,new.res,top,numtop,top.ic,top.statuses,top.res)
  {
    if(is.na(new.ic) | is.infinite(new.ic))
      return(list(numtop=numtop,top.ic=top.ic,top.statuses=top.statuses,top.res=top.res))
   
    if(numtop>0)
    {
      #show(top.ic)
      if(new.ic<top.ic[numtop])
      {
        equalstatus=function(x) 
           sum(new.status==x)==numcov      
        
        if(sum(apply(top.statuses,1,equalstatus))==0) # not in already?
        {
          #show(sprintf("%f<%f",new.ic,top.ic[numtop]))
          top.ic=c(top.ic,new.ic)
          top.statuses=rbind(top.statuses,new.status)
          numtop=length(top.ic)
          #show(numtop)
          top.res[[numtop]]=new.res
          
	  # Order:
          o=order(top.ic)
          if(numtop>top)
          {
            o=o[1:top]
            numtop=top
            #show("numtop2")
            #show(numtop)
          }
          top.ic=top.ic[o]
          top.statuses=top.statuses[o,]
          top.res2=list()
          for(i in 1:top)
            top.res2[[i]]=top.res[[o[i]]]
          top.res=top.res2

        }  
      }    
    }
    if(numtop==0)
    {
      top.res[[numtop+1]]=new.res
      top.ic=c(new.ic)
      top.statuses=rbind(new.status)
      numtop=1
      #show(numtop)
    } 

   return(list(numtop=numtop,top.ic=top.ic,top.statuses=top.statuses,top.res=top.res))
  }


  # Run one step up/down/sideways (stop when the IC doesn't improve):
  step.ic=function(status,top,numtop,top.ic,top.statuses,top.res)
  {
  res0=getres(status)  
  
  ic0=best.ic=BIC(res0)
  if(IC=="AIC")
    ic0=best.ic=AIC(res0)
  if(IC=="AICc")
    ic0=best.ic=AICc(res0)
  if(IC!="BIC" & IC!="AIC" & IC!="AICc")
    ic0=best.ic=AIC(res0,k=as.numeric(IC))
  print.srcref(sprintf("Best model: %s",getmodel(status)))
  print.srcref(sprintf("IC=%f",ic0))
  best.status=status
  if(top>1)
  {
    utop=update.top(ic0,status,res0,top,numtop,top.ic,top.statuses,top.res)
    numtop=utop$numtop
    top.ic=utop$top.ic
    top.statuses=utop$top.statuses
    top.res=utop$top.res
  }

  # Add:
  if(talkative)
    show("Add")
  for(i in which(status==0))
  {
   status1=status
   status1[i]=1
   res1=getres(status1)
   ic1=BIC(res1)
   if(IC=="AIC")
    ic1=AIC(res1)
   if(IC=="AICc")
    ic1=AICc(res1)
   if(IC!="BIC" & IC!="AIC" & IC!="AICc")
    ic1=AIC(res1,k=as.numeric(IC))
   if(talkative)
     show(c(covs[i],i,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if(ic1<best.ic)
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.ic,top.statuses,top.res)
    numtop=utop$numtop
    top.ic=utop$top.ic
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  #Remove:
  if(talkative)
    show("Remove")
  for(i in which(status==1))
  {
   status1=status
   status1[i]=0
   res1=getres(status1)
   ic1=BIC(res1)
   if(IC=="AIC")
    ic1=AIC(res1)
   if(IC=="AICc")
    ic1=AICc(res1)
   if(IC!="BIC" & IC!="AIC" & IC!="AICc")
    ic1=AIC(res1,k=as.numeric(IC))
   if(talkative)
     show(c(covs[i],i,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if(ic1<best.ic | (ic1==best.ic & sum(status1)<sum(status)))
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.ic,top.statuses,top.res)
    numtop=utop$numtop
    top.ic=utop$top.ic
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  #Replace:
  if(talkative)
    show("Replace")
  for(i in which(status==1))
   for(j in which(status==0))
  {
   status1=status
   status1[i]=0
   status1[j]=1
   res1=getres(status1)
   ic1=BIC(res1)
   if(IC=="AIC")
    ic1=AIC(res1)
   if(IC=="AICc")
    ic1=AICc(res1)
   if(IC!="BIC" & IC!="AIC" & IC!="AICc")
    ic1=AIC(res1,k=as.numeric(IC))
   if(talkative)
     show(c(covs[i],i,"->",covs[j],j,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if(ic1<best.ic | (ic1==best.ic & sum(status1)<sum(status)))
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.ic,top.statuses,top.res)
    numtop=utop$numtop
    top.ic=utop$top.ic
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  return(list(best.status=best.status,numtop=numtop,top.ic=top.ic,top.statuses=top.statuses,top.res=top.res))
  }

  status=rep(0,numcov)
  if(sum(is.na(start.status))==0 & length(start.status)==numcov)
    status=start.status
  stepres=step.ic(status,top,numtop,top.ic,top.statuses,top.res)
  status.new=stepres$best.status
  if(top>0)
  {
    numtop=stepres$numtop
    top.ic=stepres$top.ic
    top.statuses=stepres$top.statuses
    top.res=stepres$top.res
  }
  while(sum(status.new!=status)>0)
  {
    status=status.new
    stepres=step.ic(status,top,numtop,top.ic,top.statuses,top.res)
    status.new=stepres$best.status
    if(top>0)
    {
      numtop=stepres$numtop
      top.ic=stepres$top.ic
      top.statuses=stepres$top.statuses
      top.res=stepres$top.res
    }
  }
  status=status.new
  
  if(top>1)
    return(list(status=status, res=getres(status),top.res=top.res,numtop=numtop,top.statuses=top.statuses,top.ic=top.ic,top.weights=exp(-(top.ic-min(top.ic))/2)/sum(exp(-(top.ic-min(top.ic))/2))  ))
  if(!do.return.status)
    return(getres(status))
  return(list(status=status, res=getres(status)))
}




regress.functioncriterion.search=function(data,response,covs,family="normal",fun=BIC,
 do.maximize=F,offset="",strata="",use.glmmADMB=F, threshold="flexible",
 check.est=F,check.se=F, do.return.status=F, start.status=NA, talkative=T,
 top=1)
{
  numcov=length(covs)
  
  if(family=="clogit")
    library("survival")

  if(family=="ordinal")
    library("ordinal")
  
  # Assume random factors start with "(":
  isrand=substr(covs,1,1)=="("  

  numtop=0
  top.res=list()
  top.fun=rep(0,0)
  top.statuses=rep(0,0)
  if(top>1)
    do.return.status=TRUE
  
  # Create regression string:
  getmodel=function(status)
  {
   if(sum(status)==0)
    {
     ret=sprintf("%s~1",response)  
    }
   else
    {
       # only random effects?
       if(family=="ordinal" & sum(status==1)==sum(isrand[status==1])) 
       {
         ret=sprintf("%s~1+",response)
       }
       else
       {
         ret=sprintf("%s~",response)
       }
       index=which(status==1)
       ret=sprintf("%s%s",ret,covs[index[1]])
       if(length(index)>1)
        for(i in index[2:length(index)])
         ret=sprintf("%s+%s", ret, covs[i]) 
    }
   if(offset!="")
    ret=sprintf("%s+%s",ret,offset)

   if(strata!="")
    ret=sprintf("%s+strata(%s)",ret,strata)
  
   return(ret)
  }
  
  
  # Make zero model:
  if(family=="normal")
  {
    if(offset!="")
      res0=lm(as.formula(sprintf("%s~1+%s",response,offset)),data=data)
    if(offset=="")
      res0=lm(as.formula(sprintf("%s~1",response)),data=data)
  }
  if(family=="clogit")
  {
    if(strata!="")
      res0=clogit(as.formula(sprintf("%s~%s+strata(%s)",response,covs[1],strata)),data=data)
    if(strata=="")
      res0=clogit(as.formula(sprintf("%s~%s",response,covs[1])),data=data)
  }
  if(family=="ordinal")
  {
    res0=clm(as.formula(sprintf("%s~1",response)),data=data,threshold=threshold)
  }
  if(family!="normal" & family!="clogit" & family!="ordinal")
  {
    if(offset!="")
      res0=glm(as.formula(sprintf("%s~1+%s",response,offset)),data=data,family=family)
    if(offset=="")
      res0=glm(as.formula(sprintf("%s~1",response)),data=data,family=family)
  }

  # Run regression:
  getres=function(status)
  {
    str=getmodel(status)
    if(family=="normal")
    {
      if(use.glmmADMB)
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
               family="gaussian") ), 
               abort=function() {return(res0)}))
	       
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(lmer(as.formula(str), data=data, REML=F) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(lm(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }
    
    if(family=="clogit")
    {
      return(withRestarts(tryCatch(clogit(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }

    if(family=="ordinal")
    {
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(clmm(as.formula(str), data=data,threshold=threshold) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(clm(as.formula(str), data=data,threshold=threshold) ), 
             abort=function() {return(res0)}))
    }
    
    if(use.glmmADMB)
      {
        #print.srcref(sprintf("glmmadmb(%s, data=data, family=%s)", str,family))
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
             family=family) ), 
             abort=function() {return(res0)}))
      }
 	      
    if(sum(isrand[status!=0])>0)
    {
      print.srcref(sprintf("glmer(as.formula(%s), data=data, family=%s)",str,family))
      return(withRestarts(tryCatch(glmer(as.formula(str), data=data,family=family) ), 
             abort=function() {return(res0)}))
    }
    
    return(withRestarts(tryCatch(glm(as.formula(str), data=data,family=family) ), 
           abort=function() {return(res0)}))
  }
  
  update.top=function(new.fun,new.status,new.res,top,numtop,top.fun,top.statuses,top.res)
  {
    if(is.na(new.fun) | is.infinite(new.fun))
      return(list(numtop=numtop,top.ic=top.ic,top.statuses=top.statuses,top.res=top.res))
   
    if(numtop>0)
    {
      #show(top.fun)
      if(new.fun>top.fun[numtop])
      {
        equalstatus=function(x) 
           sum(new.status==x)==numcov      
        
        if(sum(apply(top.statuses,1,equalstatus))==0) # not in already?
        {
          #show(sprintf("%f>%f",new.fun,top.fun[numtop]))
          top.fun=c(top.fun,new.fun)
          top.statuses=rbind(top.statuses,new.status)
          numtop=length(top.fun)
          #show(numtop)
          top.res[[numtop]]=new.res
          
	  # Order:
          o=order(top.fun,decreasing=TRUE)
          if(numtop>top)
          {
            o=o[1:top]
            numtop=top
            #show("numtop2")
            #show(numtop)
          }
          top.fun=top.fun[o]
          top.statuses=top.statuses[o,]
          top.res2=list()
          for(i in 1:top)
            top.res2[[i]]=top.res[[o[i]]]
          top.res=top.res2

        }  
      }    
    }
    if(numtop==0)
    {
      top.res[[numtop+1]]=new.res
      top.fun=c(new.fun)
      top.statuses=rbind(new.status)
      numtop=1
      #show(numtop)
    } 

   return(list(numtop=numtop,top.fun=top.fun,top.statuses=top.statuses,top.res=top.res))
  }

  # Run one step up/down/sideways (stop when the function doesn't improve):
  step.fun=function(status,top,numtop,top.fun,top.statuses,top.res)
  {
  res0=getres(status)  
  
  ic0=best.ic=fun(res0)

  print.srcref(sprintf("Best model: %s",getmodel(status)))
  print.srcref(sprintf("function=%f",ic0))
  best.status=status
  if(top>1)
  {
    utop=update.top(ic0,status,res0,top,numtop,top.fun,top.statuses,top.res)
    numtop=utop$numtop
    top.fun=utop$top.fun
    top.statuses=utop$top.statuses
    top.res=utop$top.res
  }
  
  # Add:
  if(talkative)
    show("Add")
  for(i in which(status==0))
  {
   status1=status
   status1[i]=1
   res1=getres(status1)
   ic1=fun(res1)
   if(talkative)
     show(c(covs[i],i,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if((!do.maximize & ic1<best.ic) | (do.maximize & ic1>best.ic))
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.fun,top.statuses,top.res)
    numtop=utop$numtop
    top.fun=utop$top.fun
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  #Remove:
  if(talkative)
    show("Remove")
  for(i in which(status==1))
  {
   status1=status
   status1[i]=0
   res1=getres(status1)
   ic1=fun(res1)
   if(talkative)
     show(c(covs[i],i,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if((!do.maximize & ic1<best.ic) | (do.maximize & ic1>best.ic) | 
        (ic1==best.ic & sum(status1)<sum(status)))
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.fun,top.statuses,top.res)
    numtop=utop$numtop
    top.fun=utop$top.fun
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  #Replace:
  if(talkative)
    show("Replace")
  for(i in which(status==1))
   for(j in which(status==0))
  {
   status1=status
   status1[i]=0
   status1[j]=1
   res1=getres(status1)
   ic1=fun(res1)
   if(talkative)
     show(c(covs[i],i,"->",covs[j],j,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if((!do.maximize & ic1<best.ic) | (do.maximize & ic1>best.ic) | 
        (ic1==best.ic & sum(status1)<sum(status)))
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.fun,top.statuses,top.res)
    numtop=utop$numtop
    top.fun=utop$top.fun
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  return(list(best.status=best.status,numtop=numtop,top.fun=top.fun,top.statuses=top.statuses,top.res=top.res))
  }

  status=rep(0,numcov)
  if(sum(is.na(start.status))==0 & length(start.status)==numcov)
    status=start.status
  stepres=step.fun(status,top,numtop,top.fun,top.statuses,top.res)
  status.new=stepres$best.status
  if(top>0)
  {
    numtop=stepres$numtop
    top.fun=stepres$top.fun
    top.statuses=stepres$top.statuses
    top.res=stepres$top.res
  }
  while(sum(status.new!=status)>0)
  {
    status=status.new
    stepres=step.fun(status,top,numtop,top.fun,top.statuses,top.res)
    status.new=stepres$best.status
    if(top>0)
    {
      numtop=stepres$numtop
      top.fun=stepres$top.fun
      top.statuses=stepres$top.statuses
      top.res=stepres$top.res
    }
  }
  status=status.new
  
  if(top>1)
    return(list(status=status, res=getres(status),top.res=top.res,numtop=numtop,top.statuses=top.statuses,top.fun=top.fun ))
  if(!do.return.status)
    return(getres(status))
  return(list(status=status, res=getres(status)))
}







regress.ic.search.2step=function(data,response,covs,family="normal",
 offset="",strata="",use.glmmADMB=F, IC="BIC",threshold="flexible",
 check.est=F,check.se=F, do.return.status=F, start.status=NA, talkative=T,
 top=1)
{
  numcov=length(covs)
  
  if(IC=="AICc")
    library(MuMIn)

  if(use.glmmADMB)
    library("glmmADMB")

  if(family=="clogit")
    library("survival")

  if(family=="ordinal")
    library("ordinal")
  
  # Assume random factors start with "(":
  isrand=substr(covs,1,1)=="("  

  
  numtop=0
  top.res=list()
  top.ic=rep(0,0)
  top.statuses=rep(0,0)
  if(top>1)
    do.return.status=TRUE
  
  # Create regression string:
  getmodel=function(status)
  {
   if(sum(status)==0)
    {
     ret=sprintf("%s~1",response)  
    }
   else
    {
       # only random effects?
       if(family=="ordinal" & sum(status==1)==sum(isrand[status==1])) 
       {
         ret=sprintf("%s~1+",response)
       }
       else
       {
         ret=sprintf("%s~",response)
       }
       index=which(status==1)
       ret=sprintf("%s%s",ret,covs[index[1]])
       if(length(index)>1)
        for(i in index[2:length(index)])
         ret=sprintf("%s+%s", ret, covs[i]) 
    }
   if(offset!="")
    ret=sprintf("%s+%s",ret,offset)

   if(strata!="")
    ret=sprintf("%s+strata(%s)",ret,strata)
  
   return(ret)
  }
  
  
  # Make zero model:
  if(family=="normal")
  {
    if(offset!="")
      res0=lm(as.formula(sprintf("%s~1+%s",response,offset)),data=data)
    if(offset=="")
      res0=lm(as.formula(sprintf("%s~1",response)),data=data)
  }
  if(family=="clogit")
  {
    if(strata!="")
      res0=clogit(as.formula(sprintf("%s~%s+strata(%s)",response,covs[1],strata)),data=data)
    if(strata=="")
      res0=clogit(as.formula(sprintf("%s~%s",response,covs[1])),data=data)
  }
  if(family=="ordinal")
  {
    res0=clm(as.formula(sprintf("%s~1",response)),data=data,threshold=threshold)
  }
  if(family!="normal" & family!="clogit" & family!="ordinal")
  {
    if(offset!="")
      res0=glm(as.formula(sprintf("%s~1+%s",response,offset)),data=data,family=family)
    if(offset=="")
      res0=glm(as.formula(sprintf("%s~1",response)),data=data,family=family)
  }

  # Run regression:
  getres=function(status)
  {
    str=getmodel(status)
    if(family=="normal")
    {
      if(use.glmmADMB)
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
               family="gaussian") ), 
               abort=function() {return(res0)}))
	       
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(lmer(as.formula(str), data=data, REML=F) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(lm(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }
    
    if(family=="clogit")
    {
      return(withRestarts(tryCatch(clogit(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }

    if(family=="ordinal")
    {
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(clmm(as.formula(str), data=data,threshold=threshold) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(clm(as.formula(str), data=data,threshold=threshold) ), 
             abort=function() {return(res0)}))
    }
    
    if(use.glmmADMB)
      {
        #print.srcref(sprintf("glmmadmb(%s, data=data, family=%s)", str,family))
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
             family=family) ), 
             abort=function() {return(res0)}))
      }
 	      
    if(sum(isrand[status!=0])>0)
    {
      #print.srcref(sprintf("glmer(as.formula(%s), data=data, family=%s)",str,family))
      return(withRestarts(tryCatch(glmer(as.formula(str), data=data,family=family) ), 
             abort=function() {return(res0)}))
    }
    
    return(withRestarts(tryCatch(glm(as.formula(str), data=data,family=family) ), 
           abort=function() {return(res0)}))
  }
  
  # IC instead of stepwise hypothesis testing
  # allows for going up and down and sideways in the model search
  
  update.top=function(new.ic,new.status,new.res,top,numtop,top.ic,top.statuses,top.res)
  {
    if(is.na(new.ic) | is.infinite(new.ic))
      return(list(numtop=numtop,top.ic=top.ic,top.statuses=top.statuses,top.res=top.res))
   
    if(numtop>0)
    {
      #show(top.ic)
      if(new.ic<top.ic[numtop])
      {
        equalstatus=function(x) 
           sum(new.status==x)==numcov      
        
        if(sum(apply(top.statuses,1,equalstatus))==0) # not in already?
        {
          #show(sprintf("%f<%f",new.ic,top.ic[numtop]))
          top.ic=c(top.ic,new.ic)
          top.statuses=rbind(top.statuses,new.status)
          numtop=length(top.ic)
          #show(numtop)
          top.res[[numtop]]=new.res
          
	  # Order:
          o=order(top.ic)
          if(numtop>top)
          {
            o=o[1:top]
            numtop=top
            #show("numtop2")
            #show(numtop)
          }
          top.ic=top.ic[o]
          top.statuses=top.statuses[o,]
          top.res2=list()
          for(i in 1:top)
            top.res2[[i]]=top.res[[o[i]]]
          top.res=top.res2

        }  
      }    
    }
    if(numtop==0)
    {
      top.res[[numtop+1]]=new.res
      top.ic=c(new.ic)
      top.statuses=rbind(new.status)
      numtop=1
      #show(numtop)
    } 

   return(list(numtop=numtop,top.ic=top.ic,top.statuses=top.statuses,top.res=top.res))
  }

  # Run one step up/down/sideways (stop when the IC doesn't improve):
  step.ic=function(status,top,numtop,top.ic,top.statuses,top.res)
  {
  res0=getres(status)  

  ic0=best.ic=BIC(res0)
  if(IC=="AIC")
    ic0=best.ic=AIC(res0)
  if(IC=="AICc")
    ic0=best.ic=AICc(res0)
  if(IC!="BIC" & IC!="AIC" & IC!="AICc")
    ic0=best.ic=AIC(res0,k=as.numeric(IC))
  if(talkative)
  {
    print.srcref(sprintf("Best model: %s",getmodel(status)))
    print.srcref(sprintf("IC=%f",ic0))
  }
  best.status=status
  if(top>1)
  {
    utop=update.top(ic0,status,res0,top,numtop,top.ic,top.statuses,top.res)
    numtop=utop$numtop
    top.ic=utop$top.ic
    top.statuses=utop$top.statuses
    top.res=utop$top.res
  }

  # Add one:
  #show("Add one")
  for(i in which(status==0))
  {
   status1=status
   status1[i]=1
   res1=getres(status1)
   ic1=BIC(res1)
   if(IC=="AIC")
    ic1=AIC(res1)
   if(IC=="AICc")
    ic1=AICc(res1)
   if(IC!="BIC" & IC!="AIC" & IC!="AICc")
    ic1=AIC(res1,k=as.numeric(IC))
   #show(c(covs[i],i,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if(ic1<best.ic)
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.ic,top.statuses,top.res)
    numtop=utop$numtop
    top.ic=utop$top.ic
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  # Add two:
  #show("Add two")
  for(i in which(status==0))
   for(j in which(status==0))
    if(i<j)
  {
   status1=status
   status1[i]=1
   status1[j]=1
   res1=getres(status1)
   ic1=BIC(res1)
   if(IC=="AIC")
    ic1=AIC(res1)
   if(IC=="AICc")
    ic1=AICc(res1)
   if(IC!="BIC" & IC!="AIC" & IC!="AICc")
    ic1=AIC(res1,k=as.numeric(IC))
   #show(c(covs[i],i,covs[j],j,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if(ic1<best.ic)
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.ic,top.statuses,top.res)
    numtop=utop$numtop
    top.ic=utop$top.ic
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  #Remove one:
  #show("Remove one")
  for(i in which(status==1))
  {
   status1=status
   status1[i]=0
   res1=getres(status1)
   ic1=BIC(res1)
   if(IC=="AIC")
    ic1=AIC(res1)
   if(IC=="AICc")
    ic1=AICc(res1)
   if(IC!="BIC" & IC!="AIC" & IC!="AICc")
    ic1=AIC(res1,k=as.numeric(IC))
   #show(c(covs[i],i,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if(ic1<best.ic | (ic1==best.ic & sum(status1)<sum(status)))
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.ic,top.statuses,top.res)
    numtop=utop$numtop
    top.ic=utop$top.ic
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  #Remove two:
  #show("Remove")
  for(i in which(status==1))
   for(j in which(status==1))
    if(i<j)
  {
   status1=status
   status1[i]=0
   status1[j]=0
   res1=getres(status1)
   ic1=BIC(res1)
   if(IC=="AIC")
    ic1=AIC(res1)
   if(IC=="AICc")
    ic1=AICc(res1)
   if(IC!="BIC" & IC!="AIC" & IC!="AICc")
    ic1=AIC(res1,k=as.numeric(IC))
   #show(c(covs[i],i,covs[j],j,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if(ic1<best.ic | (ic1==best.ic & sum(status1)<sum(status)))
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.ic,top.statuses,top.res)
    numtop=utop$numtop
    top.ic=utop$top.ic
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  #Replace one:
  #show("Replace")
  for(i in which(status==1))
   for(j in which(status==0))
  {
   status1=status
   status1[i]=0
   status1[j]=1
   res1=getres(status1)
   ic1=BIC(res1)
   if(IC=="AIC")
    ic1=AIC(res1)
   if(IC=="AICc")
    ic1=AICc(res1)
   if(IC!="BIC" & IC!="AIC" & IC!="AICc")
    ic1=AIC(res1,k=as.numeric(IC))
   #show(c(covs[i],i,"->",covs[j],j,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if(ic1<best.ic | (ic1==best.ic & sum(status1)<sum(status)))
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.ic,top.statuses,top.res)
    numtop=utop$numtop
    top.ic=utop$top.ic
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  #Replace one, remove one:
  #show("Replace one, remove one")
  for(i in which(status==1))
   for(j in which(status==0))
    for(k in which(status==1))
     if(i>k)
  {
   status1=status
   status1[i]=0
   status1[j]=1
   status1[k]=0
   res1=getres(status1)
   ic1=BIC(res1)
   if(IC=="AIC")
    ic1=AIC(res1)
   if(IC=="AICc")
    ic1=AICc(res1)
   if(IC!="BIC" & IC!="AIC" & IC!="AICc")
    ic1=AIC(res1,k=as.numeric(IC))
   #show(c(covs[i],i,"->",covs[j],j,"rm",covs[k],k,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if(ic1<best.ic | (ic1==best.ic & sum(status1)<sum(status)))
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.ic,top.statuses,top.res)
    numtop=utop$numtop
    top.ic=utop$top.ic
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  #Replace one, add one:
  #show("Replace one, add one")
  for(i in which(status==1))
   for(j in which(status==0))
    for(k in which(status==0))
     if(j>k)
  {
   status1=status
   status1[i]=0
   status1[j]=1
   status1[k]=1
   res1=getres(status1)
   ic1=BIC(res1)
   if(IC=="AIC")
    ic1=AIC(res1)
   if(IC=="AICc")
    ic1=AICc(res1)
   if(IC!="BIC" & IC!="AIC" & IC!="AICc")
    ic1=AIC(res1,k=as.numeric(IC))
   #show(c(covs[i],i,"->",covs[j],j,"add",covs[k],k,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if(ic1<best.ic | (ic1==best.ic & sum(status1)<sum(status)))
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.ic,top.statuses,top.res)
    numtop=utop$numtop
    top.ic=utop$top.ic
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  #Replace two:
  #show("Replace two")
  for(i in which(status==1))
   for(j in which(status==0))
    for(k in which(status==1))
     for(l in which(status==0))
       if(i>k & j>l)
  {
   status1=status
   status1[i]=0
   status1[j]=1
   status1[k]=0
   status1[l]=1
   res1=getres(status1)
   ic1=BIC(res1)
   if(IC=="AIC")
    ic1=AIC(res1)
   if(IC=="AICc")
    ic1=AICc(res1)
   if(IC!="BIC" & IC!="AIC" & IC!="AICc")
    ic1=AIC(res1,k=as.numeric(IC))
   #show(c(covs[i],i,covs[k],k,"->",covs[j],j,covs[l],l,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if(ic1<best.ic | (ic1==best.ic & sum(status1)<sum(status)))
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.ic,top.statuses,top.res)
    numtop=utop$numtop
    top.ic=utop$top.ic
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  return(list(best.status=best.status,numtop=numtop,top.ic=top.ic,top.statuses=top.statuses,top.res=top.res))
  }

  status=rep(0,numcov)
  if(sum(is.na(start.status))==0 & length(start.status)==numcov)
    status=start.status
  stepres=step.ic(status,top,numtop,top.ic,top.statuses,top.res)
  status.new=stepres$best.status
  if(top>0)
  {
    numtop=stepres$numtop
    top.ic=stepres$top.ic
    top.statuses=stepres$top.statuses
    top.res=stepres$top.res
  }
  while(sum(status.new!=status)>0)
  {
    status=status.new
    stepres=step.ic(status,top,numtop,top.ic,top.statuses,top.res)
    status.new=stepres$best.status
    if(top>0)
    {
      numtop=stepres$numtop
      top.ic=stepres$top.ic
      top.statuses=stepres$top.statuses
      top.res=stepres$top.res
    }
  }
  status=status.new
  
  if(top>1)
    return(list(status=status, res=getres(status),top.res=top.res,numtop=numtop,top.statuses=top.statuses,top.ic=top.ic,top.weights=exp(-(top.ic-min(top.ic))/2)/sum(exp(-(top.ic-min(top.ic))/2))  ))
  if(!do.return.status)
    return(getres(status))
  return(list(status=status, res=getres(status)))
}




regress.functioncriterion.search.2step=function(data,response,covs,family="normal",
 fun=BIC, do.maximize=F, offset="",strata="",use.glmmADMB=F, threshold="flexible",
 check.est=F,check.se=F, do.return.status=F, start.status=NA, talkative=T, top=1)
{
  numcov=length(covs)
  
  if(use.glmmADMB)
    library("glmmADMB")

  if(family=="clogit")
    library("survival")

  if(family=="ordinal")
    library("ordinal")
  
  # Assume random factors start with "(":
  isrand=substr(covs,1,1)=="("  

  numtop=0
  top.res=list()
  top.fun=rep(0,0)
  top.statuses=rep(0,0)
  if(top>1)
    do.return.status=TRUE
  

  # Create regression string:
  getmodel=function(status)
  {
   if(sum(status)==0)
    {
     ret=sprintf("%s~1",response)  
    }
   else
    {
       # only random effects?
       if(family=="ordinal" & sum(status==1)==sum(isrand[status==1])) 
       {
         ret=sprintf("%s~1+",response)
       }
       else
       {
         ret=sprintf("%s~",response)
       }
       index=which(status==1)
       ret=sprintf("%s%s",ret,covs[index[1]])
       if(length(index)>1)
        for(i in index[2:length(index)])
         ret=sprintf("%s+%s", ret, covs[i]) 
    }
   if(offset!="")
    ret=sprintf("%s+%s",ret,offset)

   if(strata!="")
    ret=sprintf("%s+strata(%s)",ret,strata)
  
   return(ret)
  }
  
  
  # Make zero model:
  if(family=="normal")
  {
    if(offset!="")
      res0=lm(as.formula(sprintf("%s~1+%s",response,offset)),data=data)
    if(offset=="")
      res0=lm(as.formula(sprintf("%s~1",response)),data=data)
  }
  if(family=="clogit")
  {
    if(strata!="")
      res0=clogit(as.formula(sprintf("%s~%s+strata(%s)",response,covs[1],strata)),data=data)
    if(strata=="")
      res0=clogit(as.formula(sprintf("%s~%s",response,covs[1])),data=data)
  }
  if(family=="ordinal")
  {
    res0=clm(as.formula(sprintf("%s~1",response)),data=data,threshold=threshold)
  }
  if(family!="normal" & family!="clogit" & family!="ordinal")
  {
    if(offset!="")
      res0=glm(as.formula(sprintf("%s~1+%s",response,offset)),data=data,family=family)
    if(offset=="")
      res0=glm(as.formula(sprintf("%s~1",response)),data=data,family=family)
  }

  # Run regression:
  getres=function(status)
  {
    str=getmodel(status)
    if(family=="normal")
    {
      if(use.glmmADMB)
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
               family="gaussian") ), 
               abort=function() {return(res0)}))
	       
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(lmer(as.formula(str), data=data, REML=F) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(lm(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }
    
    if(family=="clogit")
    {
      return(withRestarts(tryCatch(clogit(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }

    if(family=="ordinal")
    {
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(clmm(as.formula(str), data=data,threshold=threshold) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(clm(as.formula(str), data=data,threshold=threshold) ), 
             abort=function() {return(res0)}))
    }
    
    if(use.glmmADMB)
      {
        #print.srcref(sprintf("glmmadmb(%s, data=data, family=%s)", str,family))
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
             family=family) ), 
             abort=function() {return(res0)}))
      }
 	      
    if(sum(isrand[status!=0])>0)
    {
      #print.srcref(sprintf("glmer(as.formula(%s), data=data, family=%s)",str,family))
      return(withRestarts(tryCatch(glmer(as.formula(str), data=data,family=family) ), 
             abort=function() {return(res0)}))
    }
    
    return(withRestarts(tryCatch(glm(as.formula(str), data=data,family=family) ), 
           abort=function() {return(res0)}))
  }
  
  update.top=function(new.fun,new.status,new.res,top,numtop,top.fun,top.statuses,top.res)
  {
    if(is.na(new.fun) | is.infinite(new.fun))
      return(list(numtop=numtop,top.ic=top.ic,top.statuses=top.statuses,top.res=top.res))
   
    if(numtop>0)
    {
      #show(top.fun)
      if(new.fun>top.fun[numtop])
      {
        equalstatus=function(x) 
           sum(new.status==x)==numcov      
        
        if(sum(apply(top.statuses,1,equalstatus))==0) # not in already?
        {
          #show(sprintf("%f>%f",new.fun,top.fun[numtop]))
          top.fun=c(top.fun,new.fun)
          top.statuses=rbind(top.statuses,new.status)
          numtop=length(top.fun)
          #show(numtop)
          top.res[[numtop]]=new.res
          
	  # Order:
          o=order(top.fun,decreasing=TRUE)
          if(numtop>top)
          {
            o=o[1:top]
            numtop=top
            #show("numtop2")
            #show(numtop)
          }
          top.fun=top.fun[o]
          top.statuses=top.statuses[o,]
          top.res2=list()
          for(i in 1:top)
            top.res2[[i]]=top.res[[o[i]]]
          top.res=top.res2

        }  
      }    
    }
    if(numtop==0)
    {
      top.res[[numtop+1]]=new.res
      top.fun=c(new.fun)
      top.statuses=rbind(new.status)
      numtop=1
      #show(numtop)
    } 

   return(list(numtop=numtop,top.fun=top.fun,top.statuses=top.statuses,top.res=top.res))
  }

  # Run one step up/down/sideways (stop when the IC doesn't improve):
  step.fun=function(status)
  {
  res0=getres(status)  

  ic0=best.ic=fun(res0)
  if(talkative)
  {
    print.srcref(sprintf("Best model: %s",getmodel(status)))
    print.srcref(sprintf("function=%f",ic0))
  }
  best.status=status
  if(top>1)
  {
    utop=update.top(ic0,status,res0,top,numtop,top.fun,top.statuses,top.res)
    numtop=utop$numtop
    top.fun=utop$top.fun
    top.statuses=utop$top.statuses
    top.res=utop$top.res
  }

  # Add one:
  #show("Add one")
  for(i in which(status==0))
  {
   status1=status
   status1[i]=1
   res1=getres(status1)
   ic1=fun(res1)
   #show(c(covs[i],i,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if((!do.maximize & ic1<best.ic) | (do.maximize & ic1>best.ic))
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.fun,top.statuses,top.res)
    numtop=utop$numtop
    top.fun=utop$top.fun
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  # Add two:
  #show("Add two")
  for(i in which(status==0))
   for(j in which(status==0))
    if(i<j)
  {
   status1=status
   status1[i]=1
   status1[j]=1
   res1=getres(status1)
   ic1=fun(res1)
   #show(c(covs[i],i,covs[j],j,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if((!do.maximize & ic1<best.ic) | (do.maximize & ic1>best.ic))
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.fun,top.statuses,top.res)
    numtop=utop$numtop
    top.fun=utop$top.fun
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  #Remove one:
  #show("Remove one")
  for(i in which(status==1))
  {
   status1=status
   status1[i]=0
   res1=getres(status1)
   ic1=fun(res1)
   #show(c(covs[i],i,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if((!do.maximize & ic1<best.ic) | (do.maximize & ic1>best.ic) |
        (ic1==best.ic & sum(status1)<sum(status)))
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.fun,top.statuses,top.res)
    numtop=utop$numtop
    top.fun=utop$top.fun
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  #Remove two:
  #show("Remove")
  for(i in which(status==1))
   for(j in which(status==1))
    if(i<j)
  {
   status1=status
   status1[i]=0
   status1[j]=0
   res1=getres(status1)
   ic1=fun(res1)
   #show(c(covs[i],i,covs[j],j,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if((!do.maximize & ic1<best.ic) | (do.maximize & ic1>best.ic) |
        (ic1==best.ic & sum(status1)<sum(status)))
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.fun,top.statuses,top.res)
    numtop=utop$numtop
    top.fun=utop$top.fun
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  #Replace one:
  #show("Replace")
  for(i in which(status==1))
   for(j in which(status==0))
  {
   status1=status
   status1[i]=0
   status1[j]=1
   res1=getres(status1)
   ic1=fun(res1)
   #show(c(covs[i],i,"->",covs[j],j,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if((!do.maximize & ic1<best.ic) | (do.maximize & ic1>best.ic) |
        (ic1==best.ic & sum(status1)<sum(status)))
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.fun,top.statuses,top.res)
    numtop=utop$numtop
    top.fun=utop$top.fun
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  #Replace one, remove one:
  #show("Replace one, remove one")
  for(i in which(status==1))
   for(j in which(status==0))
    for(k in which(status==1))
     if(i>k)
  {
   status1=status
   status1[i]=0
   status1[j]=1
   status1[k]=0
   res1=getres(status1)
   ic1=fun(res1)
   #show(c(covs[i],i,"->",covs[j],j,"rm",covs[k],k,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if((!do.maximize & ic1<best.ic) | (do.maximize & ic1>best.ic) |
        (ic1==best.ic & sum(status1)<sum(status)))
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.fun,top.statuses,top.res)
    numtop=utop$numtop
    top.fun=utop$top.fun
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  #Replace one, add one:
  #show("Replace one, add one")
  for(i in which(status==1))
   for(j in which(status==0))
    for(k in which(status==0))
     if(j>k)
  {
   status1=status
   status1[i]=0
   status1[j]=1
   status1[k]=1
   res1=getres(status1)
   ic1=fun(res1)
   #show(c(covs[i],i,"->",covs[j],j,"add",covs[k],k,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if((!do.maximize & ic1<best.ic) | (do.maximize & ic1>best.ic) |
        (ic1==best.ic & sum(status1)<sum(status)))
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.fun,top.statuses,top.res)
    numtop=utop$numtop
    top.fun=utop$top.fun
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  #Replace two:
  #show("Replace two")
  for(i in which(status==1))
   for(j in which(status==0))
    for(k in which(status==1))
     for(l in which(status==0))
       if(i>k & j>l)
  {
   status1=status
   status1[i]=0
   status1[j]=1
   status1[k]=0
   status1[l]=1
   res1=getres(status1)
   ic1=fun(res1)
   #show(c(covs[i],i,covs[k],k,"->",covs[j],j,covs[l],l,ic1))
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
   if(check.est)
    if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
   if(!is.na(ic1))
   {
    if(is.na(ic0))
    {
     best.ic=ic1
     best.status=status1
    }
    if(!is.na(ic0))
     if((!do.maximize & ic1<best.ic) | (do.maximize & ic1>best.ic) |
        (ic1==best.ic & sum(status1)<sum(status)))
    {
     best.ic=ic1
     best.status=status1
    }
   }
   if(top>1)
   {
    utop=update.top(ic1,status1,res1,top,numtop,top.fun,top.statuses,top.res)
    numtop=utop$numtop
    top.fun=utop$top.fun
    top.statuses=utop$top.statuses
    top.res=utop$top.res
   }
  }

  return(list(best.status=best.status,numtop=numtop,top.fun=top.fun,top.statuses=top.statuses,top.res=top.res))
  }

  status=rep(0,numcov)
  if(sum(is.na(start.status))==0 & length(start.status)==numcov)
    status=start.status
  stepres=step.fun(status,top,numtop,top.fun,top.statuses,top.res)
  status.new=stepres$best.status
  if(top>0)
  {
    numtop=stepres$numtop
    top.fun=stepres$top.fun
    top.statuses=stepres$top.statuses
    top.res=stepres$top.res
  }
  while(sum(status.new!=status)>0)
  {
    status=status.new
    stepres=step.fun(status,top,numtop,top.fun,top.statuses,top.res)
    status.new=stepres$best.status
    if(top>0)
    {
      numtop=stepres$numtop
      top.fun=stepres$top.fun
      top.statuses=stepres$top.statuses
      top.res=stepres$top.res
    }
  }
  status=status.new
  
  if(top>1)
    return(list(status=status, res=getres(status),top.res=top.res,numtop=numtop,top.statuses=top.statuses,top.fun=top.fun ))
  if(!do.return.status)
    return(getres(status))
  return(list(status=status, res=getres(status)))
}








regress.random.ic.search=function(data,response,covs,num.tries.per.iteration=10000,
 num.tries.before.update=100,
 family="normal", offset="",strata="",use.glmmADMB=F, IC="BIC",threshold="flexible",
 check.est=F,check.se=F, do.return.status=F, start.status=NA,talkative=F, top=1)
{
  numcov=length(covs)
  
  link=""

  if(IC=="AICc")
    library(MuMIn)

  if(use.glmmADMB)
    library("glmmADMB")

  if(family=="clogit")
    library("survival")

  if(family=="ordinal" | family=="ordinal_cauchit")
    library("ordinal")
  
  # Assume random factors start with "(":
  isrand=substr(covs,1,1)=="("  

  numtop=0
  top.res=list()
  top.ic=rep(0,0)
  top.statuses=rep(0,0)
  if(top>1)
    do.return.status=TRUE
  
  
  # Create regression string:
  getmodel=function(status)
  {
   if(sum(status)==0)
    {
     ret=sprintf("%s~1",response)  
    }
   else
    {
       # only random effects?
       if((family=="ordinal" | family=="ordinal_cauchit") & 
          sum(status==1)==sum(isrand[status==1])) 
       {
         ret=sprintf("%s~1+",response)
       }
       else
       {
         ret=sprintf("%s~",response)
       }
       index=which(status==1)
       ret=sprintf("%s%s",ret,covs[index[1]])
       if(length(index)>1)
        for(i in index[2:length(index)])
         ret=sprintf("%s+%s", ret, covs[i]) 
    }
   if(offset!="")
    ret=sprintf("%s+%s",ret,offset)

   if(strata!="")
    ret=sprintf("%s+strata(%s)",ret,strata)
  
   return(ret)
  }
  
  
  # Make zero model:
  if(family=="normal")
  {
    if(offset!="")
      res0=lm(as.formula(sprintf("%s~1+%s",response,offset)),data=data)
    if(offset=="")
      res0=lm(as.formula(sprintf("%s~1",response)),data=data)
  }
  if(family=="clogit")
  {
    if(strata!="")
      res0=clogit(as.formula(sprintf("%s~%s+strata(%s)",response,covs[1],strata)),data=data)
    if(strata=="")
      res0=clogit(as.formula(sprintf("%s~%s",response,covs[1])),data=data)
  }
  if(family=="ordinal")
  {
    res0=clm(as.formula(sprintf("%s~1",response)),data=data,threshold=threshold)
  }
  if(family=="ordinal_cauchit")
  {
    res0=clm(as.formula(sprintf("%s~1",response)),data=data,
             threshold=threshold,link="cauchit")
  }
  if(family!="normal" & family!="clogit" & family!="ordinal" & family!="ordinal_cauchit")
  {
    if(offset!="")
      res0=glm(as.formula(sprintf("%s~1+%s",response,offset)),data=data,family=family)
    if(offset=="")
      res0=glm(as.formula(sprintf("%s~1",response)),data=data,family=family)
  }

  # Run regression:
  getres=function(status)
  {
    str=getmodel(status)
    if(family=="normal")
    {
      if(use.glmmADMB)
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
               family="gaussian") ), 
               abort=function() {return(res0)}))
	       
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(lmer(as.formula(str), data=data, REML=F) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(lm(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }
    
    if(family=="clogit")
    {
      return(withRestarts(tryCatch(clogit(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }

    if(family=="ordinal")
    {
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(clmm(as.formula(str), data=data,threshold=threshold) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(clm(as.formula(str), data=data,threshold=threshold) ), 
             abort=function() {return(res0)}))
    }
    
    if(family=="ordinal_cauchit")
    {
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(clmm(as.formula(str), data=data,threshold=threshold, link="cauchit") ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(clm(as.formula(str), data=data,threshold=threshold, link="cauchit") ), 
             abort=function() {return(res0)}))
    }
    
    if(use.glmmADMB)
      {
        #print.srcref(sprintf("glmmadmb(%s, data=data, family=%s)", str,family))
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
             family=family) ), 
             abort=function() {return(res0)}))
      }
 	      
    if(sum(isrand[status!=0])>0)
    {
      print.srcref(sprintf("glmer(as.formula(%s), data=data, family=%s)",str,family))
      return(withRestarts(tryCatch(glmer(as.formula(str), data=data,family=family) ), 
             abort=function() {return(res0)}))
    }
    
    return(withRestarts(tryCatch(glm(as.formula(str), data=data,family=family) ), 
           abort=function() {return(res0)}))
  }
  
  # IC instead of stepwise hypothesis testing
  # allows for going up and down and sideways in the model search
  
  update.top=function(new.ic,new.status,new.res,top,numtop,top.ic,top.statuses,top.res)
  {
    if(is.na(new.ic) | is.infinite(new.ic))
      return(list(numtop=numtop,top.ic=top.ic,top.statuses=top.statuses,top.res=top.res))
   
    if(numtop>0)
    {
      #show(top.ic)
      if(new.ic<top.ic[numtop])
      {
        equalstatus=function(x) 
           sum(new.status==x)==numcov      
        
        if(sum(apply(top.statuses,1,equalstatus))==0) # not in already?
        {
          #show(sprintf("%f<%f",new.ic,top.ic[numtop]))
          top.ic=c(top.ic,new.ic)
          top.statuses=rbind(top.statuses,new.status)
          numtop=length(top.ic)
          #show(numtop)
          top.res[[numtop]]=new.res
          
	  # Order:
          o=order(top.ic)
          if(numtop>top)
          {
            o=o[1:top]
            numtop=top
            #show("numtop2")
            #show(numtop)
          }
          top.ic=top.ic[o]
          top.statuses=top.statuses[o,]
          top.res2=list()
          for(i in 1:top)
            top.res2[[i]]=top.res[[o[i]]]
          top.res=top.res2

        }  
      }    
    }
    if(numtop==0)
    {
      top.res[[numtop+1]]=new.res
      top.ic=c(new.ic)
      top.statuses=rbind(new.status)
      numtop=1
      #show(numtop)
    } 

   return(list(numtop=numtop,top.ic=top.ic,top.statuses=top.statuses,top.res=top.res))
  }

  # Run one step up/down/sideways (stop when the IC doesn't improve):
  step.random.ic=function(status,num.tries,top,numtop,top.ic,top.statuses,top.res)
  {
  status0=status
  res0=getres(status0)  

  ic0=best.ic=BIC(res0)
  if(IC=="AIC")
    ic0=best.ic=AIC(res0)
  if(IC=="AICc")
    ic0=best.ic=AICc(res0)
  if(IC!="BIC" & IC!="AIC" & IC!="AICc")
    ic0=best.ic=AIC(res0,k=as.numeric(IC))

  if(talkative)
  {
   print.srcref("")
   print.srcref("New iteration")
   print.srcref("Best model from previous iteration")
   print.srcref(sprintf("%s: IC=%7.2f", getmodel(status0), best.ic))
   print.srcref("")
  }
  best.status=status0
  if(top>1)
  {
    utop=update.top(ic0,status,res0,top,numtop,top.ic,top.statuses,top.res)
    numtop=utop$numtop
    top.ic=utop$top.ic
    top.statuses=utop$top.statuses
    top.res=utop$top.res
  }
  
  for(i in 1:num.tries)
   {
    if(i%%num.tries.before.update==0)
      status0=best.status
    
    status1=status0
    l=length(status0)
    num.added=sum(status0)  

    num.add=sample(0:4,1,replace=F,prob=c(1/8,1/2+1/32,1/4,1/16,1/32))
    num.rem=sample(0:4,1,replace=F,prob=c(3/4+1/64, 1/8, 1/16,1/32,1/64))
    while(num.add+num.rem==0 | num.rem>num.added | num.add>(l-num.added))
    {
     num.add=sample(0:4,1,replace=F,prob=c(1/2+1/32,1/4,1/8,1/16,1/32))
     num.rem=sample(0:4,1,replace=F,prob=c(1/2+1/32,1/4,1/8,1/16,1/32))
    }
    
    index.add=sample(which(status0==0),num.add,replace=F)
    index.rem=sample(which(status0==1),num.rem,replace=F)
    index.change=c(index.add,index.rem)
    
    status1[index.change]=1-status1[index.change]
    res1=getres(status1)
    ic1=BIC(res1)
    if(IC=="AIC")
     ic1=AIC(res1)
    if(IC=="AICc")
     ic1=AICc(res1)
    if(IC!="BIC" & IC!="AIC" & IC!="AICc")
     ic1=AIC(res1,k=as.numeric(IC))
    if(check.est)
     if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
    if(check.est)
     if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
    if(!is.na(ic1))
    {
     if(is.na(ic0))
     {
      best.ic=ic1
      best.status=status1
      if(talkative)
       print.srcref(sprintf("i=%5d %s: IC=%7.2f",i, getmodel(status1), ic1))
     }
     if(!is.na(ic0))
      if(ic1<best.ic)
      {
       best.ic=ic1
       best.status=status1
       if(talkative)
        print.srcref(sprintf("i=%5d %s: IC=%7.2f",i, getmodel(status1), ic1))
      }
    }  
    if(top>1)
    {
     utop=update.top(ic1,status1,res1,top,numtop,top.ic,top.statuses,top.res)
     numtop=utop$numtop
     top.ic=utop$top.ic
     top.statuses=utop$top.statuses
     top.res=utop$top.res
    }
   }

  return(list(best.status=best.status,numtop=numtop,top.ic=top.ic,top.statuses=top.statuses,top.res=top.res))
  }

  status=rep(0,numcov)
  if(sum(is.na(start.status))==0 & length(start.status)==numcov)
    status=start.status
  stepres=step.random.ic(status,num.tries.per.iteration,top,numtop,top.ic,top.statuses,top.res)
  status.new=stepres$best.status
  if(top>0)
  {
    numtop=stepres$numtop
    top.ic=stepres$top.ic
    top.statuses=stepres$top.statuses
    top.res=stepres$top.res
  }
  
  while(sum(status.new!=status)>0)
  {
    status=status.new
    stepres=step.random.ic(status,num.tries.per.iteration,top,numtop,top.ic,top.statuses,top.res)
    status.new=stepres$best.status
    if(top>0)
    {
      numtop=stepres$numtop
      top.ic=stepres$top.ic
      top.statuses=stepres$top.statuses
      top.res=stepres$top.res
    }
  }
  status=status.new
  
  if(top>1)
    return(list(status=status, res=getres(status),top.res=top.res,numtop=numtop,top.statuses=top.statuses,top.ic=top.ic,top.weights=exp(-(top.ic-min(top.ic))/2)/sum(exp(-(top.ic-min(top.ic))/2))  ))
  if(!do.return.status)
    return(getres(status))
  return(list(status=status, res=getres(status)))
}




regress.random.functioncriterion.search=function(data,response,covs,
 fun=BIC, do.maximize=F,
 num.tries.per.iteration=10000, num.tries.before.update=100,
 family="normal", offset="",strata="",use.glmmADMB=F, IC="BIC",threshold="flexible",
 check.est=F,check.se=F, do.return.status=F, start.status=NA,talkative=F, top=1)
{
  numcov=length(covs)
  
  if(use.glmmADMB)
    library("glmmADMB")

  if(family=="clogit")
    library("survival")

  if(family=="ordinal")
    library("ordinal")
  
  # Assume random factors start with "(":
  isrand=substr(covs,1,1)=="("  
  
  numtop=0
  top.res=list()
  top.fun=rep(0,0)
  top.statuses=rep(0,0)
  if(top>1)
    do.return.status=TRUE
  
  # Create regression string:
  getmodel=function(status)
  {
   if(sum(status)==0)
    {
     ret=sprintf("%s~1",response)  
    }
   else
    {
       # only random effects?
       if(family=="ordinal" & sum(status==1)==sum(isrand[status==1])) 
       {
         ret=sprintf("%s~1+",response)
       }
       else
       {
         ret=sprintf("%s~",response)
       }
       index=which(status==1)
       ret=sprintf("%s%s",ret,covs[index[1]])
       if(length(index)>1)
        for(i in index[2:length(index)])
         ret=sprintf("%s+%s", ret, covs[i]) 
    }
   if(offset!="")
    ret=sprintf("%s+%s",ret,offset)

   if(strata!="")
    ret=sprintf("%s+strata(%s)",ret,strata)
  
   return(ret)
  }
  

  # Make zero model:
  if(family=="normal")
  {
    if(offset!="")
      res0=lm(as.formula(sprintf("%s~1+%s",response,offset)),data=data)
    if(offset=="")
      res0=lm(as.formula(sprintf("%s~1",response)),data=data)
  }
  if(family=="clogit")
  {
    if(strata!="")
      res0=clogit(as.formula(sprintf("%s~%s+strata(%s)",response,covs[1],strata)),data=data)
    if(strata=="")
      res0=clogit(as.formula(sprintf("%s~%s",response,covs[1])),data=data)
  }
  if(family=="ordinal")
  {
    res0=clm(as.formula(sprintf("%s~1",response)),data=data,threshold=threshold)
  }
  if(family!="normal" & family!="clogit" & family!="ordinal")
  {
    if(offset!="")
      res0=glm(as.formula(sprintf("%s~1+%s",response,offset)),data=data,family=family)
    if(offset=="")
      res0=glm(as.formula(sprintf("%s~1",response)),data=data,family=family)
  }

  # Run regression:
  getres=function(status)
  {
    str=getmodel(status)
    if(family=="normal")
    {
      if(use.glmmADMB)
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
               family="gaussian") ), 
               abort=function() {return(res0)}))
	       
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(lmer(as.formula(str), data=data, REML=F) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(lm(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }
    
    if(family=="clogit")
    {
      return(withRestarts(tryCatch(clogit(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }

    if(family=="ordinal")
    {
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(clmm(as.formula(str), data=data,threshold=threshold) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(clm(as.formula(str), data=data,threshold=threshold) ), 
             abort=function() {return(res0)}))
    }
    
    if(use.glmmADMB)
      {
        #print.srcref(sprintf("glmmadmb(%s, data=data, family=%s)", str,family))
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
             family=family) ), 
             abort=function() {return(res0)}))
      }
 	      
    if(sum(isrand[status!=0])>0)
    {
      if(talkative)
       print.srcref(sprintf("glmer(as.formula(%s), data=data, family=%s)",str,family))
      return(withRestarts(tryCatch(glmer(as.formula(str), data=data,family=family) ), 
             abort=function() {return(res0)}))
    }
    
    return(withRestarts(tryCatch(glm(as.formula(str), data=data,family=family) ), 
           abort=function() {return(res0)}))
  }
  
  update.top=function(new.fun,new.status,new.res,top,numtop,top.fun,top.statuses,top.res)
  {
    if(is.na(new.fun) | is.infinite(new.fun))
      return(list(numtop=numtop,top.ic=top.ic,top.statuses=top.statuses,top.res=top.res))
   
    if(numtop>0)
    {
      #show(top.fun)
      if(new.fun>top.fun[numtop])
      {
        equalstatus=function(x) 
           sum(new.status==x)==numcov      
        
        if(sum(apply(top.statuses,1,equalstatus))==0) # not in already?
        {
          #show(sprintf("%f>%f",new.fun,top.fun[numtop]))
          top.fun=c(top.fun,new.fun)
          top.statuses=rbind(top.statuses,new.status)
          numtop=length(top.fun)
          #show(numtop)
          top.res[[numtop]]=new.res
          
	  # Order:
          o=order(top.fun,decreasing=TRUE)
          if(numtop>top)
          {
            o=o[1:top]
            numtop=top
            #show("numtop2")
            #show(numtop)
          }
          top.fun=top.fun[o]
          top.statuses=top.statuses[o,]
          top.res2=list()
          for(i in 1:top)
            top.res2[[i]]=top.res[[o[i]]]
          top.res=top.res2

        }  
      }    
    }
    if(numtop==0)
    {
      top.res[[numtop+1]]=new.res
      top.fun=c(new.fun)
      top.statuses=rbind(new.status)
      numtop=1
      #show(numtop)
    } 

   return(list(numtop=numtop,top.fun=top.fun,top.statuses=top.statuses,top.res=top.res))
  }

  # Run one step up/down/sideways (stop when the IC doesn't improve):
  step.random.fun=function(status,num.tries,top,numtop,top.fun,top.statuses,top.res)
  {
  status0=status
  res0=getres(status0)  

  ic0=best.ic=fun(res0)

  if(talkative)
  {
   print.srcref("")
   print.srcref("New iteration")
   print.srcref("Best model from previous iteration")
   print.srcref(sprintf("%s: function=%7.2f", getmodel(status0), best.ic))
   print.srcref("")
  }
  best.status=status0
  if(top>1)
  {
    utop=update.top(ic0,status,res0,top,numtop,top.fun,top.statuses,top.res)
    numtop=utop$numtop
    top.fun=utop$top.fun
    top.statuses=utop$top.statuses
    top.res=utop$top.res
  }
  
  for(i in 1:num.tries)
   {
    if(i%%num.tries.before.update==0)
      status0=best.status
    
    status1=status0
    l=length(status0)
    num.added=sum(status0)  

    num.add=sample(0:4,1,replace=F,prob=c(1/8,1/2+1/32,1/4,1/16,1/32))
    num.rem=sample(0:4,1,replace=F,prob=c(3/4+1/64, 1/8, 1/16,1/32,1/64))
    while(num.add+num.rem==0 | num.rem>num.added | num.add>(l-num.added))
    {
     num.add=sample(0:4,1,replace=F,prob=c(1/2+1/32,1/4,1/8,1/16,1/32))
     num.rem=sample(0:4,1,replace=F,prob=c(1/2+1/32,1/4,1/8,1/16,1/32))
    }
    
    index.add=sample(which(status0==0),num.add,replace=F)
    index.rem=sample(which(status0==1),num.rem,replace=F)
    index.change=c(index.add,index.rem)
    
    status1[index.change]=1-status1[index.change]
    res1=getres(status1)
    ic1=fun(res1)
    if(check.est)
     if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
    if(check.est)
     if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
    if(!is.na(ic1))
    {
     if(is.na(ic0))
     {
      best.ic=ic1
      best.status=status1
      if(talkative)
       print.srcref(sprintf("i=%5d %s: IC=%7.2f",i, getmodel(status1), ic1))
     }
     if(!is.na(ic0))
      if((!do.maximize & ic1<best.ic) | (do.maximize & ic1>best.ic) | 
         (ic1==best.ic & sum(status1)<sum(status0)))
      {
       best.ic=ic1
       best.status=status1
       if(talkative)
        print.srcref(sprintf("i=%5d %s: IC=%7.2f",i, getmodel(status1), ic1))
      }
    } 
    if(top>1)
    {
     utop=update.top(ic1,status1,res1,top,numtop,top.fun,top.statuses,top.res)
     numtop=utop$numtop
     top.fun=utop$top.fun
     top.statuses=utop$top.statuses
     top.res=utop$top.res
    } 
   }

  return(list(best.status=best.status,numtop=numtop,top.fun=top.fun,top.statuses=top.statuses,top.res=top.res))
  }

  status=rep(0,numcov)
  if(sum(is.na(start.status))==0 & length(start.status)==numcov)
    status=start.status
  stepres=step.random.fun(status,num.tries.per.iteration,
    top,numtop,top.fun,top.statuses,top.res)
  status.new=stepres$best.status
  if(top>0)
  {
    numtop=stepres$numtop
    top.fun=stepres$top.fun
    top.statuses=stepres$top.statuses
    top.res=stepres$top.res
  }
  while(sum(status.new!=status)>0)
  {
    status=status.new
    stepres=step.random.fun(status,num.tries.per.iteration,
      top,numtop,top.fun,top.statuses,top.res)
    status.new=stepres$best.status
    if(top>0)
    {
      numtop=stepres$numtop
      top.fun=stepres$top.fun
      top.statuses=stepres$top.statuses
      top.res=stepres$top.res
    }
  }
  status=status.new
  
  if(top>1)
    return(list(status=status, res=getres(status),top.res=top.res,numtop=numtop,top.statuses=top.statuses,top.fun=top.fun ))
  if(!do.return.status)
    return(getres(status))
  return(list(status=status, res=getres(status)))
}







regress.mcmc.ic.search=function(data,response,covs,num.mcmc=10000,
 family="normal", offset="",strata="",use.glmmADMB=F,
 IC="BIC",threshold="flexible", check.est=F, check.se=F,
 do.return.status=F, start.status=NA, talkative=F, top=1)
{
  numcov=length(covs)
  
  if(IC=="AICc")
    library(MuMIn)

  if(use.glmmADMB)
    library("glmmADMB")

  if(family=="clogit")
    library("survival")

  if(family=="ordinal")
    library("ordinal")
  
  # Assume random factors start with "(":
  isrand=substr(covs,1,1)=="("  

  numtop=0
  top.res=list()
  top.ic=rep(0,0)
  top.statuses=rep(0,0)
  if(top>1)
    do.return.status=TRUE
  
  # Create regression string:
  getmodel=function(status)
  {
   if(sum(status)==0)
    {
     ret=sprintf("%s~1",response)  
    }
   else
    {
       # only random effects?
       if(family=="ordinal" & sum(status==1)==sum(isrand[status==1])) 
       {
         ret=sprintf("%s~1+",response)
       }
       else
       {
         ret=sprintf("%s~",response)
       }
       index=which(status==1)
       ret=sprintf("%s%s",ret,covs[index[1]])
       if(length(index)>1)
        for(i in index[2:length(index)])
         ret=sprintf("%s+%s", ret, covs[i]) 
    }
   if(offset!="")
    ret=sprintf("%s+%s",ret,offset)

   if(strata!="")
    ret=sprintf("%s+strata(%s)",ret,strata)
  
   return(ret)
  }
  
  
  # Make zero model:
  if(family=="normal")
  {
    if(offset!="")
      res0=lm(as.formula(sprintf("%s~1+%s",response,offset)),data=data)
    if(offset=="")
      res0=lm(as.formula(sprintf("%s~1",response)),data=data)
  }
  if(family=="clogit")
  {
    if(strata!="")
      res0=clogit(as.formula(sprintf("%s~%s+strata(%s)",response,covs[1],strata)),data=data)
    if(strata=="")
      res0=clogit(as.formula(sprintf("%s~%s",response,covs[1])),data=data)
  }
  if(family=="ordinal")
  {
    res0=clm(as.formula(sprintf("%s~1",response)),data=data,threshold=threshold)
  }
  if(family!="normal" & family!="clogit" & family!="ordinal")
  {
    if(offset!="")
      res0=glm(as.formula(sprintf("%s~1+%s",response,offset)),data=data,family=family)
    if(offset=="")
      res0=glm(as.formula(sprintf("%s~1",response)),data=data,family=family)
  }

  # Run regression:
  getres=function(status)
  {
    str=getmodel(status)
    if(family=="normal")
    {
      if(use.glmmADMB)
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
               family="gaussian") ), 
               abort=function() {return(res0)}))
	       
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(lmer(as.formula(str), data=data, REML=F) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(lm(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }
    
    if(family=="clogit")
    {
      return(withRestarts(tryCatch(clogit(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }

    if(family=="ordinal")
    {
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(clmm(as.formula(str), data=data,threshold=threshold) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(clm(as.formula(str), data=data,threshold=threshold) ), 
             abort=function() {return(res0)}))
    }
    
    if(use.glmmADMB)
      {
        #print.srcref(sprintf("glmmadmb(%s, data=data, family=%s)", str,family))
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
             family=family) ), 
             abort=function() {return(res0)}))
      }
 	      
    if(sum(isrand[status!=0])>0)
    {
      if(talkative)
       print.srcref(sprintf("glmer(as.formula(%s), data=data, family=%s)",str,family))
      return(withRestarts(tryCatch(glmer(as.formula(str), data=data,family=family) ), 
             abort=function() {return(res0)}))
    }
    
    return(withRestarts(tryCatch(glm(as.formula(str), data=data,family=family) ), 
           abort=function() {return(res0)}))
  }
  
  # IC instead of stepwise hypothesis testing
  # allows for going up and down and sideways in the model search
  
  update.top=function(new.ic,new.status,new.res,top,numtop,top.ic,top.statuses,top.res)
  {
    if(is.na(new.ic) | is.infinite(new.ic))
      return(list(numtop=numtop,top.ic=top.ic,top.statuses=top.statuses,top.res=top.res))
   
    if(numtop>0)
    {
      #show(top.ic)
      if(new.ic<top.ic[numtop])
      {
        equalstatus=function(x) 
           sum(new.status==x)==numcov      
        
        if(sum(apply(top.statuses,1,equalstatus))==0) # not in already?
        {
          #show(sprintf("%f<%f",new.ic,top.ic[numtop]))
          top.ic=c(top.ic,new.ic)
          top.statuses=rbind(top.statuses,new.status)
          numtop=length(top.ic)
          #show(numtop)
          top.res[[numtop]]=new.res
          
	  # Order:
          o=order(top.ic)
          if(numtop>top)
          {
            o=o[1:top]
            numtop=top
            #show("numtop2")
            #show(numtop)
          }
          top.ic=top.ic[o]
          top.statuses=top.statuses[o,]
          top.res2=list()
          for(i in 1:top)
            top.res2[[i]]=top.res[[o[i]]]
          top.res=top.res2

        }  
      }    
    }
    if(numtop==0)
    {
      top.res[[numtop+1]]=new.res
      top.ic=c(new.ic)
      top.statuses=rbind(new.status)
      numtop=1
      #show(numtop)
    } 

   return(list(numtop=numtop,top.ic=top.ic,top.statuses=top.statuses,top.res=top.res))
  }


  # Run 'num.tries' mcmc step:
  step.mcmc.ic=function(status,num.tries,top,numtop,top.ic,top.statuses,top.res)
  {
  status0=status
  res0=getres(status0)  
  best.status=status0

  ic0=best.ic=BIC(res0)
  if(IC=="AIC")
    ic0=best.ic=AIC(res0)
  if(IC=="AICc")
    ic0=best.ic=AICc(res0)
  if(IC!="BIC" & IC!="AIC" & IC!="AICc")
    ic0=best.ic=AIC(res0,k=as.numeric(IC))
  if(top>1)
  {
    utop=update.top(ic0,status,res0,top,numtop,top.ic,top.statuses,top.res)
    numtop=utop$numtop
    top.ic=utop$top.ic
    top.statuses=utop$top.statuses
    top.res=utop$top.res
  }
  
  for(i in 1:num.tries)
   {
    status1=status0
    l=length(status0)
    num.added=sum(status0) 
    
    num.change=sample(1:5,1,replace=F,prob=c(1/2+1/32,1/4,1/8,1/16,1/32))
    index.change=sample(1:length(status),num.change,replace=F)
    
    status1[index.change]=1-status1[index.change]
    res1=getres(status1)
    ic1=BIC(res1)
    if(IC=="AIC")
     ic1=AIC(res1)
    if(IC=="AICc")
     ic1=AICc(res1)
    if(IC!="BIC" & IC!="AIC" & IC!="AICc")
     ic1=AIC(res1,k=as.numeric(IC))
    if(check.est)
     if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
    if(check.est)
     if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
    if(!is.na(ic1))
    {
     if(is.na(ic0))
     {
      ic0=best.ic=ic1
      status0=best.status=status1
      if(talkative)
       print.srcref(sprintf("i=%5d %s: IC=%7.2f",i, getmodel(status1), ic1))
     }
     if(!is.na(ic0))
     {
      if(log(runif(1)) < -(ic1-ic0)/2)
      {
       ic0=ic1
       status0=status1
       if(ic1<best.ic)
       {
        best.ic=ic1
        best.status=status1
        if(talkative)
         print.srcref(sprintf("i=%5d %s: IC=%7.2f",i, getmodel(status1), ic1))
       }
      }
     }
    } 
    if(top>1)
    {
     utop=update.top(ic1,status1,res1,top,numtop,top.ic,top.statuses,top.res)
     numtop=utop$numtop
     top.ic=utop$top.ic
     top.statuses=utop$top.statuses
     top.res=utop$top.res
    }
   }

  return(list(best.status=best.status,numtop=numtop,top.ic=top.ic,top.statuses=top.statuses,top.res=top.res))
  }

  status=rep(0,numcov)
  if(sum(is.na(start.status))==0 & length(start.status)==numcov)
    status=start.status
  stepres=step.mcmc.ic(status,num.mcmc,top,numtop,top.ic,top.statuses,top.res)
  status.new=stepres$best.status
  if(top>0)
  {
    numtop=stepres$numtop
    top.ic=stepres$top.ic
    top.statuses=stepres$top.statuses
    top.res=stepres$top.res
  }
  status=status.new
  
  if(top>1)
    return(list(status=status, res=getres(status),top.res=top.res,numtop=numtop,top.statuses=top.statuses,top.ic=top.ic,top.weights=exp(-(top.ic-min(top.ic))/2)/sum(exp(-(top.ic-min(top.ic))/2))  ))
  if(!do.return.status)
    return(getres(status))
  return(list(status=status, res=getres(status)))
}



regress.mcmc.functioncriterion.search=function(data,response,covs,num.mcmc=10000,
 fun=BIC, do.maximize=F,
 family="normal", offset="",strata="",use.glmmADMB=F, threshold="flexible",
 check.est=F, check.se=F, do.return.status=F, start.status=NA, talkative=F, top=1)
{
  numcov=length(covs)
  
  if(use.glmmADMB)
    library("glmmADMB")

  if(family=="clogit")
    library("survival")

  if(family=="ordinal")
    library("ordinal")
  
  # Assume random factors start with "(":
  isrand=substr(covs,1,1)=="("  

  numtop=0
  top.res=list()
  top.fun=rep(0,0)
  top.statuses=rep(0,0)
  if(top>1)
    do.return.status=TRUE
  
  
  # Create regression string:
  getmodel=function(status)
  {
   if(sum(status)==0)
    {
     ret=sprintf("%s~1",response)  
    }
   else
    {
       # only random effects?
       if(family=="ordinal" & sum(status==1)==sum(isrand[status==1])) 
       {
         ret=sprintf("%s~1+",response)
       }
       else
       {
         ret=sprintf("%s~",response)
       }
       index=which(status==1)
       ret=sprintf("%s%s",ret,covs[index[1]])
       if(length(index)>1)
        for(i in index[2:length(index)])
         ret=sprintf("%s+%s", ret, covs[i]) 
    }
   if(offset!="")
    ret=sprintf("%s+%s",ret,offset)

   if(strata!="")
    ret=sprintf("%s+strata(%s)",ret,strata)
  
   return(ret)
  }
  

  # Make zero model:
  if(family=="normal")
  {
    if(offset!="")
      res0=lm(as.formula(sprintf("%s~1+%s",response,offset)),data=data)
    if(offset=="")
      res0=lm(as.formula(sprintf("%s~1",response)),data=data)
  }
  if(family=="clogit")
  {
    if(strata!="")
      res0=clogit(as.formula(sprintf("%s~%s+strata(%s)",response,covs[1],strata)),data=data)
    if(strata=="")
      res0=clogit(as.formula(sprintf("%s~%s",response,covs[1])),data=data)
  }
  if(family=="ordinal")
  {
    res0=clm(as.formula(sprintf("%s~1",response)),data=data,threshold=threshold)
  }
  if(family!="normal" & family!="clogit" & family!="ordinal")
  {
    if(offset!="")
      res0=glm(as.formula(sprintf("%s~1+%s",response,offset)),data=data,family=family)
    if(offset=="")
      res0=glm(as.formula(sprintf("%s~1",response)),data=data,family=family)
  }

  # Run regression:
  getres=function(status)
  {
    str=getmodel(status)
    if(family=="normal")
    {
      if(use.glmmADMB)
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
               family="gaussian") ), 
               abort=function() {return(res0)}))
	       
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(lmer(as.formula(str), data=data, REML=F) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(lm(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }
    
    if(family=="clogit")
    {
      return(withRestarts(tryCatch(clogit(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }

    if(family=="ordinal")
    {
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(clmm(as.formula(str), data=data,threshold=threshold) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(clm(as.formula(str), data=data,threshold=threshold) ), 
             abort=function() {return(res0)}))
    }
    
    if(use.glmmADMB)
      {
        #print.srcref(sprintf("glmmadmb(%s, data=data, family=%s)", str,family))
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
             family=family) ), 
             abort=function() {return(res0)}))
      }
 	      
    if(sum(isrand[status!=0])>0)
    {
      if(talkative)
       print.srcref(sprintf("glmer(as.formula(%s), data=data, family=%s)",str,family))
      return(withRestarts(tryCatch(glmer(as.formula(str), data=data,family=family) ), 
             abort=function() {return(res0)}))
    }
    
    return(withRestarts(tryCatch(glm(as.formula(str), data=data,family=family) ), 
           abort=function() {return(res0)}))
  }
  
  update.top=function(new.fun,new.status,new.res,top,numtop,top.fun,top.statuses,top.res)
  {
    if(is.na(new.fun) | is.infinite(new.fun))
      return(list(numtop=numtop,top.ic=top.ic,top.statuses=top.statuses,top.res=top.res))
   
    if(numtop>0)
    {
      #show(top.fun)
      if(new.fun>top.fun[numtop])
      {
        equalstatus=function(x) 
           sum(new.status==x)==numcov      
        
        if(sum(apply(top.statuses,1,equalstatus))==0) # not in already?
        {
          #show(sprintf("%f>%f",new.fun,top.fun[numtop]))
          top.fun=c(top.fun,new.fun)
          top.statuses=rbind(top.statuses,new.status)
          numtop=length(top.fun)
          #show(numtop)
          top.res[[numtop]]=new.res
          
	  # Order:
          o=order(top.fun,decreasing=TRUE)
          if(numtop>top)
          {
            o=o[1:top]
            numtop=top
            #show("numtop2")
            #show(numtop)
          }
          top.fun=top.fun[o]
          top.statuses=top.statuses[o,]
          top.res2=list()
          for(i in 1:top)
            top.res2[[i]]=top.res[[o[i]]]
          top.res=top.res2

        }  
      }    
    }
    if(numtop==0)
    {
      top.res[[numtop+1]]=new.res
      top.fun=c(new.fun)
      top.statuses=rbind(new.status)
      numtop=1
      #show(numtop)
    } 

   return(list(numtop=numtop,top.fun=top.fun,top.statuses=top.statuses,top.res=top.res))
  }

  # Run 'num.tries' mcmc step:
  step.mcmc.fun=function(status,num.tries,top,numtop,top.fun,top.statuses,top.res)
  {
  status0=status
  res0=getres(status0)  
  best.status=status0

  ic0=best.ic=fun(res0)
  if(top>1)
  {
    utop=update.top(ic0,status,res0,top,numtop,top.fun,top.statuses,top.res)
    numtop=utop$numtop
    top.fun=utop$top.fun
    top.statuses=utop$top.statuses
    top.res=utop$top.res
  }
  
  for(i in 1:num.tries)
   {
    status1=status0
    l=length(status0)
    num.added=sum(status0) 
    
    num.change=sample(1:5,1,replace=F,prob=c(1/2+1/32,1/4,1/8,1/16,1/32))
    index.change=sample(1:length(status),num.change,replace=F)
    
    status1[index.change]=1-status1[index.change]
    res1=getres(status1)
    ic1=fun(res1)
    if(check.est)
     if(sum(is.na(summary(res1)$coef[,1]))>0)
      ic1=NA
    if(check.est)
     if(sum(is.na(summary(res1)$coef[,2]))>0)
      ic1=NA
    if(!is.na(ic1))
    {
     if(is.na(ic0))
     {
      ic0=best.ic=ic1
      status0=best.status=status1
      if(talkative)
       print.srcref(sprintf("i=%5d %s: IC=%7.2f",i, getmodel(status1), ic1))
     }
     if(!is.na(ic0))
     {
      if(log(runif(1)) < -(ic1-ic0)/2)
      {
       ic0=ic1
       status0=status1
       if((!do.maximize & ic1<best.ic) | (do.maximize & ic1>best.ic) |
          (ic1==best.ic & sum(status1)<sum(status0)))
       {
        best.ic=ic1
        best.status=status1
        if(talkative)
         print.srcref(sprintf("i=%5d %s: IC=%7.2f",i, getmodel(status1), ic1))
       }
      }
     }
    } 
    if(top>1)
    {
     utop=update.top(ic1,status1,res1,top,numtop,top.fun,top.statuses,top.res)
     numtop=utop$numtop
     top.fun=utop$top.fun
     top.statuses=utop$top.statuses
     top.res=utop$top.res
    }
   }

  return(list(best.status=best.status,numtop=numtop,top.fun=top.fun,top.statuses=top.statuses,top.res=top.res))
  }

  status=rep(0,numcov)
  if(sum(is.na(start.status))==0 & length(start.status)==numcov)
    status=start.status
  stepres=step.mcmc.ic(status,num.mcmc,top,numtop,top.fun,top.statuses,top.res)
  status.new=stepres$best.status
  if(top>0)
  {
    numtop=stepres$numtop
    top.fun=stepres$top.fun
    top.statuses=stepres$top.statuses
    top.res=stepres$top.res
  }
  status=status.new
  
  if(top>1)
    return(list(status=status, res=getres(status),top.res=top.res,numtop=numtop,top.statuses=top.statuses,top.fun=top.fun ))
  if(!do.return.status)
    return(getres(status))
  return(list(status=status, res=getres(status)))
}





regress.hypo.search=function(data,response,covs,family="normal",offset="",strata="",use.glmmADMB=F,threshold="flexible", bonferroni=FALSE)
{
  numcov=length(covs)
  
  if(use.glmmADMB)
    library("glmmADMB")

  if(family=="clogit")
    library("survival")

  # Assume random factors start with "(":
  isrand=substr(covs,1,1)=="("  

  
  # Create regression string:
  getmodel=function(status)
  {
   if(sum(status)==0)
    {
     ret=sprintf("%s~1",response)  
    }
   else
    {
       # only random effects?
       if(family=="ordinal" & sum(status==1)==sum(isrand[status==1])) 
       {
         ret=sprintf("%s~1+",response)
       }
       else
       {
         ret=sprintf("%s~",response)
       }
       index=which(status==1)
       ret=sprintf("%s%s",ret,covs[index[1]])
       if(length(index)>1)
        for(i in index[2:length(index)])
         ret=sprintf("%s+%s", ret, covs[i]) 
    }
   if(offset!="")
    ret=sprintf("%s+%s",ret,offset)

   if(strata!="")
    ret=sprintf("%s+strata(%s)",ret,strata)

   show(ret)
   return(ret)
  }
  

  # Make zero model:
  if(family=="normal")
  {
    if(offset!="")
      res0=lm(as.formula(sprintf("%s~1+%s",response,offset)),data=data)
    if(offset=="")
      res0=lm(as.formula(sprintf("%s~1",response)),data=data)
  }
  if(family=="clogit")
  {
    if(strata!="")
      res0=clogit(as.formula(sprintf("%s~%s+strata(%s)",response,covs[1],strata)),data=data)
    if(strata=="")
      res0=clogit(as.formula(sprintf("%s~%s",response,covs[1])),data=data)
  }
  if(family=="ordinal")
  {
    res0=clm(as.formula(sprintf("%s~1",response)),data=data,threshold=threshold)
  }
  if(family!="normal" & family!="clogit" & family!="ordinal")
  {
    if(offset!="")
      res0=glm(as.formula(sprintf("%s~1+%s",response,offset)),data=data,family=family)
    if(offset=="")
      res0=glm(as.formula(sprintf("%s~1",response)),data=data,family=family)
  }

  # Run regression:
  getres=function(status)
  {
    str=getmodel(status)
    if(family=="normal")
    {
      if(use.glmmADMB)
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
               family="gaussian") ), 
               abort=function() {return(res0)}))
	       
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(lmer(as.formula(str), data=data, REML=F) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(lm(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }
    
    if(family=="clogit")
    {
      return(withRestarts(tryCatch(clogit(as.formula(str), data=data) ), 
             abort=function() {return(res0)}))
    }

    if(family=="ordinal")
    {
      if(sum(isrand[status!=0])>0)
        return(withRestarts(tryCatch(clmm(as.formula(str), data=data,threshold=threshold) ), 
               abort=function() {return(res0)}))

      return(withRestarts(tryCatch(clm(as.formula(str), data=data,threshold=threshold) ), 
             abort=function() {return(res0)}))
    }
    
    if(use.glmmADMB)
      {
        #print.srcref(sprintf("glmmadmb(%s, data=data, family=%s)", str,family))
        return(withRestarts(tryCatch(glmmadmb(as.formula(str), data=data, 
             family=family) ), 
             abort=function() {return(res0)}))
      }
 	      
    if(sum(isrand[status!=0])>0)
    {
      #print.srcref(sprintf("glmer(as.formula(%s), data=data, family=%s)",str,family))
      return(withRestarts(tryCatch(glmer(as.formula(str), data=data,family=family) ), 
             abort=function() {return(res0)}))
    }
    
    return(withRestarts(tryCatch(glm(as.formula(str), data=data,family=family) ), 
           abort=function() {return(res0)}))
  }
  
  # Run one step up/down/sideways (stop when the p-values doesn't improve):
  step.hypo=function(status)
  {
  res0=getres(status)  

  show("")
  show("new iteration")
  show(getmodel(status))
  show("")
  best.status=status
  best.p=1
  
  # Add:
  show("Add")
  for(i in which(status==0))
  {
   status1=status
   status1[i]=1
   res1=getres(status1)
   p.value=anova(res0,res1,test="Chisq")$Pr[2]

   show(c(covs[i],i,p.value))
   if(!is.na(p.value))
   {
   
     if((bonferroni==FALSE & p.value<best.p & p.value<0.05) |
        (bonferroni==TRUE & p.value<best.p & p.value<0.05/sum(status==0)))
    {
     best.p=p.value
     best.status=status1
    }
   }
  }

  #Remove:
  show("Remove")
  for(i in which(status==1))
  {
   status1=status
   status1[i]=0
   res1=getres(status1)
   p.value=anova(res1,res0,test="Chisq")$Pr[2]
   best.p=0

   show(c(covs[i],i,p.value))
   if(!is.na(p.value))
   {
    if(p.value>0.05 & p.value>best.p)
    {
     best.p=p.value
     best.status=status1
    }
   }
  }

  return(best.status)
  }

  status=rep(0,numcov)
  status.new=step.hypo(status)
  while(sum(status.new!=status)>0)
  {
    status=status.new
    status.new=step.hypo(status)
  }
  status=status.new

  show(getmodel(status))
  return(getres(status))
}

