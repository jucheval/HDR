require(abind)

#####################################################################################
#######  Auxiliary functions
#####################################################################################
Delta_function<-function(input_intensity,L1_norm,refractory_period)
{
  return((1+input_intensity*refractory_period-L1_norm)^{2} + 4*input_intensity*L1_norm*refractory_period)
}

theoretical_steady_activity<-function(input_intensity,L1_norm,refractory_period)
{
  if (refractory_period<0){warning("refractory period should be non negative")}
  #if (L1_norm<0){warning("connectivity should be non negative")}
  #if (input_intensity<0){warning("input intensity should be non negative")}
  if (L1_norm==0){return(1/(refractory_period+1/input_intensity))}
  else
  {
    if (refractory_period==0)
    {
      if (L1_norm>=1){warning("if refractory period is null then connectivity must be less than one")}
      else{return(input_intensity/(1-L1_norm))}
    }
    else
    {
      return(1/(refractory_period) - (1+L1_norm+input_intensity*refractory_period - sqrt(Delta_function(input_intensity,L1_norm,refractory_period)))/(2*L1_norm*refractory_period))
    }
  }
}

theoretical_sensitivity<-function(input_intensity,L1_norm,refractory_period)
{
  if (refractory_period<=0){warning("refractory period should be positive")}
  #if (L1_norm<=0){warning("connectivity should be positive")}
  if (input_intensity<=0){warning("input intensity should be positive")}
  return(-1/(2*L1_norm) + (1+input_intensity*refractory_period+L1_norm)/(2*L1_norm*sqrt(Delta_function(input_intensity,L1_norm,refractory_period))))
}

theoretical_sensitivity_derivative<-function(input_intensity,L1_norm,refractory_period)
{
  if (refractory_period<=0){warning("refractory period should be positive")}
  if (L1_norm<=0){warning("connectivity should be positive")}
  if (input_intensity<=0){warning("input intensity should be positive")}
  return((2*Delta_function(input_intensity,L1_norm,refractory_period)^{3/2} - 2*(1+input_intensity*refractory_period)*Delta_function(input_intensity,L1_norm,refractory_period) - 2*L1_norm*((input_intensity*refractory_period+L1_norm)^{2}-1))/(4*L1_norm^{2}*Delta_function(input_intensity,L1_norm,refractory_period)^{3/2}))
}

simul_Hawkes_expo_with_age<-function(nb_neuron,refractory_period,intensity_function,L1_norm,time_constant,age_init,X_init = 0,nb_points=10000)
{
  # Fonction qui simule un processus de Hawkes (dependant de l'age) multivarie de dimension n=nb_neuron 
  # Les processus N^1,...,N^n ont pour intensite lambda^i_t = intensity_function(S^i_t,X(t))
  # où S^i_t est l'age du neurone i au temps t
  # et X(t) = int_0^t L1_norm*time_constant*exp(-time_constant*(t-t')) [n^{-1} \sum N^i(dt')]
  # On s'arrete une fois qu'on a atteint nb_points points (pour le processus dominant).
  t=0
  X=X_init
  age=age_init
  lambda_majorant=intensity_function(Inf,X)
  points=rep(0,nb_points)
  type=rep(0,nb_points)
  intensity=rep(0,nb_points)
  for (k in 1:nb_points)
  {
    spike_possible=which(age>refractory_period-1/intensity_function(Inf,0)) # On considere que seuls les neurones qui vont etre hors periode refractaire dans un petit moment peuvent spiker
    if (length(spike_possible)==0)
    {
      waiting_time=min(refractory_period-age)
      t=t+waiting_time
      age=age+waiting_time
      X=exp(-time_constant*waiting_time)*X
      lambda_majorant=intensity_function(Inf,X)
    }
    else
    {
      dominant_ISI=rexp(1,lambda_majorant)/length(spike_possible) # ISI pour un des neurones qui peut spiker
      if (dominant_ISI==0){break} # Pour sortir d'un cas de blow-up
      t=t+dominant_ISI
      age=age+dominant_ISI
      new_X=exp(-time_constant*dominant_ISI)*X
      process_number=sample(spike_possible,1)
      if (runif(1)<intensity_function(age[process_number],new_X)/lambda_majorant)    # Si ca passe, alors le processus numero process_number a un spike. Sinon, il ne spike pas a cause de la periode refractaire ou a cause de la majoration.
      {
        points[k]=t
        type[k]=process_number
        age[process_number]=0
        X=new_X+L1_norm*time_constant/nb_neuron
        lambda_majorant=intensity_function(Inf,X)
      }
      else
      {
        X=new_X
        lambda_majorant=intensity_function(Inf,X)
      } 
    }
  }
  good_index=which(points!=0)
  return(list(spike_train=points[good_index],type=type[good_index],intensity=intensity[good_index]))
}

steady_activity<-function(nb_neuron,vec_input_intensity,vec_refractory_period,vec_L1_norm,time_constant,nb_points=100000)
{
  # Fonction qui calcule l'intensite de l'etat stable pour une collection de parametres a chercher dans :
  # vec_inpu_intensity qui contient plusieurs valeurs d'intensite spontanee
  # vec_refractory_period qui contient plusieurs valeurs de periodes refractaires
  # vec_L1_norm qui contient plusieurs valeurs de force de connectivite
  steady_activities=array(0,c(length(vec_input_intensity),length(vec_refractory_period),length(vec_L1_norm)))
  spike_length=array(0,c(length(vec_input_intensity),length(vec_refractory_period),length(vec_L1_norm)))
  for (index_input in 1:length(vec_input_intensity))
  {
    input=vec_input_intensity[index_input]
    for (index_refrac in 1:length(vec_refractory_period))
    {
      refractory_period=vec_refractory_period[index_refrac]
      intensity_function=function(s,x){(input+x)*(s>refractory_period)}
      for (index_L1 in 1:length(vec_L1_norm))
      {
        L1_norm=vec_L1_norm[index_L1]
        a_infty=theoretical_steady_activity(input,L1_norm,refractory_period)
        X_init=L1_norm*a_infty # steady value of X
        is_refrac=rbinom(nb_neuron,1,a_infty*refractory_period) # tells if the neuron is initially in refractory period or not
        age_init=runif(nb_neuron,min = 0,max = refractory_period)*is_refrac + (refractory_period+1)*(1-is_refrac) # steady value of the age
        result=simul_Hawkes_expo_with_age(nb_neuron,refractory_period,intensity_function,L1_norm,time_constant,age_init,X_init,nb_points)
        points=result$spike_train
        steady_activities[index_input,index_refrac,index_L1]=length(points)/(max(points)*nb_neuron)
        spike_length[index_input,index_refrac,index_L1]=length(points)
      }
    }
  }
  return(list(steady_activity=steady_activities,number_of_spikes=spike_length))
}


#####################################################################################
#######  Calls and plots
#####################################################################################

# parameters
vec_input=exp(log(10)*seq(-6,2,length=30))
vec_refractory_period=5
vec_L1_norm=seq(0,2,length=7)
Nneur=1000
time_constant=0.005
nb_points=10000

# simulation
test=steady_activity(Nneur,vec_input,vec_refractory_period,vec_L1_norm,time_constant,nb_points)

# plot of activity in log/log scale
#pdf("activity_log_log.pdf",width=7,height=6)
plot(vec_input,rep(0,length(vec_input)),log="xy",ylim = c(10^(-6),1),axes=FALSE,ann=FALSE)
axis(cex.axis=1.2,1, at=10^(-6:2), lab=c(expression(10^-3),expression(10^-2), expression(10^-1), expression(1), expression(10),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))
axis(cex.axis=1.2,2, at=10^(-6:0), lab=c(expression(10^-2), expression(10^-1), expression(1), expression(10),expression(10^2),expression(10^3)))
box()
title(cex.lab=1.2,xlab = expression(paste(mu, ' (Hz)',sep="")),ylab = expression(paste(a[infinity], ' (Hz)')))
for (index_L1 in 1:length(vec_L1_norm))
{
  points(vec_input,test$steady_activity[,1,index_L1],cex=0.7)
  lines(vec_input,theoretical_steady_activity(vec_input,vec_L1_norm[index_L1],vec_refractory_period))
}
#dev.off()

# plot of activity in log/linear scale
#pdf("activity_log_linear.pdf",width=7,height=6)
plot(vec_input,rep(-1,length(vec_input)),log="x",ylim = c(10^(-6),0.2),axes=FALSE,ann=FALSE)
axis(cex.axis=1.2,1, at=10^(-6:2), lab=c(expression(10^-3),expression(10^-2), expression(10^-1), expression(1), expression(10),expression(10^2),expression(10^3),expression(10^4),expression(10^5)))
axis(cex.axis=1.2,2, at=seq(0,0.2,length=5), lab=seq(0,200,length=5))
box()
title(cex.lab=1.2,xlab = expression(paste(mu, ' (Hz)',sep="")),ylab = expression(paste(a[infinity], ' (Hz)')))
for (index_L1 in 1:length(vec_L1_norm))
{
  points(vec_input,test$steady_activity[,1,index_L1],cex=0.7)
  lines(vec_input,theoretical_steady_activity(vec_input,vec_L1_norm[index_L1],vec_refractory_period))
}
#dev.off()

# plot of sensitivity
beta_vec=c(0,0.01,0.2,0.5,0.6)
#pdf("sensitivity_test.pdf",width=7,height=6)
for (beta_index in 1:length(beta_vec))
{
  x=seq(0.01,2,length=1000)
  y=theoretical_sensitivity(beta_vec[beta_index],x,1)
  plot(x,y,ylim=c(0,min(max(y),500)+0.1*min(beta_index,4)),type='l',axes=FALSE,ann=FALSE,col=beta_index)
  axis(2,at=min(max(y),510),labels=floor(10*min(max(y),500))/10,col=beta_index,las=1)
  if (beta_index==2)
  {lines(c(0.9731099,0.9731099),c(0,max(y)),col='red')}
  par(new=T)
}
box()
title(cex.lab=1.2,xlab = expression(alpha),ylab = expression(sigma))
axis(cex.axis=1.2,1,seq(0,2,length=5))
axis(2,0,las=1)
legend('topright',c(expression(beta== 0),expression(beta== 0.01),expression(beta== 0.2),expression(beta== 0.5),expression(beta== 0.6)),lty=1,col=1:length(beta_vec))
#dev.off()

x_test=seq(0.01,2,length=100000)
y_test=y=theoretical_sensitivity(0.01,x_test,1)
approx_alpha_m=x_test[which(y_test==max(y_test))]
# 0.9731099
