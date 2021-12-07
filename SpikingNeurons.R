library(ggplot2)
library(ggridges)
library(IDPmisc)
library(GA)
library(GenSA)
library(viridis)
library(tuneR)
#install.packages("Rcpp")
library(Rcpp)
#install.packages("parallel")
#install.packages("doParallel")
library(parallel)
library(doParallel)
library(methods)
library(cowplot)

setwd("C:/Users/rapey/Desktop/Biomatematyka i Teoria Gier/R_IF/data")

jakDlugoMaDzialac = 60*60*10 # W sekundach, 60*60*2 oznacza 2 godziny
models = c("IF", "LIF", "QIF","EIF")[c(1,2,3,4)]

stepsize = 0.01

# _____________________ CONSTANTS AND CONSTANT FUNCTIONS ----

HHGain = readRDS( "HHgain.RData")

# HODGKIN HUXLEY CONSTANTS
C = 1
gna = 120
gk = 36
gl = 0.3 # m.mho/cm^2
Vna = -115 # 50 # -115 # mV 
Vk = 12 # -77 # 12 # mV
Vl = -10.613 # -54.387 # -10.613 # mV 

# HODGKIN HUXLEY INITIAL VALUES
V = 0.01
n = 0.3177323
m = 0.05295508
h = 0.5959924


alphanFunc = function( V ){
  0.01*(V + 10)/( exp((V + 10)/(10)) -1 ) 
}
betanFunc = function( V ){
  0.125*exp( V / 80 )
}
alphamFunc = function( V ){
  0.1*( V + 25 )/( exp( (V + 25)/10 ) - 1 )
}
betamFunc = function(V){
  4*exp(V/18)
}
alphahFunc = function(V){
  0.07*exp(V/20)
}
betahFunc = function(V){
  1/( exp( (V + 30)/10 ) +1 ) 
}

fn = function(V,n,m,h){
  alphanFunc(V)*(1 - n) - betanFunc(V)*n
}
fm = function(V,n,m,h){
  alphamFunc(V)*(1 - m) - betamFunc(V)*m
}
fh = function(V,n,m,h){
  alphahFunc(V)*(1 - h) - betahFunc(V)*h
}
Itot = function(I,V,n,m,h){
  (-I - gk*(n**4)*( V - Vk ) - gna*(m**3)*h*(V - Vna) - gl*(V - Vl ))/C
}

f =function(I,V,n,m,h){
  c( Itot(I,V,n,m,h), fn(V,n,m,h), fm(V,n,m,h), fh(V,n,m,h) )
}
HHFunction = f

alpha = function(V){ # gauusian curve 
  exp( (-(V*0.2-1)**2)*10)
}

IFFunction = function(I, V, C ){
  # FUNCTION F FOR PERFECT INTEGRATE AND FIRE
  I/C
}

LIFFunction = function(I, V, Vr, C, R){
  # FUNCTION F FOR LEAKY INTEGRATE AND FIRE  
  I/C - (V - Vr)/(R*C)
}

EIFFunction = function(I, V, Vr, C, R, delta, theta){
  # FUNCTION F FOR Exponential INTEGRATE AND FIRE  
  result = I/C - (V - Vr)/(R*C) + delta*exp( (V - theta)/delta )/(R*C)
  return( min( result, 10000 ) ) # numbers would get to big and break algorithm
}

QIFFunction = function(I, V, Vr, C, R, theta){
  # FUNCTION F FOR Quadratic INTEGRATE AND FIRE  
  I/C + (V - Vr)*(V - theta)/(R*C)
}

aIFFunction = function(I, V, w, a, Vr, C){
  # FUNCTION F FOR ADAPTIVE PERFECT INTEGRATE AND FIRE  
  c( I/C - w/C, a*(V - Vr) - w )
}

aLIFFunction = function(I, V, w, a, b, Vr, C){
  # FUNCTION F FOR LEAKY INTEGRATE AND FIRE  
  c( I/C - w/(C), a*b*(V - Vr) - a*w )
}

aEIFFunction = function(I, V, w, a, b, Vr, C, R, delta, theta){
  # FUNCTION F FOR Exponential INTEGRATE AND FIRE  
  dV = I/C - (V - Vr)/(C*R) - w/(C) + delta*exp( (V - theta)/delta )/(C)
  dV = min( dV, 10000 ) # numbers would get to big and break algorithm
  dw =  a*b*(V - Vr) - a*w 
  dw = min( dw, 10000 )
  c( dV , dw )
}

aQIFFunction = function(I, V, w, a, b,Vr, C, R, theta){
  # FUNCTION F FOR Quadratic INTEGRATE AND FIRE  
  dV =  I/C + (V - Vr)*(V - theta)/(R*C) - w
  dV = min( dV, 10000 ) # numbers would get to big and break algorithm
  dw = a*b*(V - Vr) - a*w
  dw = min( dw, 10000 )
  c( dV , dw  )
}


GainLIF = function(I, V0 = -12, Vr = 0, threshold = 10.14, R = 1, C = 1, refractoryPeriod = 10){
  freq = 1/( refractoryPeriod -R*C* log( 1 - (V0 - threshold)/( V0 - Vr - R*I ) ) )
  if( !is.nan(freq) ){
    return(freq)
  }else{
    return(0)
  }
}


GainIF = function(I, V0 = -12, threshold = 10.14, C = 1, refractoryPeriod = 10){
  return(1/( refractoryPeriod + C*(threshold - V0)/(I) ))
}


modelTemplate = list( name = NA_character_,
                      parameters = c(   "C"  ,    "V0"  ,   "R"  ,    "Vr"  ,  "delta"  ,   "theta"  ,   "spike"   , "refractoryPeriod", "threshold" , "a"   ,  "b" , "d" ),
                      lowerBound = c(  0.617 ,   -20    ,  0.01  ,    -10   ,     0.01  ,     6      ,     50      ,        0          ,       2     ,   0   ,   0  , -3  ),
                      upperBound = c(   5    ,     0    ,    10  ,     10   ,     10    ,     12     ,     100     ,        15         ,      20     ,   1   ,   2  , 10  ),
                      default    = c(   1    ,    -11   ,    5   ,      0   ,     2     ,     10     ,     100     ,        0          ,     10.14   , 0.02  ,   1  ,  8  ), 
                      subPars    = c(   "C"  ,              "R"  ,    "Vr"                                                                           , "a"   ,  "b"       ), 
                      overPars   = c(             "V0"  ,                      "delta"  ,   "theta"  ,   "spike"   , "refractoryPeriod", "threshold"                , "d" ),
                      dimensions = 1,
                      dV = NULL )

names(modelTemplate$lowerBound) = modelTemplate$parameters
names(modelTemplate$upperBound) = modelTemplate$parameters
names(modelTemplate$default) = modelTemplate$parameters

# ___________________________ CLASSES  ----

# ___________________________ CLASSES - NEURON MODEL ----

setClass( "NeuronModel" ,
  slots = c(
    name       = "character" ,
    parameters = "character" ,
    subPars    = "character" ,
    overPars   = "character" ,
    upperBound = "numeric"   ,
    lowerBound = "numeric"   ,
    default    = "numeric"   ,
    dV         = "function"  
  ),
  prototype = list(
    name       = NA_character_ ,
    parameters = NA_character_ ,
    subPars    = NA_character_ ,
    overPars   = NA_character_ ,
    upperBound = NA_real_      ,
    lowerBound = NA_real_      ,
    default    = NA_real_      ,
    dV         = NULL            
  ),
  validity = function( object ){
    numOfParameters = length( object@parameters )
    if( length( object@upperBound      ) != numOfParameters  ){ return("length of @upperBound is not equal to length of @parameters") }
    if( length( object@lowerBound      ) != numOfParameters  ){ return("length of @lowerBound is not equal to length of @parameters") }
    if( length( object@default         ) != numOfParameters & length( object@default )  != 0 ){ return("length of @default is not equal to length of @parameters") }
    return(TRUE)
  }
)



setMethod( f = "initialize", signature = "NeuronModel",
           definition = function( .Object,  name , parameters , subPars = NA , overPars = NA , upperBound = NA   , lowerBound = NA , default = NA , dV ){
             .Object@name = name
             .Object@parameters = parameters
             .Object@dV = dV
             
             if( is.na( subPars  ) ){ .Object@subPars = intersect( modelTemplate$subPars, parameters ) }else{ .Object@subPars = subPars }
             if( is.na( overPars ) ){ .Object@overPars = intersect( modelTemplate$overPars, parameters ) }else{ .Object@overPars = overPars }
             if( is.na(upperBound) ){ .Object@upperBound = modelTemplate$upperBound[ intersect( modelTemplate$parameters, parameters ) ] }else{ .Object@upperBound = upperBound }
             if( is.na(lowerBound) ){ .Object@lowerBound = modelTemplate$lowerBound[ parameters ] }else{ .Object@lowerBound = lowerBound }
             
             .Object@default = modelTemplate$default[ intersect( modelTemplate$parameters, parameters ) ] 
             if( !is.na(default) ){ 
               for( par in names(default) ){ .Object@default[ par ] = default[par]  }
             }
             
             
             validObject(.Object)
             return(.Object)
           }
)


HH   = new( "NeuronModel", name = "HH"   , dV = HHFunction   , parameters = c( "C"                                                                                                        ) )

IF   = new( "NeuronModel", name = "IF"   , dV = IFFunction   , parameters = c( "C" , "V0" , "R" ,                                      "refractoryPeriod" , "threshold"                   ) )
LIF  = new( "NeuronModel", name = "LIF"  , dV = LIFFunction  , parameters = c( "C" , "V0" , "R" , "Vr" ,                               "refractoryPeriod" , "threshold"                   ) )
QIF  = new( "NeuronModel", name = "QIF"  , dV = QIFFunction  , parameters = c( "C" , "V0" , "R" , "Vr" ,           "theta" , "spike" , "refractoryPeriod"                                 ) )
EIF  = new( "NeuronModel", name = "EIF"  , dV = EIFFunction  , parameters = c( "C" , "V0" , "R" , "Vr" , "delta" , "theta" , "spike" , "refractoryPeriod"                                 ) )
aIF  = new( "NeuronModel", name = "aIF"  , dV = aIFFunction  , parameters = c( "C" , "V0" , "R" ,                                      "refractoryPeriod" , "threshold" , "a" , "b" , "d" ) )
aLIF = new( "NeuronModel", name = "aLIF" , dV = aLIFFunction , parameters = c( "C" , "V0" , "R" , "Vr" ,                               "refractoryPeriod" , "threshold" , "a" , "b" , "d" ) )
aQIF = new( "NeuronModel", name = "aQIF" , dV = aQIFFunction , parameters = c( "C" , "V0" , "R" , "Vr" ,           "theta" , "spike" , "refractoryPeriod" ,               "a" , "b" , "d" ) )
aEIF = new( "NeuronModel", name = "aEIF" , dV = aEIFFunction , parameters = c( "C" , "V0" , "R" , "Vr" , "delta" , "theta" , "spike" , "refractoryPeriod" ,               "a" , "b" , "d" ) )
allModels = c(HH , IF , LIF , QIF , EIF , aIF , aLIF , aQIF , aEIF )
SubModels = list( EIF = LIF, aEIF = aLIF , QIF = LIF, aQIF = aLIF)

# ___________________________ CLASSES CURRENT AND CURRENT TYPE  ----

setClass( "CurrentType",
          slots = c(
            name       = "character" ,
            parameters = "character" ,
            default    = "numeric"
          )
          )


setMethod(f = "initialize", signature = "CurrentType",
          definition = function( .Object, name, parameters, default ){
            .Object@name = name
            .Object@parameters = parameters
            .Object@default = default
            names(.Object@default) = parameters
            validObject(.Object)
            return(.Object)
          }
)


setClass( "Current",
  slots = c(
    I          = "numeric"     ,
    time       = "numeric"     ,
    stepsize   = "numeric"     ,
    type       = "CurrentType" ,
    parameters = "numeric" 
  )
)

setMethod(f = "initialize", signature = "Current",
          definition = function( .Object, I = NA , time = 200, stepsize = 0.01, type, parameters = NA){
            .Object@time = time
            .Object@stepsize = stepsize
            .Object@type = type
            
            if(!is.na(parameters)){
              .Object@parameters = parameters
            }else{
              .Object@parameters = type@default
            }
            
            names(.Object@parameters) = type@parameters
            
            if( type@name == "constant" ){
              .Object@I = rep( parameters["I"] , time/stepsize )
            }
            if( type@name == "preset"   ){
              .Object@I = I
            }
            
            if( type@name == "weiner" ){
              J = 0
              I = numeric()
              # Weiner Process
              for(t in 2:(time/stepsize) ){
                J = J + rnorm(1)*.Object@parameters[["sigma"]]
                # LIMITING THE CURRENT
                if( !is.nan(.Object@parameters[["lowerLimit"]]) & is.numeric(.Object@parameters[["lowerLimit"]]) ){
                  if(J < .Object@parameters[["lowerLimit"]]){J = .Object@parameters[["lowerLimit"]]}
                }
                
                if( !is.nan(.Object@parameters[["upperLimit"]]) & is.numeric(.Object@parameters[["upperLimit"]]) ){
                  if(J > .Object@parameters[["upperLimit"]]){J = .Object@parameters[["upperLimit"]]}
                }
                I = c(I,J)
              }
              .Object@I = I
            }
            if( type@name == "whiteNoise" ){
              .Object@I = rnorm( time/stepsize, mean = .Object@parameters["mu"]  , sd = .Object@parameters["sigma"])
            }
            if( type@name == "pinkNoise" ){
              .Object@I = noise(kind = "pink", duration = time/stepsize , samp.rate = 1/stepsize)@left*.Object@parameters["sigma"] + .Object@parameters["mu"]
            }
            if( type@name == "Orstein-Uhlenbeck" ){
              dW = rnorm( time/stepsize ,0, .Object@parameters[["sigma"]]*sqrt(stepsize) )
              I = .Object@parameters[["nu"]]
              for(i in 2:(time/stepsize) ){
                I[i] = I[i-1] - .Object@parameters[["lambda"]]*( I[i-1] - .Object@parameters[["nu"]] ) + dW[i]
              }
              .Object@I = I
            }
          validObject(.Object)
          return(.Object) 
            
          }
          )

# ___________________________ CLASSES  PARAMETER POINTS  ----

setClass("ParameterPoint",
         slots = c(
           model  = "NeuronModel" ,
           values = "numeric"     ,
           score  = "numeric"
         ),
         validity = function( object ){
           if( !all(names(object@values) %in%  object@model@parameters) ){ return( "parameters don't fit the model" ) }
           return(TRUE)
         }
)

setMethod(f = "initialize", signature = "ParameterPoint",
          definition = function( .Object, model, values , score = NA_real_){
            .Object@values  = values
            .Object@model   = model
            .Object@score   = score
            validObject(.Object)
            return(.Object)
          }
)

# ___________________________ CLASSES  SOLUTION   ----

setClass("Solution", 
         slots = c(
           V          =     "numeric" ,
           spikes     =     "numeric" ,
           current    =     "Current" ,
           model      = "NeuronModel" ,
           parameters = "ParameterPoint"
         )
)

setMethod(f = "initialize", signature = "Solution",
          definition = function( .Object, V = NA, spikes = NA, current , model , parameters ){
            .Object@current    = current
            .Object@model      = model
            .Object@parameters = parameters
            .Object@V          = V
            .Object@spikes     = spikes
            validObject(.Object)
            return(.Object)
          }
          )

# ___________________________ CLASSES  SCORE AND SCORE MEASURE ----

setClass( "ScoreMeasure",
          slots = c(
            name        = "character" ,
            maxScore    = "numeric"   ,
            parameters  = "character" ,
            default     = "numeric"   
          )
)

setClass("Score",
         slots = c(
           value   = "numeric"      ,
           measure = "ScoreMeasure" ,
           parameters = "numeric"
         )
)

# ___________________________ SOLVE FUNCTION  ----


solve = function( current ,  parameterPoint = NULL , model = NULL, refractoryPeriod = 0,
                     blockSpikes = FALSE){
  
  if(is.null(model) & !is.null(parameterPoint) ){ model = parameterPoint@model }
  
  for(par in model@parameters)
    assign(par, model@default[par])
  
  if(is.null(parameterPoint)){ parameterPoint = new("ParameterPoint", model = model, values = model@default) }
  
  parameters = parameterPoint@values
  # assigning all the parameters from the argument list(or vector) "parameters"
  if(length(parameters) != 0){
    for(par in names(parameters))
      assign(par, parameters[[par]] )
  }
  
  if( class(current) == "numeric" ){ current = new("Current", parameters = c("I" = current), type = constantCurrent, stepsize = 0.01) }
  if( class(current) == "Current" ){ I = current@I }
  P = 0
  n = 0.3177323
  m = 0.05295508
  h = 0.5959924
  
  V = 0.01
  stepsize = current@stepsize
  spikeTimes = numeric()
  steps = length(I)
  refractoryPeriodLength = refractoryPeriod/stepsize
  
  # INTEGRATE AND FIRE 
  if(model@name == "IF"){
    
    for(t in 2:(steps) ){
      
      if( t - max(spikeTimes[length(spikeTimes)], -Inf) > refractoryPeriodLength ){
        
        k1 = IFFunction(I[t], V[t-1] , C)
        k2 = IFFunction(I[t], V[t-1] + k1[1]*(stepsize/2), C)
        k3 = IFFunction(I[t], V[t-1] + k2[1]*(stepsize/2), C)
        k4 = IFFunction(I[t], V[t-1] + k3[1]*(stepsize  ), C)
        V[t] = V[t-1] + stepsize*(k1 + 2*k2 + 2*k3 + k4)/6
        
        if(V[t] > threshold){
          V[t] = V0 
          spikeTimes[length(spikeTimes) + 1] = t
        }
      }else{
        V[t] = V[t-1]
      }
    }
  }
  
  
  # ADAPTIVE INTEGRATE AND FIRE 
  if(model@name == "aIF"){
    w = 0
    for(t in 2:( length(I) ) ){
      
      if( t - max(spikeTimes[length(spikeTimes)], -Inf) > refractoryPeriodLength ){
        
        k1 = aIFFunction(I[t], V[t-1] , w[t-1], a, Vr, C)
        k2 = aIFFunction(I[t], V[t-1] + k1[1]*(stepsize/2) , w[t-1] + k1[2]*(stepsize/2) , a , Vr , C)
        k3 = aIFFunction(I[t], V[t-1] + k2[1]*(stepsize/2) , w[t-1] + k2[2]*(stepsize/2) , a , Vr , C)
        k4 = aIFFunction(I[t], V[t-1] + k3[1]*(stepsize  ) , w[t-1] + k3[2]*(stepsize  ) , a , Vr , C)
        V[t] = V[t-1] + stepsize*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6
        w[t] = w[t-1] + stepsize*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6
        
        if(V[t] > threshold){
          V[t] = V0 
          w[t] = d
          spikeTimes[length(spikeTimes) + 1] = t
        }
      }else{
        w[t] = w[t-1]
        V[t] = V[t-1]
      }
    }
  }
  
  # LEAKY INTEGRATE AND FIRE 
  if(model@name == "LIF" ){
    for(t in 2:(steps) ){
      
      if( t - max(spikeTimes[length(spikeTimes)], -Inf) > refractoryPeriodLength ){
        
        k1 = LIFFunction(I[t], V[t-1] , Vr, C, R)
        k2 = LIFFunction(I[t], V[t-1] + k1[1]*(stepsize/2), Vr, C, R)
        k3 = LIFFunction(I[t], V[t-1] + k2[1]*(stepsize/2), Vr, C, R)
        k4 = LIFFunction(I[t], V[t-1] + k3[1]*(stepsize  ), Vr, C, R)
        V[t] = V[t-1] + stepsize*(k1 + 2*k2 + 2*k3 + k4)/6
        if(V[t] > threshold & !blockSpikes){
          V[t] = V0 
          spikeTimes[length(spikeTimes) + 1] = t
        }
      }else{
        V[t] = V[t-1]
      }
      
    }
  }
  
  # LEAKY INTEGRATE AND FIRE 
  if(model@name == "aLIF" ){
    w = 0
    for(t in 2:(steps) ){
      
      if( t - max(spikeTimes[length(spikeTimes)], -Inf) > refractoryPeriodLength ){
        
        k1 = aLIFFunction(I[t], V[t-1] , w[t-1] , a , b , Vr , C )
        k2 = aLIFFunction(I[t], V[t-1] + k1[1]*(stepsize/2), w[t-1] + k1[2]*(stepsize/2) , a , b , Vr , C )
        k3 = aLIFFunction(I[t], V[t-1] + k2[1]*(stepsize/2), w[t-1] + k2[2]*(stepsize/2) , a , b , Vr , C )
        k4 = aLIFFunction(I[t], V[t-1] + k3[1]*(stepsize  ), w[t-1] + k3[2]*(stepsize  ) , a , b , Vr , C )
        V[t] = V[t-1] + stepsize*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6
        w[t] = w[t-1] + stepsize*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6
        if( is.na(V[t]) ){ V[t] = threshold + 1 }
        if(V[t] > threshold & !blockSpikes){
          V[t] = V0 
          w[t] = d
          spikeTimes[length(spikeTimes) + 1] = t
        }
      }else{
        w[t] = w[t-1]
        V[t] = V[t-1]
      }
      
    }
  }
  
  # EXPONENTIAL INTEGRATE AND FIRE
  if(model@name == "EIF" ){
    for(t in 2:(steps) ){
      
      if( t - max(spikeTimes[length(spikeTimes)], -Inf) > refractoryPeriodLength ){
        k1 = EIFFunction(I[t], V[t-1] , Vr, C, R, delta, theta)
        k2 = EIFFunction(I[t], V[t-1] + k1[1]*(stepsize/2),  Vr, C, R, delta, theta)
        k3 = EIFFunction(I[t], V[t-1] + k2[1]*(stepsize/2),  Vr, C, R, delta, theta)
        k4 = EIFFunction(I[t], V[t-1] + k3[1]*(stepsize  ),  Vr, C, R, delta, theta)
        V[t] = V[t-1] + stepsize*(k1 + 2*k2 + 2*k3 + k4)/6
        #if( is.na(V[t]) ){ print(V); print( list(k1=k1,k2=k2,k3=k3,k4=k4,v=V[t-1],I= I[t], V0 = V0, R= R, C=C,delta=delta, theta=theta,V0 = V0) ) }
        if( is.na(V[t]) ){ V[t] = spike + 1 }
        if(V[t] > spike){
          V[t] = V0 
          spikeTimes[length(spikeTimes) + 1] = t
        }
      }else{
        V[t] = V[t-1]
      }
      
    }
  }
  
  # EXPONENTIAL INTEGRATE AND FIRE
  if(model@name == "aEIF" ){
    w = 0
    for(t in 2:(steps) ){
      
      if( t - max(spikeTimes[length(spikeTimes)], -Inf) > refractoryPeriodLength ){
        k1 = aEIFFunction(I[t], V[t-1] , w[t-1], a, b , Vr, C, R, delta, theta)
        k2 = aEIFFunction(I[t], V[t-1] + k1[1]*(stepsize/2), w[t-1] + k1[2]*(stepsize/2), a, b , Vr, C, R, delta, theta)
        k3 = aEIFFunction(I[t], V[t-1] + k2[1]*(stepsize/2), w[t-1] + k2[2]*(stepsize/2), a, b , Vr, C, R, delta, theta)
        k4 = aEIFFunction(I[t], V[t-1] + k3[1]*(stepsize  ), w[t-1] + k3[2]*(stepsize  ), a, b , Vr, C, R, delta, theta)
        V[t] = V[t-1] + stepsize*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6
        w[t] = w[t-1] + stepsize*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6 
        #if( is.na(V[t]) ){ print(V); print( list(k1=k1,k2=k2,k3=k3,k4=k4,v=V[t-1],I= I[t], EIFV0 = V0, EIFR= R, EIFC=C,EIFdelta=delta, EIFtheta=theta,EIFV0 = V0) ) }
        if( is.na(V[t]) ){ V[t] = spike + 1 }
        if(V[t] > spike){
          V[t] = V0 
          w[t] = w[t] + d
          spikeTimes[length(spikeTimes) + 1] = t
        }
      }else{
        w[t] = w[t-1]
        V[t] = V[t-1]
      }
      
    }
  }
  
  
  # QUADRATIC INTEGRATE AND FIRE 
  if(model@name == "QIF" ){
    for(t in 2:(steps) ){
      
      if( t - max(spikeTimes[length(spikeTimes)], -Inf) > refractoryPeriodLength ){
        k1 = QIFFunction(I[t], V[t-1] , Vr, C, R , theta )
        k2 = QIFFunction(I[t], V[t-1] + k1[1]*(stepsize/2), Vr, C, R , theta)
        k3 = QIFFunction(I[t], V[t-1] + k2[1]*(stepsize/2), Vr, C, R , theta)
        k4 = QIFFunction(I[t], V[t-1] + k3[1]*(stepsize  ), Vr, C, R , theta)
        V[t] = V[t-1] + stepsize*(k1 + 2*k2 + 2*k3 + k4)/6
        if(V[t] > spike){
          V[t] = V0 
          spikeTimes[length(spikeTimes) + 1] = t
        }
      }else{
        V[t] = V[t-1]
      }
      
    }
  }
  
  
  # QUADRATIC INTEGRATE AND FIRE 
  if(model@name == "aQIF" ){
    w = 0
    for(t in 2:(steps) ){
      if( t - max(spikeTimes[length(spikeTimes)], -Inf) > refractoryPeriodLength ){
        k1 = aQIFFunction(I[t], V[t-1] , w[t-1] , a , b , Vr, C, R , theta )
        k2 = aQIFFunction(I[t], V[t-1] + k1[1]*(stepsize/2), w[t-1] + k1[2]*(stepsize/2) , a , b , Vr, C, R , theta)
        k3 = aQIFFunction(I[t], V[t-1] + k2[1]*(stepsize/2), w[t-1] + k2[2]*(stepsize/2) , a , b , Vr, C, R , theta)
        k4 = aQIFFunction(I[t], V[t-1] + k3[1]*(stepsize  ), w[t-1] + k3[2]*(stepsize  ) , a , b , Vr, C, R , theta)
        V[t] = V[t-1] + stepsize*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6
        w[t] = w[t-1] + stepsize*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6
        if( is.na(V[t]) ){ V[t] = spike + 1 }
        if(V[t] > spike){
          w[t] = d
          V[t] = V0 
          spikeTimes[length(spikeTimes) + 1] = t
        }
      }else{
        V[t] = V[t-1]
        w[t] = w[t-1]
      }
      
    }
  }
  
  
  # HODGKIN-HUXLEY - 4-TH ORDER RUNGE-KUTTA METHOD
  if(model@name == "HH" ){
    for(t in 2:(steps) ){
      if(blockSpikes){ {m[t-1] = 0} }
      
      k1 = f(I[t], V[t-1] , n[t-1] , m[t-1] , h[t-1] )
      k2 = f(I[t], V[t-1] + k1[1]*(stepsize/2), n[t-1] + k1[2]*(stepsize/2), m[t-1] + k1[3]*(stepsize/2), h[t-1] + k1[4]*(stepsize/2))
      k3 = f(I[t], V[t-1] + k2[1]*(stepsize/2), n[t-1] + k2[2]*(stepsize/2), m[t-1] + k2[3]*(stepsize/2), h[t-1] + k2[4]*(stepsize/2))
      k4 = f(I[t], V[t-1] + k3[1]*(stepsize  ), n[t-1] + k3[2]*(stepsize  ), m[t-1] + k3[3]*(stepsize  ), h[t-1] + k3[4]*(stepsize  ))
      
      V[t] = V[t-1] + stepsize*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6
      n[t] = n[t-1] + stepsize*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6
      m[t] = m[t-1] + stepsize*(k1[3] + 2*k2[3] + 2*k3[3] + k4[3])/6
      h[t] = h[t-1] + stepsize*(k1[4] + 2*k2[4] + 2*k3[4] + k4[4])/6
      
    }
    
    if( any(-V > 40) ){ 
      spikeTimes = peaks(-V, thr = 50)["x"][[1]]
    }else{
      spikeTimes = numeric()
    } 
    V = -V
  }
  
  return( new("Solution", model = model, current = current, parameters = parameterPoint , V = V, spikes = spikeTimes) )
}



# ___________________________ CLASSES ALGORITHM AND ALGORITHM TYPE ----

setClass("AlgorithmType",
         slots = c(
          name       = "character" ,
          parameters = "character" ,
          default    = "numeric"   ,
          run        = "function"  
         )
)

setClass("Algorithm",
         slots = c(
          type            = "AlgorithmType" ,
          parameterVals   = "numeric"       ,
          estimatedTime   = "numeric"       ,
          fixedParameters = "numeric"       
         )
)

# ___________________________ CLASSES ESTIMATE ----

setClass ("Estimate",
          slots = c(
            parameterPoint = "ParameterPoint" , 
            algorithm      = "Algorithm"      ,
            scoreTr        = "numeric"        ,
            scoreTe        = "numeric"        ,
            GA             = "ga"
          )
)

setMethod(f = "print", signature = "Estimate",
          definition = function(x){
            scoreTr   = round( x@scoreTr  , 2 )
            scoreTe   = round( x@scoreTe  , 2 )
            vals      = x@parameterPoint@values
            model     = x@parameterPoint@model@name
            algorithm = paste( "[",  x@algorithm@type@name , ": p=" , x@algorithm@parameterVals[["popSize"]], " m=" , x@algorithm@parameterVals[["maxIter"]], " r=" , x@algorithm@parameterVals[["run"]] , "]", sep = "" )
            print(paste(  "[" , model , " |Tr:" ,  scoreTr , " |Te:" , scoreTe , " |vals=", vals , "]" , sep = ""  ))
          }
          )

setMethod(f = "show", signature = "Estimate",
          definition = function(x){
            scoreTr   = round( x@scoreTr  , 2 )
            scoreTe   = round( x@scoreTe  , 2 )
            vals      = round(x@parameterPoint@values , 2)
            model     = x@parameterPoint@model@name
            algorithm = paste( "[",  x@algorithm@type@name , ": p=" , x@algorithm@parameterVals[["popSize"]], " m=" , x@algorithm@parameterVals[["maxIter"]], " r=" , x@algorithm@parameterVals[["run"]] , "]", sep = "" )
            print(paste(  "[" , model , " |Tr:" ,  scoreTr , " |Te:" , scoreTe , " |algorithm:" , algorithm ,"]" , sep = ""  ))
            print(vals)
          }
)

# ___________________________ CLASSES PROBLEM ----

setClass("Problem",
         slots = c(
           dataTrOver = "Solution"  ,
           dataTrSub  = "Solution"  ,
           dataTeOver = "Solution"  ,
           dataTeSub  = "Solution"  ,
           mode       = "character" ,
           results    = "list"      
         )
)

setMethod(f = "initialize", signature = "Problem",
          definition = function( .Object, dataTrOver, dataTrSub , dataTeOver , dataTeSub, mode = NA_character_){
            .Object@dataTrOver  =  dataTrOver
            .Object@dataTrSub  = dataTrSub 
            .Object@dataTeOver  =  dataTeOver
            .Object@dataTeSub  = dataTeSub 
            .Object@mode  = mode
            
            .Object@results = list()
            for( model in allModels ){
              .Object@results[[ model@name ]] = list( "sub" = list() , "over" = list() )
            }
            
            validObject(.Object)
            return(.Object)
          }
)


# ___________________________ FITNESS FUNCTIONS AND CURVE CALCULATION ----

calculateGain = function( parameterPoint , algorithm ,  model , fast = FALSE ){
  
  if(class(parameterPoint) == "numeric"){ parameterPoint = new("ParameterPoint", model = model, values = parameterPoint) }
  
  #default values:
  resolution = 0.1  
  minCurrent = 1
  maxCurrent = 15
  
  for( par in names(algorithm@parameterVals)){
    assign( par , algorithm@parameterVals[[par]] )
  }
  currentVals = (minCurrent/resolution):(maxCurrent/resolution)*resolution
  gain = numeric()
  
  if( algorithm@parameterVals["noise"] == T ){
    
    for( i in 1:length(currentVals) ){
      freq = numeric() 
      
      for( loop in 1:5 ){
        
        current = new("Current", type = whiteNoise , time = 200, stepsize = 0.01 , 
                      parameters = c( "mu" = currentVals[i], "sigma" = algorithm@parameterVals[["sigma"]] )  )
        sol = solve( current = current, model = model , parameterPoint = parameterPoint  )
        howManySpikes = length(sol@spikes)
        if( howManySpikes > 2 ){ 
          freq = c( freq , (howManySpikes - 1)/(  sol@spikes[howManySpikes] - sol@spikes[2] )/current@stepsize )
        }else{ freq = c( freq, 0 ) }
      }
      freq = mean( freq )
      curr = as.character(i)
      gain = c( gain , c( curr = freq)  )
    }
    
  }
  if( ( model@name != "LIF" & model@name != "IF" ) & algorithm@parameterVals["noise"] == F ){
    
    if( fast ){
      for( i in 1:length(currentVals) ){
        current = new("Current", type = constantCurrent , parameters = c("I" = currentVals[i]) , time = 50, stepsize = 0.01  )
        sol = solve( current = current, model = model , parameterPoint = parameterPoint )
        howManySpikes = length(sol@spikes)
        if( howManySpikes >= 2 ){ 
          freq = (howManySpikes - 1)/(  sol@spikes[howManySpikes] - sol@spikes[1] )/current@stepsize
        }else{ freq = 0 }
        curr = as.character(i)
        gain = c( gain , c( curr = freq)  )
      }
    }else{
      for( i in 1:length(currentVals) ){
        current = new("Current", type = constantCurrent , parameters = c("I" = currentVals[i]) , time = 200, stepsize = 0.01  )
        sol = solve( current = current, model = model , parameterPoint = parameterPoint  )
        howManySpikes = length(sol@spikes)
        if( howManySpikes > 2 ){ 
          freq = (howManySpikes - 2)/(  sol@spikes[howManySpikes] - sol@spikes[2] )/current@stepsize
        }else{ freq = 0 }
        
        curr = as.character(i)
        gain = c( gain , c( curr = freq)  )
      }
    }
    
  }
  if( model@name == "LIF" &  algorithm@parameterVals["noise"] == F ){
    for( i in 1:length(currentVals) ){
      for(par in names( model@default          ) ){ assign( par , model@default[par]         ) }
      for(par in names( parameterPoint@values  ) ){ assign( par , parameterPoint@values[par] ) }
      freq = GainLIF(I = currentVals[i] , V0 = V0, Vr = Vr , threshold = threshold , R = R, C = C, refractoryPeriod = refractoryPeriod )
      gain = c( gain ,  freq )
    }
  }
  if( model@name == "IF" &  algorithm@parameterVals["noise"] == F ){
    for( i in 1:length(currentVals) ){
      for(par in names( model@default          ) ){ assign( par , model@default[par]         ) }
      for(par in names( parameterPoint@values  ) ){ assign( par , parameterPoint@values[par] ) }
      freq = GainIF(I = currentVals[i] , V0 = V0, threshold = threshold , C = C, refractoryPeriod = refractoryPeriod )
      curr = as.character(i)
      gain = c( gain , c( curr = freq)  )
    }
  }
  return(gain)
  
}

fitnessAll  = function( parameterVector, problem , model , mode = "over", testData = FALSE , parameterNames = NA, fixedParameters = NA, subModel = F){
  if(mode == "sub"  ){ blockSpikes = T ; data = problem@dataTrSub  ; scoreMeasure = integral    ; if(testData){ data = problem@dataTeSub  } }
  if(mode == "over" ){ blockSpikes = F ; data = problem@dataTrOver ; scoreMeasure = coincidence ; if(testData){ data = problem@dataTeOver } }
  
  
  if( !any(is.na(parameterNames)) ){names(parameterVector) = parameterNames}
  if( !any(is.na(fixedParameters)) ){ 
    parameterVector = c( parameterVector , fixedParameters[ setdiff( names(fixedParameters) , parameterNames ) ] ) 
    parameterVector = parameterVector[ intersect( names(parameterVector) , model@parameters  ) ]
  }
  
  parameterVector = c( parameterVector , model@default[ setdiff( model@parameters , names(parameterVector) ) ]  )
  
  parameterPoint = new("ParameterPoint", values = parameterVector, model = model)
  solIF = solve( current = data@current , parameterPoint = parameterPoint , model = model , blockSpikes = blockSpikes )
  
  evaluate( scoreMeasure , solIF = solIF, solHH = data )@value
}
gainFitness = function( parameterVector, problem , model , algorithm , parameterNames = NA , fixedParameters = NA, fast = FALSE){
  
  if( length(parameterNames) == 0 ){ parameterNames = NA }
  if( !any(is.na(parameterNames)) ){ names(parameterVector) = parameterNames }
  if( !any(is.na( fixedParameters )) ){
    parameterVector = c( parameterVector , fixedParameters[ intersect( names(fixedParameters) , model@parameters ) ] )
  }
  
  #default values:
  resolution = 0.1  
  minCurrent = 1 
  maxCurrent = 15
  sigma = 11
  for( par in names(algorithm@parameterVals)){
    assign( par , algorithm@parameterVals[[par]] )
  }
  
  parameterPoint = new("ParameterPoint", values = parameterVector , model = model)
  gain = calculateGain( parameterPoint = parameterPoint , algorithm = algorithm , model = model , fast = fast )
  currentVals = (minCurrent/resolution):(maxCurrent/resolution)*resolution
  score = 0
  for( i in 1:length(currentVals) ){
    score = score - ( 100*abs( gain[ i ] - HHGain[ currentVals[i]*10 , sigma*2 ] )**2  )
  }
  score = score/length(currentVals)
  return(score)
}

# ________________________ OBJECT DEFINITONS
# ________________________ OBJECT DEFINITONS ALGORITHM TYPES ----

GAAlgorithm      = new("AlgorithmType", name = "GA"   , parameters = c("maxIter", "popSize", "run"), default = c( "maxIter" = 25, "popSize" = 50, run = 10))
GainFitAlgorithm = new("AlgorithmType", name = "gain" , 
                       parameters = c( "noise" , "resolution" , "minCurrent", "maxCurrent", "sigma" , "maxIter", "popSize", "run") ,
                       default    = c( noise = F, resolution = 0.1 , minCurrent = 4, maxCurrent = 15, sigma = 11 , "maxIter" = 25, "popSize" = 20, run = 10) )


# ________________________ OBJECT DEFINITONS MEASURES ----

coincMeasure = function( solIF, solHH, parameters = NA ){
  
  for( par in names(parameters) ){
    assign( par, parameters[par] )
  }
  
  spikesIF = solIF@spikes
  spikesHH = solHH@spikes
  time     = solIF@current@time
  stepsize = solIF@current@stepsize
  
  if( delta > 4 || delta < 1 ){ warning( paste( "Delta =", delta, "it is advised that delta be within [1,4]" ) ) }
  worstScore = -1
  
  if( length(spikesIF) == 0 & length(spikesHH) == 0 ){return(0)}
  
  vsim = length(spikesIF)/time
  
  Ncoinc = 0 # the number of coincidental spikes (spikes of the HH model that have an IF model close)
  for( i in 1:length(spikesHH) ){
    Ncoinc = Ncoinc + any( abs(spikesIF - spikesHH[i]) < delta/stepsize )
  }
  
  Ndata = length(spikesHH)
  NcoincPoiss = 2*vsim*delta*Ndata
  Nif = length(spikesIF)
  N = 1 - 2*(Ndata/time)*delta
  
  lambda = (1/N)*( Ncoinc - NcoincPoiss )/( Ndata/2 + Nif/2 )
  names(lambda) = NULL
  
  if( is.na(lambda) | is.nan(lambda) | is.null(lambda) | is.infinite(lambda) ){ lambda = worstScore }
  
  if( lambda > 1 | lambda < -1 ){
    warning(paste( "something's wrong, the score func is lambda =" , lambda ))
    setwd("C:/Users/rapey/Desktop/Biomatematyka i Teoria Gier/R_IF/data")
    report = list( Ncoinc = Ncoinc, NcoincPoiss = NcoincPoiss, Ndata = Ndata, Nif = Nif, N = N, spikesIF = spikesIF, spikesHH = spikesHH, lambda = lambda)
    saveRDS(report, file = paste("bug", runif(1), sep = "" ))
  }
  
  lambda
}
intMeasure = function( solIF, solHH , parameters = NA ){
  VIF = solIF@V
  V   = solHH@V
  if( length(solHH@spikes) != 0 ){
    V = ifelse( V>10.14 , 10.14 , V )
  }
  score = sum( abs(V - VIF) )/length(V)
  return(-score)
}

coincidence = new( "ScoreMeasure", name = "coincidence", maxScore = 1, parameters = c("delta"), default = c("delta" = 2) )
integral    = new( "ScoreMeasure", name = "integral"   , maxScore = 0, parameters = NA_character_, default = NA_real_    )
evaluate = function( .Object , solIF, solHH, parameters = NA_real_){
  
  for(par in names(.Object@default)){
    assign( par, .Object@default[par] )
  }
  
  for(par in names(parameters)){
    assign( par, parameters[par] )
  }
  
  if( is.na(parameters) ){ parameters = .Object@default }
  
  if( .Object@name == "coincidence" ){
    scoreValue = coincMeasure( solIF = solIF, solHH = solHH, parameters = parameters )
  }
  if( .Object@name == "integral" ){
    scoreValue = intMeasure( solIF = solIF, solHH = solHH, parameters = parameters )
  }
  return( new("Score", value = scoreValue, measure = .Object, parameters = parameters) )
}


# ________________________ OBJECT DEFINITONS CURRENTS ----

pinkNoise       = new("CurrentType", name = "pinkNoise"  , 
                      parameters = c("mu", "sigma") , 
                      default    = c( 2.5,   12.5 ) )
whiteNoise      = new("CurrentType", name = "whiteNoise" , 
                      parameters = c("mu", "sigma") , 
                      default    = c( 5  ,   5    ) )
constantCurrent = new("CurrentType", name = "constant", 
                      parameters = c("I"), 
                      default    = c(5) )
weiner          = new("CurrentType", name = "weiner",
                      parameters = c("sigma", "lowerLimit" , "upperLimit"),
                      default    = c(0.2    ,     0        ,      10     )
)
preset          = new("CurrentType", name = "preset",
                      parameters = character(),
                      default    = numeric()   )
OUprocess       = new("CurrentType", name = "Orstein-Uhlenbeck",
                      parameters = c("lambda", "nu", "sigma")  ,
                      default    = c(  0.001  ,  9  ,   0.5    ))




# ________________________ OBJECT DEFINITONS ALGORITHMS ----

GA1 = new("Algorithm", type = GAAlgorithm, parameterVals = c("popSize" = 100, "maxIter" = 25, "run" = 6), fixedParameters = NA_real_ )
GA2 = new("Algorithm", type = GAAlgorithm, parameterVals = c("popSize" = 100, "maxIter" = 25, "run" = 6), fixedParameters = c( V0 = -10 , Vr = 0, theta = 10.14 , spike = 100 , threshold = 10.14, refractoryPeriod = 0 ))
GA2shortNoRp = new("Algorithm", type = GAAlgorithm, parameterVals = c("popSize" = 100, "maxIter" = 25, "run" = 3), fixedParameters = c( V0 = -10 , Vr = 0, theta = 10.14 , spike = 100 , threshold = 10.14, refractoryPeriod = 0 ))
GA2Toy = new("Algorithm", type = GAAlgorithm, parameterVals = c("popSize" = 20, "maxIter" = 25, "run" = 6), fixedParameters = c( V0 = -10 , Vr = 0, theta = 10.14 , spike = 100 , threshold = 10.14, refractoryPeriod = 0 ))
GAtoy = new("Algorithm", type = GAAlgorithm, parameterVals = c("popSize" = 20, "maxIter" = 20, "run" = 6), fixedParameters = NA_real_)

Gain1         = new("Algorithm", type = GainFitAlgorithm , parameterVals = c( noise = F, resolution = 0.2 , minCurrent = 1, maxCurrent = 20, sigma = 11 , "maxIter" = 25, "popSize" = 30, run = 6) )
Gain1         = new("Algorithm", type = GainFitAlgorithm , parameterVals = c( noise = F, resolution = 0.2 , minCurrent = 1, maxCurrent = 20, sigma = 0.5 , "maxIter" = 25, "popSize" = 30, run = 6) ,
                    fixedParameters = c("V0" = -12))
Gain2         = new("Algorithm", type = GainFitAlgorithm , parameterVals = c( noise = F, resolution = 0.2 , minCurrent = 1, maxCurrent = 20, sigma = 0.5 , "maxIter" = 25, "popSize" = 50, run = 6) ,
                    fixedParameters = c("V0" = -12, "refractoryPeriod" = 0))
Gaintoy      = new("Algorithm", type = GainFitAlgorithm , parameterVals = c( noise = F, resolution = 0.3  , minCurrent = 1, maxCurrent = 15, sigma = 0.5 , "maxIter" = 25, "popSize" = 10, run = 10 ) )
Gaintoyfixed = new("Algorithm", type = GainFitAlgorithm , parameterVals = c( noise = F, resolution = 0.2  , minCurrent = 1, maxCurrent = 20, sigma = 0.5 , "maxIter" = 25, "popSize" = 15, run = 4 ) ,
                   fixedParameters = c( "V0" = -12 ))
GainNoise    = new("Algorithm", type = GainFitAlgorithm , parameterVals = c( noise = T, resolution = 0.1 , minCurrent = 1, maxCurrent = 15, sigma = 11 , "maxIter" = 25, "popSize" = 50, run = 5) )

GainFinal    = new("Algorithm", type = GainFitAlgorithm , parameterVals = c( noise = F, resolution = 0.1 , minCurrent = 0.1, maxCurrent = 30 , sigma = 0.5 ))

# ________________________ FUNCTIONS FOR RUNNING ALGORITHMS ----

runGA = function( algorithm , problem , model , mode = "both", shrinkSpaceFactor = 0.1 ){
  
  parametersEstimated = setdiff( model@parameters , names(algorithm@fixedParameters) )
  fixedParameters = algorithm@fixedParameters
  print( parametersEstimated  )
  
  if( (mode == "both" | mode == "sub") & model@name != "IF" ){
    parametersEstimated = setdiff( parametersEstimated , model@overPars )
    print(parametersEstimated)
    
    GAsub = ga( type = "real-valued", fitness = fitnessAll , problem = problem , model = model , mode = "sub" ,
                parameterNames = parametersEstimated, fixedParameters = algorithm@fixedParameters ,
                popSize = algorithm@parameterVals["popSize"], 
                maxiter = algorithm@parameterVals["maxIter"], 
                run     = algorithm@parameterVals["run"], parallel = FALSE,
                lower   = model@lowerBound[ parametersEstimated ] ,
                upper   = model@upperBound[ parametersEstimated ] )
    
    solution = summary(GAsub)$solution
    for(i in 1:length(solution[,1]) ){
      values =  c( solution[i,] , fixedParameters)
      values = values[ intersect( names(values) , model@parameters ) ]
      parameterPoint = new( "ParameterPoint", model = model, values =  values )
      scoreTe = fitnessAll( values , problem = problem, model = model, mode = "sub", testData = TRUE )
      estimate = new( "Estimate",
                      parameterPoint = parameterPoint ,
                      algorithm = algorithm ,
                      scoreTr = summary(GAsub)$fitness,
                      scoreTe = scoreTe               ,
                      GA      = GAsub
      )  
      resLength = length(problem@results[[model@name]][["sub"]])
      problem@results[[ model@name ]][[ "sub" ]][[ resLength + 1 ]] = estimate
    }
  }
  
  if( (mode == "both" | mode == "over") & !is.na( shrinkSpaceFactor ) ){
    
    subPars = problem@results[[ model@name ]][[ "sub" ]]
    if( model@name %in% names(SubModels) ){
      subPars = c(subPars, problem@results[[ SubModels[[model@name]]@name ]][[ "sub" ]]  ) }
    estimate = subPars[[1]]
    for( est in subPars ){
      if( est@scoreTe > estimate@scoreTe ){ estimate = est }
    }
    parameterSpace = (model@upperBound[ model@subPars ] - model@lowerBound[ model@subPars ] )/2
    lowerBound = c( model@lowerBound[ model@overPars ] , estimate@parameterPoint@values[model@subPars] - parameterSpace*shrinkSpaceFactor )
    upperBound = c( model@upperBound[ model@overPars ] , estimate@parameterPoint@values[model@subPars] + parameterSpace*shrinkSpaceFactor )
    
    lowerBound = lowerBound[parametersEstimated]
    upperBound = upperBound[parametersEstimated]
    
    
    GAover = ga( type = "real-valued", fitness = fitnessAll, problem = problem, model = model, mode = "over",
                 parameterNames = names(lowerBound), fixedParameters = algorithm@fixedParameters ,
                 popSize = algorithm@parameterVals["popSize"], 
                 maxiter = algorithm@parameterVals["maxIter"], 
                 run     = algorithm@parameterVals["run"], parallel = FALSE,
                 lower   = lowerBound ,
                 upper   = upperBound )
    
    solution = summary( GAover )$solution
    for(i in 1:length(solution[,1]) ){
      values =  c( solution[i,] , fixedParameters)
      values = values[ intersect( names(values) , model@parameters ) ]
      parameterPoint = new( "ParameterPoint", model = model, values = values )
      scoreTe = fitnessAll( values , problem = problem, model = model, mode = "over", testData = TRUE )
      estimate = new( "Estimate",
                      parameterPoint = parameterPoint ,
                      algorithm = algorithm ,
                      scoreTr = summary(GAover)$fitness,
                      scoreTe = scoreTe ,
                      GA = GAover
      )  
      resLength = length(problem@results[[model@name]][["over"]])
      problem@results[[ model@name ]][[ "over" ]][[ resLength + 1 ]] = estimate 
    }
  }
  
  if( (mode == "all") ){
    GAover = ga( type = "real-valued", fitness = fitnessAll, problem = problem, model = model, mode = "over",
                 parameterNames = parametersEstimated , fixedParameters = algorithm@fixedParameters ,
                 popSize = algorithm@parameterVals["popSize"], 
                 maxiter = algorithm@parameterVals["maxIter"], 
                 run     = algorithm@parameterVals["run"], parallel = FALSE,
                 lower   = model@lowerBound[ parametersEstimated ] ,
                 upper   = model@upperBound[ parametersEstimated ] )
    
    solution = summary( GAover )$solution
    for(i in 1:length(solution[,1]) ){
      values =  c( solution[i,] , fixedParameters)
      values = values[ intersect( names(values) , model@parameters ) ]
      parameterPoint = new( "ParameterPoint", model = model, values = values )
      scoreTe = fitnessAll( values , problem = problem, model = model, mode = "over", testData = TRUE )
      estimate = new( "Estimate",
                      parameterPoint = parameterPoint ,
                      algorithm = algorithm ,
                      scoreTr = summary(GAover)$fitness,
                      scoreTe = scoreTe , 
                      GA = GAover
      )  
      resLength = length(problem@results[[model@name]][["over"]])
      problem@results[[ model@name ]][[ "over" ]][[ resLength + 1 ]] = estimate 
    }
  }
  
  return(problem)
}

runGainFit = function( algorithm , problem , model ){
  evaluatedParameters = setdiff( model@parameters , names(algorithm@fixedParameters) )
  fixedParameters = algorithm@fixedParameters
  
  GA = ga( type = "real-valued", fitness = gainFitness , problem = problem , model = model , algorithm = algorithm,
           parameterNames = intersect(model@parameters,evaluatedParameters), fixedParameters = algorithm@fixedParameters, maxFitness = 0 ,
           popSize = algorithm@parameterVals["popSize"], 
           maxiter = algorithm@parameterVals["maxIter"], 
           run     = algorithm@parameterVals["run"], parallel = FALSE,
           lower   = model@lowerBound[evaluatedParameters] ,
           upper   = model@upperBound[evaluatedParameters] )
  
  solution = summary(GA)$solution
  for(i in 1:length(solution[,1]) ){
    values = solution[i,]
    values = c(values , fixedParameters[ intersect( names(fixedParameters) , model@parameters ) ] )
    parameterPoint = new( "ParameterPoint", model = model, values = values )
    scoreTe = fitnessAll( values , problem = problem, model = model, mode = "over", testData = TRUE )
    estimate = new( "Estimate",
                    parameterPoint = parameterPoint ,
                    algorithm = algorithm ,
                    scoreTr = summary(GA)$fitness,
                    scoreTe = scoreTe
    )  
    resLength = length(problem@results[[model@name]][["over"]])
    problem@results[[ model@name ]][[ "over" ]][[ resLength + 1 ]] = estimate
  }
  return(problem)
}


# ________________________ FUNCTIONS for plotting and saving ----

viridisShift = function(n){
  fl =  floor(n)
  rest = n - fl
  if(rest == 0){ return( viridis(fl) ) }
  scale = viridis( n/rest )
  return( scale[ c(1 , fl/rest)] )
}

safeSaveRDS = function( object, file ){
  
  if( file.exists(file) ){
    fileName = unlist(strsplit(file, split = "[.]"))[1]
    fileExt =unlist(strsplit(file, split = "[.]"))[2]
    i = 1
    while( file.exists( paste(fileName,"(", i ,")", ".", fileExt, sep = "") ) ){
      i = i + 1
    }
    file = paste(fileName,"(", i ,")", "." , fileExt, sep = "")
  }
  
  saveRDS( object, file )
  
}

plotSolutions = function( results, from = NA, to = NA ,
                          showPeaks = TRUE,
                          showThreshold = TRUE,
                          threshold = 10.14,
                          showLegend = FALSE,
                          showCurrent = TRUE,
                          current = NaN,
                          verticalShift = 0,
                          colourShift = 0.2,
                          colourScale = viridisShift){
  
  
  if( class(results) == "Solution" ){
    numOfResults = 1
    results = list(results)
  }else{
    numOfResults = length(results)
  }
  
  if(length(verticalShift) == 1){
    verticalShift = verticalShift*(0:(numOfResults-1) )
  }
  
  if(!is.na(current)){
    commonCurrent = T
    I = current@I
  }else{
    
    if( showCurrent ){
      
      commonCurrent = T
      I = results[[1]]@current@I
      for( i in 1:numOfResults){
        commonCurrent = commonCurrent & all(results[[i]]@current@I == I)
      }
      
    }else{
      
      commonCurrent = F
      showCurrent = F
      
    }
  }
  
  if(is.na(threshold)){
    threshold =  results[[1]]@model@default["threshold"]
  }
  if( !is.na(results[[1]]@parameters@values["threshold"]) ){ threshold = results[[1]]@parameters@values["threshold"] }
  Time= results[[1]]@current@time
  stepsize = results[[1]]@current@stepsize
  
  
  if( is.na(from) ){ from = 1             }else{ from = from/stepsize }
  if( is.na( to ) ){ to   = Time/stepsize }else{  to  =   to/stepsize }
  
  data = data.frame()
  for(i in 1:numOfResults){
    data1 = data.frame( time = (1:(results[[i]]@current@time/results[[i]]@current@stepsize)*results[[i]]@current@stepsize)[from:to] )
    data1["I"] = results[[i]]@current@I[from:to]
    data1["V"] = results[[i]]@V[from:to] + verticalShift[i]
    data1["model"] = results[[i]]@model@name
    data1["sol"] = i
    data = rbind(data,data1)
  }
  data$sol = as.factor(data$sol)
  
  neurons = ggplot(data = data) + 
    geom_line( aes(x = time, y = V, group = sol, colour = sol ), size = 1  ) +
    xlab(element_blank()) + 
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
    ylab("-E (mV)") +
    scale_color_manual(values = colourScale(numOfResults+colourShift))
  
  if( showThreshold ){ neurons = neurons + geom_hline(yintercept =  threshold, linetype = "dashed") }
  
  if( verticalShift[length(verticalShift)] > 5 ){ neurons = neurons + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank() )   }
  
  if(showLegend){
    neurons = neurons + scale_colour_manual(values=colourScale(numOfResults), 
                                            name="Model",
                                            breaks= levels(as.factor(data$sol)) ,
                                            labels= unique( data$model ) ) + 
      theme( legend.position = c(1,1), legend.justification = c(1,1))
  }else{
    neurons = neurons + theme(legend.position="none")
  }
  
  #geom_line(aes(x = V1*stepsize, y = VIF), size = 1, color = "brown3" ) +
  
  if(showPeaks){
    
    spikes = ggplot( data = data, aes(x = time, y = V) )+ geom_blank() +
      xlim(0, Time/stepsize )
    
    for(i in 1:numOfResults){
      spikes = spikes + geom_vline( xintercept =  results[[i]]@spikes, size = 1, color = colourScale(numOfResults+colourShift)[i])
    }
    
    spikes = spikes +
      xlab(element_blank())+
      ylab("Spikes")+
      theme(axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank() )
  }
  
  if(showCurrent){
    
  
      
      current = ggplot(data = data ) + 
        geom_line(aes(x = time, y = I, group = sol, color = sol), size = 1 ) + 
        xlab("time (ms)") +
        ylab("I (nA)") +
        theme(panel.grid.minor.y = element_blank(), legend.position="none") +
        scale_color_manual(values = rep(viridis(1) , numOfResults))
      
    
  }
  
  if(showCurrent){
    if(showPeaks){
      cowplot::plot_grid(neurons,spikes, current, nrow = 3, rel_heights = c(3,0.5,1), align = "v") 
    }else{
      cowplot::plot_grid(neurons, current, nrow = 2, rel_heights = c(3,1,1), align = "v") 
    }
  }else{
    if(showPeaks){
      cowplot::plot_grid(neurons,spikes, nrow = 2, rel_heights = c(3,0.5), align = "v") 
    }else{
      neurons
    }
  }
  
  #plot(1:( T /stepsize),spikeHH)
  #plot(1:(T/stepsize),spikeIF)
  #plot(1:(T/stepsize),spikeHH - spikeIF)
  #plot(triedAs,triedResets, col = rgb(sqrt(triedscores/min(triedscores)),0,0))
  
  
}

plotCompare = function( pars , current , model , blockSpikes = F, fast = T){
  
  if( class(pars) == "numeric"  ){ pars = new("ParameterPoint", model = model , values = pars) }
  if( class(pars) == "Estimate" ){ pars = pars@parameterPoint }
  solHH = solve( current, model = HH , blockSpikes = blockSpikes )
  solIF = solve(current , pars , model = model , blockSpikes = blockSpikes )
  print( coincMeasure(solIF, solHH, parameters = c(delta = 2) )  ) 
  plotSolutions( list( solIF ,  solHH ) , showPeaks = !blockSpikes)
  
}

plotEstimations = function( model , problem , i=1 , mode = "over" , data = "Tr" , from = NA , to = NA, fast = T){
  estimate = problem@results[[model@name]][[mode]][[i]]
  if( estimate@algorithm@type@name == "GA" ){
    
    if( mode == "over" & data == "Tr"   ){ dataHH = problem@dataTrOver }
    if( mode == "over" & data == "Te"   ){ dataHH = problem@dataTeOver }
    if( mode == "sub"  & data == "Tr"   ){ dataHH = problem@dataTrSub  }
    if( mode == "sub"  & data == "Te"   ){ dataHH = problem@dataTeSub  }
    if(                  data == "rand" ){ current = new("Current", type = pinkNoise, time = to) }else{ current = dataHH@current }
    solIF = solve( current = current , parameterPoint = estimate@parameterPoint , model = model , blockSpikes = (mode == "sub") )
    if( data == "rand" ){ dataHH = solve( current = current , model = HH , blockSpikes = (mode == "sub") ) }
    if( mode == "over" ){ print( coincMeasure(solIF = solIF, dataHH , parameters = c(delta = 2) ) ) }
    if( mode == "sub"  ){ print(   intMeasure(solIF = solIF, dataHH ) ) }
    plotSolutions( list( solIF , dataHH) , from = from, to = to)
  }
  if( estimate@algorithm@type@name == "gain" ){
    
    if( data == "Te" | data == "rand" ){
      if( mode == "over" ){ dataHH = problem@dataTeOver }
      if( mode == "sub"  ){ dataHH = problem@dataTeSub  }
      if( data == "rand" ){ current = new("Current", type = pinkNoise, time = to) }else{ current = dataHH@current }
      solIF = solve( current = current , parameterPoint = estimate@parameterPoint , model = model , blockSpikes = (mode == "sub") )
      if( data == "rand" ){ dataHH = solve( current = current , model = HH , blockSpikes = (mode == "sub") )  }
      if( mode == "over" ){ print( coincMeasure(solIF = solIF, dataHH , parameters = c(delta = 2) ) ) }
      if( mode == "sub"  ){ print(   intMeasure(solIF = solIF, dataHH ) ) }
      plotSolutions( list( solIF , dataHH) , from = from, to = to)
    }else{
      gain = calculateGain( estimate@parameterPoint , estimate@algorithm , estimate@parameterPoint@model , fast = fast )
      sigma = estimate@algorithm@parameterVals[["sigma"]]
      plot( 1:300/10  , HHGain[,sigma*2] , type = "l")
      resolution = estimate@algorithm@parameterVals["resolution"]
      ts = (estimate@algorithm@parameterVals["minCurrent"]/resolution ):( estimate@algorithm@parameterVals["maxCurrent"]/resolution )*resolution
      lines( ts , gain  )
    }
    
  }
}


# ________________________ OBJECT DEFINITONS PROBLEMS ----

ITrOver  = new("Current", type = OUprocess, time = 3000)
ITrSub   = new("Current", type = OUprocess, time = 3000)
ITeOver  = new("Current", type = OUprocess, time = 1000)
ITeSub   = new("Current", type = OUprocess, time = 1000)

dataTrOver = solve( current = ITrOver , model = HH , blockSpikes = F)
dataTrSub  = solve( current = ITrSub  , model = HH , blockSpikes = T)
dataTeOver = solve( current = ITeOver , model = HH , blockSpikes = F)
dataTeSub  = solve( current = ITeSub  , model = HH , blockSpikes = T)

problem = new("Problem", dataTrOver = dataTrOver, dataTrSub = dataTrSub, dataTeOver = dataTeOver , dataTeSub = dataTeSub )
problemPinkNoise = new("Problem", dataTrOver = dataTrOver, dataTrSub = dataTrSub, dataTeOver = dataTeOver , dataTeSub = dataTeSub )

ITrOverToy  = new("Current", type = OUprocess, time = 300)
ITrSubToy   = new("Current", type = OUprocess, time = 300)
ITeOverToy  = new("Current", type = OUprocess, time = 300)
ITeSubToy   = new("Current", type = OUprocess, time = 300)

dataTrOverToy = solve( current = ITrOverToy , model = HH , blockSpikes = F)
dataTrSubToy  = solve( current = ITrSubToy  , model = HH , blockSpikes = T)
dataTeOverToy = solve( current = ITeOverToy , model = HH , blockSpikes = F)
dataTeSubToy  = solve( current = ITeSubToy  , model = HH , blockSpikes = T)

problemToy = new("Problem", dataTrOver = dataTrOverToy, dataTrSub = dataTrSubToy, dataTeOver = dataTeOverToy , dataTeSub = dataTeSubToy )
problemToy2 = new("Problem", dataTrOver = dataTrOverToy, dataTrSub = dataTrSubToy, dataTeOver = dataTeOverToy , dataTeSub = dataTeSubToy )

ITrOverToyOU  = new("Current", type = OUprocess , time = 300)
ITrSubToyOU   = new("Current", type = OUprocess , time = 300)
ITeOverToyOU  = new("Current", type = OUprocess , time = 300)
ITeSubToyOU   = new("Current", type = OUprocess , time = 300)

dataTrOverToyOU = solve( current = ITrOverToyOU , model = HH , blockSpikes = F)
dataTrSubToyOU  = solve( current = ITrSubToyOU  , model = HH , blockSpikes = T)
dataTeOverToyOU = solve( current = ITeOverToyOU , model = HH , blockSpikes = F)
dataTeSubToyOU  = solve( current = ITeSubToyOU  , model = HH , blockSpikes = T)

problemToyOU = new("Problem", dataTrOver = dataTrOverToyOU , dataTrSub = dataTrSubToyOU, dataTeOver = dataTeOverToyOU , dataTeSub = dataTeSubToyOU )


# _______________________ RUNNING ALGORITHMS ----



problem = runGainFit( Gain1 , problem , LIF )
safeSaveRDS( problem, "finalProblem.RData")
problem = runGainFit( Gain1 , problem , aEIF )
safeSaveRDS( problem, "finalProblem.RData")
problem = runGainFit( Gain1 , problem , EIF )
safeSaveRDS( problem, "finalProblem.RData")
problem = runGainFit( Gain1 , problem , aQIF )
safeSaveRDS( problem, "finalProblem.RData")
problem = runGainFit( Gain1 , problem , QIF )
safeSaveRDS( problem, "finalProblem.RData")
problem = runGainFit( Gain1 , problem , aLIF )
safeSaveRDS( problem, "finalProblem.RData")
problem = runGainFit( Gain1 , problem , aIF )
safeSaveRDS( problem, "finalProblem.RData")

problem = runGainFit( Gain2 , problem , aEIF )
safeSaveRDS( problem, "finalProblem.RData")
problem = runGainFit( Gain2 , problem , EIF )
safeSaveRDS( problem, "finalProblem.RData")
problem = runGainFit( Gain2 , problem , aQIF )
safeSaveRDS( problem, "finalProblem.RData")
problem = runGainFit( Gain2 , problem , QIF )
safeSaveRDS( problem, "finalProblem.RData")
problem = runGainFit( Gain2 , problem , aLIF )
safeSaveRDS( problem, "finalProblem.RData")
problem = runGainFit( Gain2 , problem , LIF )
safeSaveRDS( problem, "finalProblem.RData")
problem = runGainFit( Gain2 , problem , aIF )
safeSaveRDS( problem, "finalProblem.RData")


problem =      runGA( GA1          , problem = problem, mode = "both", shrinkSpaceFactor = 0.2, model = aEIF )
safeSaveRDS( problem, "problem.RData" )
problem =      runGA( GA1          , problem = problem, mode = "both", shrinkSpaceFactor = 0.2, model = EIF )
safeSaveRDS( problem, "problem.RData" )
problem =      runGA( GA1          , problem = problem, mode = "both", shrinkSpaceFactor = 0.2, model = aLIF )
safeSaveRDS( problem, "problem.RData" )
problem =      runGA( GA1          , problem = problem, mode = "both", shrinkSpaceFactor = 0.2, model = LIF )
safeSaveRDS( problem, "problem.RData" )
problem =      runGA( GA1          , problem = problem, mode = "both", shrinkSpaceFactor = 0.2, model = aQIF )
safeSaveRDS( problem, "problem.RData" )
problem =      runGA( GA1          , problem = problem, mode = "both", shrinkSpaceFactor = 0.2, model = QIF )

problem =   runGainFit( Gaintoyfixed , problem , aEIF )
problem =   runGainFit( Gaintoy , problem , LIF )

safeSaveRDS( problemPinkNoise, "problem1.RData" )
problemPinkNoise = runGainFit( GainNoise , problem = problemPinkNoise , model = QIF )
safeSaveRDS( problemPinkNoise, "problemPinkNoise.RData" )
problemPinkNoise = runGainFit( GainNoise , problem = problemPinkNoise , model = EIF )
safeSaveRDS( problemPinkNoise, "problemPinkNoise.RData" )


problemToy =      runGA( GAtoy   , problem = problemToy, model = aEIF , mode = "both" , shrinkSpaceFactor = 0.1 )
problemToy =      runGA( GAtoy   , problem = problemToy, model = aEIF , mode = "both" , shrinkSpaceFactor = 0.2 )
problemToy =      runGA( GA2Toy  , problem = problemToy, model = aEIF , mode = "all"  , shrinkSpaceFactor = 0.1 )
problemToy = runGainFit( Gaintoy , problem = problemToy, model = LIF  )
problemToy = runGainFit( Gaintoyfixed , problem = problemToy, model = EIF  )

problemToyOU =      runGA( GA2Toy  , problem = problemToyOU, model = aEIF , mode = "all"  , shrinkSpaceFactor = 0.1 )

problemToy2 = runGA( GA1 , problemToy2 , aEIF , mode = "sub"  )
# _______________________ VIEWING RESULTS ----




I = new("Current", type = pinkNoise, time = 100, parameter = c( mu = 2.5, sigma = 12.5))
pars2 = new("ParameterPoint", model = aEIF, values = c("C" = 1 , "R" = 5, Vr = 0 , a = 0.025, b= 1, d = 6, delta = 1 , V0 = -10, theta = 6 ,"refractoryPeriod" = 0))
sol = solve( I, model = aQIF , parameterPoint = pars2 ,blockSpikes = F )
solHH = solve( I, model = HH  , blockSpikes =  F)
plotSolutions( list( sol, solHH) )

I = rep(0,20000)
I[5001:15000] = 4
I = new("Current", type = preset , I = I)
I = new("Current", type = OUprocess, time = 100 )
pars2 = new("ParameterPoint", model = aEIF, values = c("C" = 0.3, R = 4 , Vr = 0 , a = 0.005, b= 8, d = 6, delta = 2 , V0 = -12, theta = 10.14, spike = 100 ,"refractoryPeriod" = 0))
plotCompare( c("C" = 0.3, R = 4 , Vr = 0 , a = 0.005, b= 8, d = 6, delta = 2 , V0 = -12, theta = 10.14, spike = 100 ,"refractoryPeriod" = 0) , I , aEIF)


I = new("Current", type = pinkNoise, time = 100)
solHH = solve( I, model = HH  , blockSpikes =  F)
pars1 = new("ParameterPoint", model = aLIF, values = c("C" = 3.7 , "R" = 8, threshold = 6.15,  Vr = 0 , a = 0.02, b= 0.2, d = 1 , V0 = -10,"refractoryPeriod" = 0))
sol = solve( I, model = aLIF , parameterPoint = pars1 ,blockSpikes = F )
coincMeasure( sol , solHH , c(delta = 2) )
plotSolutions( list( sol, solHH) )


I = new("Current", type = pinkNoise, time = 100)
solHH = solve( I, model = HH  , blockSpikes =  F)
pars1 = new("ParameterPoint", model = LIF, values = c("C" = 3.7 , "R" = 8, threshold = 6.15,  Vr = 0 ,V0 = -10,"refractoryPeriod" = 8))
sol = solve( I, model = LIF , parameterPoint = pars1 ,blockSpikes = F )
plotSolutions( list( sol, solHH) )


I = new("Current", type = pinkNoise, time = 100)
solHH = solve( I, model = HH  , blockSpikes =  F)
pars1 = new("ParameterPoint", model = aEIF, values = aaaa )
sol = solve( I, model = aEIF , parameterPoint = pars1 ,blockSpikes = F )
coincMeasure( sol , solHH , c(delta = 2) )
plotSolutions( list( sol, solHH) )


I = problemToy@dataTrOver@current
solHH = problemToy@dataTrOver
aaaa = problemToy@results$aEIF$over[[33]]@parameterPoint
sol = solve( I , model = aEIF , parameterPoint = aaaa )
coincMeasure( sol , solHH , c(delta = 2) )
plotSolutions( list( sol, solHH) )

plotEstimations( aEIF , problemToy , i = 33 , mode = "over", data = "Te" , from = NA)
plotEstimations( aEIF , problemToyOU , i = 1 , mode = "over", data = "Tr" , from = NA)



plot(1:300/10, HHGain[,1], type = "l")
for( i in 1:3){
  j = i - 2
  pars = new("ParameterPoint", model = aEIF, values = c("C" = 3 , "R" = 8 + j, Vr = 0 , a = 0.2, b= 1.1, d = 5, delta = 1 , V0 = -10, theta = 10.14, spike = 65 ,"refractoryPeriod" = 7))
  gain = calculateGain( pars, Gaintoy , aEIF )
  lines(  (1/0.3):(15/0.3)*0.3 , gain , col = viridis(10)[i])
}


for( i in 1:3){
  j = i - 2
  pars = new("ParameterPoint", model = aEIF, values = c("C" = 3 + j , "R" = 8, Vr = 0 , a = 0.2, b= 1.1, d = 5, delta = 1 , V0 = -10, theta = 10.14, spike = 65 ,"refractoryPeriod" = 7))
  gain = calculateGain( pars, Gaintoy , aEIF )
  lines(  (1/0.3):(15/0.3)*0.3 , gain , col = viridis(10)[3+i])
}


for( i in 1:3){
  j = i - 2
  pars = new("ParameterPoint", model = aEIF, values = c("C" = 3 , "R" = 8, Vr = 0 + j , a = 0.2, b= 1.1, d = 5, delta = 1 , V0 = -10, theta = 10.14, spike = 65 ,"refractoryPeriod" = 7))
  gain = calculateGain( pars, Gaintoy , aEIF )
  lines(  (1/0.3):(15/0.3)*0.3 , gain , col = viridis(10)[6+i])
}

I = new("Current", type = weiner, time = 100)
solHH = solve( I, model = HH  , blockSpikes =  F)
sol = solve( I, model = aEIF , parameterPoint = pars ,blockSpikes = F )
coincMeasure( sol , solHH , c(delta = 2) )
plotSolutions( list( sol, solHH) )


plot( 1:300/10 , HHGain[,1] , type = "l" )
pars3 = new("ParameterPoint", model = QIF, values = c("C" = 4.34 , R = 9.25 , Vr = 6.71 ,V0 = -12, theta = 6.88, spike = 100 ,"refractoryPeriod" = 4.46))
gain = calculateGain( pars3 , Gaintoy , QIF , fast = T)
lines( (1/0.3):(15/0.3)*0.3 , gain  )


i = 3
model = QIF
current = new("Current", type = pinkNoise , time = 200, parameters = c("mu" = 2.5, "sigma" = 12.5))
resIF = solve( current , parameterPoint = problemPinkNoise@results[[model@name]]$over[[i]]@parameterPoint , model = aQIF)
problemPinkNoise@results[[model@name]]$over[[i]]@algorithm@type@name
problemPinkNoise@results[[model@name]]$over[[i]]@scoreTr
problemPinkNoise@results[[model@name]]$over[[i]]@scoreTe
resHH = solve( current , model = HH)
evaluate( coincidence , resIF , resHH )@value
plotSolutions( list( resIF, resHH ) )


problemPinkNoise = runGA( GA1 , problem = problemPinkNoise, model = aLIF , mode = "both" , shrinkSpaceFactor = 0.1 )
problemPinkNoise = runGA( GA1 , problem = problemPinkNoise, model = EIF  , mode = "both" , shrinkSpaceFactor = 0.1 )

problemPinkNoise = readRDS("problemPinkNoise(7).RData")


#problemPinkNoise = runGainFit( GainNoise , problem = problemPinkNoise , model = IF )
#safeSaveRDS( problemPinkNoise, "problemPinkNoise.RData" )

problemPinkNoise = runGainFit( GainNoise , problem = problemPinkNoise , model = LIF )

safeSaveRDS( problemPinkNoise, "problemPinkNoise.RData" )


plotEstimations( aQIF , problem1 , i = 1 , mode = "over", data = "Tr" , from = NA , to = 200)


# _____________________________  PLOTTING RESULTS ----

problem2 = readRDS( "problem16.RData")
problem3 = readRDS("finalProblem5.rdata")
problem4 = readRDS("finalProblem10.rdata")

parsIF       = problem4@results$IF$over[[2]]@parameterPoint     
parsLIF      = problem4@results$LIF$over[[2]]@parameterPoint  
parsQIF      = problem4@results$QIF$over[[2]]@parameterPoint  
parsQIFGain  = problem4@results$QIF$over[[2]]@parameterPoint 
parsEIF      = problem4@results$EIF$over[[2]]@parameterPoint  
parsaEIFGain = problem4@results$aEIF$over[[2]]@parameterPoint 
parsaEIF     = problem4@results$aEIF$over[[2]]@parameterPoint 

parsIF = new("ParameterPoint", values = c("C" = 1/0.15, V0 = -10), model = IF)
parsaIF = aIF@default 
parsLIF = new("ParameterPoint", values = c(C = 4.649357 , R = 4.528402 , threshold = 15.616254 , V0 = -12 ), model = LIF)
parsaLIF = aLIF@default 

parsQIF = problem2@results$QIF$over[[3]]@parameterPoint
parsQIFGain = new("ParameterPoint", values = c("C" = 4.34 , R = 9.25 , Vr = 6.71 ,V0 = -12, theta = 6.88, spike = 100 ,"refractoryPeriod" = 4.46), model = QIF)

parsaQIF = aQIF@default
parsEIF = problem2@results$EIF$over[[2]]@parameterPoint


parsaEIFGain = problem@results$aEIF$over[[1]]@parameterPoint
parsaEIF = problem2@results$aEIF$over[[100]]@parameterPoint

# HH ----

setwd("C:/Users/rapey/Desktop/Biomatematyka i Teoria Gier/R_IF/Figures")

I01 = new("Current", type = constantCurrent , parameters = c( I = 6.21  ) , time = 500)
I02 = new("Current", type = constantCurrent , parameters = c( I = 6.23  ) , time = 500)
I03 = new("Current", type = constantCurrent , parameters = c( I = 6.25  ) , time = 500)
I04 = new("Current", type = constantCurrent , parameters = c( I = 6.259 ) , time = 500)
Is = c(6.21, 6.23, 6.25, 6.259)
HH1 = solve( I01 , model = HH )
HH2 = solve( I02 , model = HH )
HH3 = solve( I03 , model = HH )
HH4 = solve( I04 , model = HH )
a = plotSolutions( list( HH1, HH2 , HH3 , HH4 ), verticalShift = (Is)*2000 , showPeaks = F , showThreshold = F, showCurrent = F , colourScale = function(n){ rep(viridis(1),n) } )
a
ggsave( filename = "HHSpikeTrains.jpeg", plot = a, device = "jpeg" )

I1 = new("Current", type = preset, I = c( rep(2     ,  500) ,
                                          rep(0     , 3000) ,
                                          rep(2.347 ,  500) ,
                                          rep(0     , 3000) ,
                                          rep(2.348 ,  500) ,
                                          rep(0     , 3000) ,
                                          rep(4     , 3500) ,
                                          rep(0     , 2000) , 
                                          rep(7     , 4000) ) , time = 200)


plotSolutions( solve(I1, model = HH), showPeaks = F) 
ggsave(filename = "FigureHHStep.jpeg", plot = plotSolutions( solve(I1, model = HH), showPeaks = F) , device = "jpeg")

I2 = new("Current", type = OUprocess , time = 200, parameters = c(  lambda =  0.001  , nu =  10  , sigma =  0.5    ))
safeSaveRDS( I2 , "Cuurent1.RData" )

plotSolutions( solve(I2, model = HH), showPeaks = F, showThreshold = F)
ggsave( filename = "FigureHHOUP.jpeg", plot = plotSolutions( solve(I2, model = HH), showPeaks = F) , device = "jpeg" )

# IF ----




I3 = new("Current", type = preset, I = c( rep(0     , 1000) , 
                                          rep(2     , 1000) ,
                                          rep(0     , 5000) , 
                                          rep(3     , 2500) ,
                                          rep(4.5     , 2500) ,
                                          rep(0     , 1000) ,
                                          rep(7     , 7000)
                                          ), time = 200)
  
plotSolutions( solve( I3, parsIF , model = IF ) )
ggsave( filename = "IFsteps.jpeg", plot = plotSolutions( solve( I3, parsIF , model = IF ) )  , device = "jpeg" )

plotCompare( parsIF, I2 , IF  )
ggsave( filename = "IFOUP.jpeg", plot = plotCompare( parsIF, I2 , IF  ) , device = "jpeg" )


# LIF ----

nonSpikingValueLIF = (parsLIF@values["threshold"])/parsLIF@values["R"]
I4 = new("Current", type = preset, I = c( rep(nonSpikingValueLIF    , 20000) , 
                                          rep(0     , 20000) 
                                                              ), time = 400)

plotSolutions( solve( I3, parsLIF , model = LIF ) )
ggsave( filename = "LIFsteps.jpeg", plot = plotSolutions( solve( I3, parsEIF , model = LIF ) )  , device = "jpeg" )

plotSolutions( solve( I4, parsLIF , model = LIF ) )
ggsave( filename = "LIFNoSpike.jpeg", plot =plotSolutions( solve( I4, parsLIF , model = LIF ) )  , device = "jpeg" )

plotCompare( parsLIF, I2 , LIF  )
ggsave( filename = "LIFOUP.jpeg", plot = plotCompare( parsIF, I2 , IF  ) , device = "jpeg" )

# EIF ----


plotSolutions( solve( I3, parsEIF , model = EIF ) )
ggsave( filename = "EIFsteps.jpeg"  , plot = plotSolutions( solve( I3, parsEIF , model = EIF ) )  , device = "jpeg" )


plotCompare( parsEIF, I2 , EIF  )
ggsave( filename = "EIFOUP.jpeg"    , plot = plotCompare( parsEIF, I2 , EIF  ) , device = "jpeg" )

parsEIFsub = problem2@results$EIF$sub[[1]]@parameterPoint

plotCompare( parsEIFsub, I2 , LIF , blockSpikes = T)
ggsave( filename = "EIFsubThreshold.jpeg", plot = plotCompare( parsEIFsub, I2 , EIF , blockSpikes = T) , device = "jpeg" )

# QIF ----


plotSolutions( solve( I3, parsQIF , model = QIF ) )
ggsave( filename = "QIFsteps.jpeg"  , plot = plotSolutions( solve( I3, parsQIF , model = QIF ) )  , device = "jpeg" )


plotCompare( parsQIF, I2 , QIF  )
ggsave( filename = "QIFOUP.jpeg"    , plot = plotCompare( parsQIF, I2 , QIF  ) , device = "jpeg" )

parsQIFsub = problem2@results$QIF$sub[[1]]@parameterPoint

plotCompare( parsQIFsub, I2 , QIF , blockSpikes = T)
ggsave( filename = "QIFsubThreshold.jpeg", plot = plotCompare( parsQIF, I2 , QIF , blockSpikes = T) , device = "jpeg" )

parsQIF = new("ParameterPoint", values = c("C" = 4.34 , R = 9.25 , Vr = 6.71 ,V0 = -12, theta = 6.88, spike = 100 ,"refractoryPeriod" = 4.46), model = QIF)

# aEIF ---- 


plotSolutions( solve( I3, parsaEIF , model = aEIF ) )
ggsave( filename = "aEIFsteps.jpeg"  , plot = plotSolutions( solve( I3, parsaEIF , model = aEIF ) )  , device = "jpeg" )


plotCompare( parsaEIF, I2 , aEIF  )
ggsave( filename = "aEIFOUP.jpeg"    , plot = plotCompare( parsaEIF, I2 , aEIF  ) , device = "jpeg" )

parsaEIFsub = problem2@results$aEIF$sub[[1]]@parameterPoint

plotCompare( parsaEIFsub, I2 , aEIF , blockSpikes = T)
ggsave( filename = "aEIFsubThreshold.jpeg", plot = plotCompare( parsaEIFsub, I2 , aEIF , blockSpikes = T) , device = "jpeg" )


# ALL ----

I5 = new("Current", type = preset, I = c( rep(13,2000) ,
                                          rep(13,2000) ,
                                          rep(0 ,2000) ,
                                          rep(10,2000) ,
                                          rep(9 ,2000) ,  
                                          rep(0 ,2000) , 
                                          rep(7 ,2000) ,
                                          rep(6 ,2000) ,
                                          rep(5 ,2000) , 
                                          rep(0 ,2000)  ) , time = 200)

solHH   = solve( I5 ,                model = HH    )
solIF   = solve( I5 , parsIF       , model = IF    )
solLIF  = solve( I5 , parsLIF      , model = LIF   )
solQIF  = solve( I5 , parsQIF      , model = QIF   )
solaEIF = solve( I5 , parsaEIFGain , model = aEIF  )

plotSolutions( list( solHH , solaEIF , solQIF , solLIF , solIF ) , verticalShift = 120 , showPeaks = F , showThreshold = F , colourScale = function(n){ rep(viridis(1),n) } )

ggsave( filename = "allSteps.jpeg", plot = plotSolutions( list( solHH , solaEIF , solQIF , solLIF , solIF ) , verticalShift = 120 , showPeaks = F , showThreshold = F , colourScale = function(n){ rep(viridis(1),n) } )  , device = "jpeg" )

I6 = new("Current", type = OUprocess, time = 200, parameters = c(lambda = 0.001, nu = 7 , sigma = 0.5))


solHH   = solve( I6 ,            model = HH    )
solIF   = solve( I6 , parsIF   , model = IF    )
solLIF  = solve( I6 , parsLIF  , model = LIF   )
solQIF  = solve( I6 , parsQIF  , model = QIF   )
solaEIF = solve( I6 , parsaEIF , model = aEIF  )

plotSolutions( list( solHH , solaEIF , solQIF , solLIF , solIF ), verticalShift = 120 , showPeaks = F , showThreshold = F , colourScale = function(n){ rep(viridis(1),n) } )

ggsave( filename = "allOUP.jpeg", plot =  plotSolutions( list( solHH , solaEIF , solQIF , solLIF , solIF ) , verticalShift = 120 , showPeaks = F , showThreshold = F , colourScale = function(n){ rep(viridis(1),n) } ), device = "jpeg" )

# Gain ----

gainIF   = calculateGain( parsIF   , GainFinal , IF   , fast = T )
gainLIF  = calculateGain( parsLIF  , GainFinal , LIF  , fast = T )
gainQIF  = calculateGain( parsQIF  , GainFinal , QIF  , fast = T )
gainEIF  = calculateGain( parsEIF  , GainFinal , EIF  , fast = T )
gainaEIF = calculateGain( parsaEIF , GainFinal , aEIF , fast = T )
plot( 1:300/10  , HHGain[,1] , type = "l" )
lines( 1:300/10 , gainIF   , col = viridis(1) )
lines( 1:300/10 , gainLIF  , col = viridis(2) )
lines( 1:300/10 , gainQIF  , col = viridis(3) )
lines( 1:300/10 , gainEIF  , col = viridis(4) )
lines( 1:300/10 , gainaEIF , col = viridis(5) )

t = 1:300/10
data = as.data.frame( cbind( t , HH = HHGain[,1], HH2 = HHGain[,22] , IF =  gainIF, LIF = gainLIF, QIF =  gainQIF, EIF =  gainEIF, aEIF = gainaEIF ) )
ggplot( data = data , aes(x = t) ) +
  geom_line( aes( y = HH   ) , size = 1 , colour = viridis(7)[1] ) + 
  geom_line( aes( y = HH2  ) , size = 1 , colour = viridis(7)[2] ) + 
  geom_line( aes( y = IF   ) , size = 1 , colour = viridis(7)[3] ) + 
  geom_line( aes( y = LIF  ) , size = 1 , colour = viridis(7)[4] ) + 
  geom_line( aes( y = QIF  ) , size = 1 , colour = viridis(7)[5] ) + 
  geom_line( aes( y = EIF  ) , size = 1 , colour = viridis(7)[6] ) +
  geom_line( aes( y = aEIF ) , size = 1 , colour = viridis(7)[7] ) +
  xlab( label = "I(nA)") + 
  ylab( label = "frequency (Hz)")+ 
  theme( legend.position = c(1,1), legend.justification = c(1,1))

data = as.data.frame( rbind( cbind(t , HHGain[,1] , "HH") , cbind(t , HHGain[,22], "HH noise") , cbind(t, gainIF, "IF") , cbind(t,gainLIF , "LIF"), cbind(t, gainQIF, "QIF"), cbind(t,gainEIF, "EIF"), cbind(t,gainaEIF, "aEIF") ) )
data$V3 = as.factor(data$V3)
ggplot( data = data ) + 
  geom_line(  aes(x = t, y = V2, group = V3 , colour = V3) , size = 1 ) + 
  theme( legend.position = c(1,1), legend.justification = c(1,1))


neurons = ggplot(data = data) + 
  geom_line( aes(x = time, y = V, group = sol, colour = sol ), size = 1  ) +
  xlab(element_blank()) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  ylab("-E (mV)") +
  scale_color_manual(values = colourScale(numOfResults+colourShift))

neurons = neurons + scale_colour_manual(values=colourScale(numOfResults), 
                                        name="Model",
                                        breaks= levels(as.factor(data$sol)) ,
                                        labels= unique( data$model ) ) + 
  theme( legend.position = c(1,1), legend.justification = c(1,1))
