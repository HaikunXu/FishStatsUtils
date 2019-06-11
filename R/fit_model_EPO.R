
#' Fit VAST to data
#'
#' \code{fit_model} fits a spatio-temporal model to data
#'
#' This function is the user-interface for the functions that determine the extrapolation-grid, define spatial objects, assemble data, build model, and estimate parameters.
#'
#' @param settings Output from \code{make_settings}
#' @inheritParams make_extrapolation_info
#' @inheritParams make_spatial_info
#' @inheritParams VAST::make_data
#' @inheritParams VAST::make_model
#' @inheritParams TMBhelper::fit_tmb
#' @param extrapolation_args tagged list of optional arguments to pass to \code{FishStatsUtils::make_extrapolation_info}
#' @param spatial_args tagged list of optional arguments to pass to \code{FishStatsUtils::make_spatial_info}
#' @param optimize_args tagged list of optional arguments to pass to \code{TMBhelper::Optimize}
#' @param model_args tagged list of optional arguments to pass to \code{VAST::make_model}
#' @param run_model Boolean indicating whether to run the model or simply return the inputs and built TMB object
#' @param ... additional parameters to pass to \code{VAST::make_data}
#'
#' @return Returns a tagged list of internal objects, the TMB object, and slot \code{parameter_estimates} containing the MLE estimates
#'
#' @family wrapper functions
#' @seealso \code{?VAST} for general documentation, \code{?make_settings} for generic settings, \code{?fit_model} for model fitting, and \code{?plot_results} for generic plots
#'
#' @examples
#' \dontrun{
#' # Load packages
#' library(TMB)
#' library(VAST)
#'
#' # load data set
#' # see `?load_example` for list of stocks with example data
#' # that are installed automatically with `FishStatsUtils`.
#' example = load_example( data_set="EBS_pollock" )
#'
#' # Make settings
#' settings = make_settings( n_x=50, Region=example$Region, purpose="index",
#'   strata.limits=example$strata.limits )
#'
#' # Run model
#' fit = fit_model( "settings"=settings, "Lat_i"=example$sampling_data[,'Lat'],
#'   "Lon_i"=example$sampling_data[,'Lon'], "t_i"=example$sampling_data[,'Year'],
#'   "c_i"=rep(0,nrow(example$sampling_data)), "b_i"=example$sampling_data[,'Catch_KG'],
#'   "a_i"=example$sampling_data[,'AreaSwept_km2'], "v_i"=example$sampling_data[,'Vessel'] )
#'
#' # Plot results
#' plot_results( settings=settings, fit=fit )
#' }
#'
#' @export
fit_model_EPO = function( settings, zone, Data, t_iz, c_iz, b_i, a_i,
  v_i=rep(0,length(b_i)), working_dir=paste0(getwd(),"/"),
  Xconfig_zcp=NULL, X_gtp=NULL, X_itp=NULL, Q_ik=NULL, newtonsteps=1,
  extrapolation_args=list(), spatial_args=list(), optimize_args=list(), model_args=list(),
  silent=TRUE, run_model=TRUE, ... ){

  # Local function -- combine two lists
  combine_lists = function( default, input ){
    output = default
    for( i in seq_along(input) ){
      if( names(input)[i] %in% names(default) ){
        output[[names(input)[i]]] = input[[i]]
      }else{
        output = c( output, input[i] )
      }
    }
    return( output )
  }

  Lat_i <- Data_Geostat[,'Lat']
  Lon_i <- Data_Geostat[,'Lon']
  
  # Assemble inputs
  data_frame = data.frame( "Lat_i"=Lat_i, "Lon_i"=Lon_i, "a_i"=a_i, "v_i"=v_i, "b_i"=b_i )
  # Decide which years to plot
  year_labels = seq( min(t_iz), max(t_iz) )
  years_to_plot = which( unique(t_iz) %in% sort(unique(t_iz)))

  # Save record
  dir.create(working_dir, showWarnings=FALSE, recursive=TRUE)
  save( settings, file=file.path(working_dir,"Record.RData"))
  capture.output( settings, file=file.path(working_dir,"Record.txt"))

  # Build extrapolation grid
  message("\n### Making extrapolation-grid")

  lat_North <- Data_Geostat$Lat > 0
  Data_Geostat_North <- Data_Geostat[lat_North,]
  Extrapolation_List_North = make_extrapolation_info(Region=Region, zone = zone, strata.limits=strata.limits, observations_LL=Data_Geostat_North[,c('Lat','Lon')], grid_dim_km=c(25,25) )
  
  lat_South <- Data_Geostat$Lat < 0
  Data_Geostat_South <- Data_Geostat[lat_South,]
  Extrapolation_List_South = make_extrapolation_info(Region=Region, zone = zone, strata.limits=strata.limits, observations_LL=Data_Geostat_South[,c('Lat','Lon')], grid_dim_km=c(25,25) )
  
  Extrapolation_List_South$Data_Extrap$N_km <- Extrapolation_List_South$Data_Extrap$N_km - 10000
  
  a_el <- rbind(Extrapolation_List_North$a_el,Extrapolation_List_South$a_el)
  Data_Extrap <- rbind(Extrapolation_List_North$Data_Extrap,Extrapolation_List_South$Data_Extrap)
  zone <- Extrapolation_List_North$zone
  flip_around_dateline <- Extrapolation_List_North$flip_around_dateline
  Area_km2_x <- c(Extrapolation_List_North$Area_km2_x,Extrapolation_List_South$Area_km2_x)
  
  extrapolation_list <- list("a_el" = a_el, "Data_Extrap" = Data_Extrap, "zone" = zone,
                             "flip_around_dateline" = flip_around_dateline,
                             "Area_km2_x" = Area_km2_x)

  spatial_list = make_spatial_info_EPO( grid_size_km=grid_size_km, n_x=n_x, Method=Method, Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=extrapolation_list, DirPath=DateFile, Save_Results=FALSE )
  plot_data(Extrapolation_List=extrapolation_list, Spatial_List=spatial_list, Data_Geostat=Data_Geostat )
  

  # Build data
  message("\n### Making data object") # VAST::
  data_list = VAST::make_data("Version"=settings$Version, "FieldConfig"=settings$FieldConfig, "OverdispersionConfig"=settings$OverdispersionConfig,
    "RhoConfig"=settings$RhoConfig, "ObsModel"=settings$ObsModel, "c_iz"=c_iz, "b_i"=b_i, "a_i"=a_i, "v_i"=v_i,
    "s_i"=spatial_list$knot_i-1, "t_iz"=t_iz, "spatial_list"=spatial_list, "Options"=settings$Options, "Aniso"=settings$use_anisotropy,
    Xconfig_zcp=Xconfig_zcp, X_gtp=X_gtp, X_itp=X_itp, Q_ik=Q_ik, ... )

  # Build object
  message("\n### Making TMB object")
  model_args = combine_lists( input=model_args, default=list("TmbData"=data_list, "RunDir"=working_dir, "Version"=settings$Version,
    "RhoConfig"=settings$RhoConfig, "loc_x"=spatial_list$loc_x, "Method"=spatial_list$Method) )
  tmb_list = do.call( what=VAST::make_model, args=model_args )  # VAST::
  if(silent==TRUE) tmb_list$Obj$env$beSilent()

  # Run the model or optionally don't
  if( run_model==FALSE ){
    # Build and output
    Return = list("data_frame"=data_frame, "extrapolation_list"=extrapolation_list, "spatial_list"=spatial_list,
      "data_list"=data_list, "tmb_list"=tmb_list, "year_labels"=year_labels, "years_to_plot"=years_to_plot,
      "settings"=settings, "extrapolation_args"=extrapolation_args, "model_args"=model_args)
    return(Return)
  }

  # Optimize object
  message("\n### Estimating parameters")
  optimize_args_phase1 = combine_lists( default=optimize_args, input=list(obj=tmb_list$Obj, lower=tmb_list$Lower, upper=tmb_list$Upper,
    savedir=working_dir, getsd=FALSE, newtonsteps=0, bias.correct=FALSE, quiet=TRUE,
    control=list(eval.max=10000,iter.max=10000,trace=1), loopnum=2) )
  parameter_estimates = do.call( what=TMBhelper::fit_tmb, args=optimize_args_phase1 )

  # Check fit of model (i.e., evidence of non-convergence based on bounds, approaching zero, etc)
  if(exists("check_fit")){
    problem_found = VAST::check_fit( parameter_estimates )
    if( problem_found==TRUE ){
      message("\n")
      stop("Please change model structure to avoid problems with parameter estimates and then re-try\n", call.=FALSE)
    }
  }

  # Restart estimates after checking parameters
  optimize_args_phase2 = combine_lists( input=optimize_args, default=list(obj=tmb_list$Obj, lower=tmb_list$Lower, upper=tmb_list$Upper,
    savedir=working_dir, bias.correct=settings$bias.correct, newtonsteps=newtonsteps,
    bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct=settings$vars_to_correct),
    control=list(eval.max=10000,iter.max=10000,trace=1), loopnum=1) )
  optimize_args_phase2 = combine_lists( input=list(startpar=parameter_estimates$par), default=optimize_args_phase2 )
  parameter_estimates = do.call( what=TMBhelper::fit_tmb, args=optimize_args_phase2 )

  # Extract standard outputs
  Report = tmb_list$Obj$report()
  ParHat = tmb_list$Obj$env$parList( parameter_estimates$par )

  # Build and output
  Return = list("data_frame"=data_frame, "extrapolation_list"=extrapolation_list, "spatial_list"=spatial_list,
    "data_list"=data_list, "tmb_list"=tmb_list, "parameter_estimates"=parameter_estimates, "Report"=Report,
    "ParHat"=ParHat, "year_labels"=year_labels, "years_to_plot"=years_to_plot, "settings"=settings,
    "extrapolation_args"=extrapolation_args, "model_args"=model_args, "optimize_args"=optimize_args)
  return( Return )
}
