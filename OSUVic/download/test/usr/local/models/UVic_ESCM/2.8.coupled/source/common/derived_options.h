! Do not edit the source file derived_options.h unless you are willing to
! remake all code. The mk script may not delete and remake all files that
! depend on derived options so files can become inconsistent.

#if defined pressure_gradient_average && !defined fourth_order_window
# define fourth_order_window
#endif
#if defined biharmonic && !defined fourth_order_window
# define fourth_order_window
#endif
#if defined fourth_order_tracer_advection && !defined fourth_order_window
# define fourth_order_window
#endif
#if defined fct && !defined fourth_order_window
# define fourth_order_window
#endif
#if defined quicker && !defined fourth_order_window
# define fourth_order_window
#endif

#if defined uvic_embm && !defined uvic_mom && !defined uvic_replacst
# define uvic_replacst
#endif

#if defined uvic_npzd || defined uvic_carbon_14 || defined uvic_o2
# if !defined source_term
#  define source_term
# endif
#endif

#if defined uvic_ice_cpts10 || defined uvic_ice_cpts5 || defined uvic_ice_cpts3 || defined uvic_ice_cpts1
# if !defined uvic_ice_cpts
#  define uvic_ice_cpts
# endif
#endif
#if defined uvic_ice_cpts
# if !defined uvic_ice_cpts10 && !defined uvic_ice_cpts5 && !defined uvic_ice_cpts3 && !defined uvic_ice_cpts1
#  define uvic_ice_cpts1
# endif
#endif

#if defined uvic_cfc
# if !defined uvic_cfc11
#  define uvic_cfc11
# endif
# if !defined uvic_cfc12
#  define uvic_cfc12
# endif
#endif

#if defined uvic_embm_CO2_lin && !defined uvic_embm_co2_lin
# define uvic_embm_co2_lin
#endif
#if defined uvic_embm_CO2_exp && !defined uvic_embm_co2_exp
# define uvic_embm_co2_exp
#endif

#if defined uvic_embm_crops && !defined uvic_embm_cropdata
# define uvic_embm_cropdata
#endif
#if defined uvic_embm_cropdata && !defined uvic_embm_crops
# define uvic_embm_crops
#endif

#if defined uvic_embm_crops_transient && !defined uvic_embm_cropdata_transient
# define uvic_embm_cropdata_transient
#endif
#if defined uvic_embm_cropdata_transient && !defined uvic_embm_crops_transient
# define uvic_embm_crops_transient
#endif

#if defined uvic_embm_crops_transient && !defined uvic_embm_crops
# define uvic_embm_crops
#endif
#if defined uvic_embm_cropdata_transient && !defined uvic_embm_cropdata
# define uvic_embm_cropdata
#endif
#if defined uvic_embm_solardata_transient && !defined uvic_embm_solardata
# define uvic_embm_solardata
#endif
#if defined uvic_embm_co2data_transient && !defined uvic_embm_co2data
# define uvic_embm_co2data
#endif
#if defined uvic_embm_icedata_transient && !defined uvic_embm_icedata
# define uvic_embm_icedata
#endif

#if defined llnl_plume_heat && !defined llnl_plume
# define llnl_plume
#endif
#if defined llnl_plume_salt && !defined llnl_plume
# define llnl_plume
#endif
#if defined llnl_plume_brine && !defined llnl_plume
# define llnl_plume
#endif
#if defined llnl_plume_all_heat && !defined llnl_plume
# define llnl_plume
#endif
#if defined llnl_plume_all_salt && !defined llnl_plume
# define llnl_plume
#endif
#if defined biharmonic_rm || defined redi_diffusion || defined gent_mcwilliams || defined save_density_terms
# if !defined isopycmix && !defined isoneutralmix
#  define isopycmix
# endif
# if !defined dm_taper
#  define gkw_taper
# endif
#else
# if !defined dm_taper && defined isopycmix
#  define gkw_taper
# endif
#endif
#if defined redi_diffusion
# if !defined full_tensor && !defined small_tensor
#  define small_tensor
# endif
#endif
#if defined gent_mcwilliams
# if !defined gm_advect && !defined gm_skew
#  define gm_skew
# endif
#endif
#if !defined sixth_order_window && !defined fourth_order_window && !defined second_order_window
# if defined biharmonic_rm || defined isotropic_mixed || defined parallel_1d || defined bbl_ag
#  if !defined fourth_order_window
#   define fourth_order_window
#  endif
# elif defined fourth_order_tracer_advection || defined fct || defined quicker
#  if !defined fourth_order_window
#   define fourth_order_window
#  endif
# elif defined tracer_horz_biharmonic || defined velocity_horz_biharmonic
#  if !defined fourth_order_window
#   define fourth_order_window
#  endif
# else
#  define second_order_window
# endif
#endif
