# library_for_pin_depin
## description:
- this is just a library for me to record the codes used to finish my Ph.D project.
- the repository has a series of workflows to operate hoomd-blue; 
- and has a library to proceed trajectories of particles in 2D and get their structural or dynamic properties

## log:
- 20230504 add file_for_CUDA.test_cuda
- 20230424 edit workflow_analysis.show_disp_field.get_bicolor_disp.
- 20230423 add workflow_part.saveIndexCN346PCairoSeed. 
- 20230422 add workflow_analysis.show_cairo_order_parameter.
- 20230420 add proceed_file to operate files and make directories.
- 20230420 fix get_conditional_bonds_and_simplices, return count_polygon_relative
- 20230410 edit get_conditional_bonds_and_simplices, return count_polygon_relative
- 20230408 edit bond_plot_module init(ax=none)
- 20230404 add points_analysis_2d.count_polygon_n
- 20230329 add workflow_analysis.show_waiting_time_interstitial_motion
	   add points_analysis_2d.dynamical_facilitation_module.scan_displacement_t
- 20230328 add workflow_analysis.optimize_polyfit
	   add workflow_analysis.show_waiting_time_brownian
	   add workflow_analysis.show_waiting_time_dynamical_facilitation
- 20230324 add workflow_analysis.show_tuned_image to add line over image.
	   add bond_plot_module_for_image in points_analysis_2d
- 20230320 add workflow_uniform_continue in symmetry_transformation.simple_simulation
	   add workflow_analysis.get_msd_from_gsd()
	   add workflow_analysis.get_displacement_field_normal()
- 20230314 add workflow_analysis to record codes.
- 20230306 plot scatter vacancy and occupation in dynamical_facilitation_module.plot_reference_occupation
	   add dynamical_facilitation_module.plot_reference_positions_waiting_time_colored
- 20230225 set honeycomb3-8-12 as facilitation reference, part1 as traps.
- 20230223 add probability distribution of dynamical facilitation(ongoing)
	   add get_displacement_1D_overlap in data_analysis_cycle
	   add dynamical_facilitation_module in points_analysis_2d
- 20230222 add function_plot including example_plot and unit_transform
- 20230216 add workflow_mysql_to_data_pin_hex_to_honeycomb_part_klt_2m_random_oop
	   add mysql_data_processor in getDataAndScatter
	   add get_string_like_motion_rank in data_analysis_workflow, 
	   add plot_string_like_motion_rank in dynamic_points_analysis_2d
	   add get_displacement_field_xy_rank in displacemnt_field_2D
- 20230206 seperate draw_bonds_conditional_bond() into two parts( listed as follows) in bond_plot_module in points_analysis_2D
	   part1, restrict_axis_property()
	   part2, draw_points_with_conditional_bond()
- 20230131 add workflow_liquid in workflow_part
	   add init_cut in dynamic_points_analysis_2d.plot_bond_neighbor_change_oop
	   add __init_state_pin_from_gsd etc in pin_seed_oop
	   add lattice_constant in coordination_number module in save_from_gsd 
- 20230129 edit example in workflow_mysql_to_data_pin_hex_to_honeycomb_part_klt_2m()
	   pip install latex, code added in bond_plot_module in points_analysis_2D
- 20230119 edit workflow_part.workflow_simu_to_mysql_pin_hex_to_honeycomb_part_oop_klt_2m()
	   edit particle_tracking
- 20230118 edit data_analysis_cycle.data_analysis_workflow.get_displacment_field() 
	   generalize displacment_field_2D.get_displacement_field_xy() to any two frames.
	   add savetxt function in plot_hist_neighbor_change_event() in points_analysis_2D.dynamic_points_analysis_2d
	   add check mode in workflow_part.workflow_simu_to_mysql_pin_hex_to_honeycomb_part_oop_klt_2m()
	   add zorder to rank plot object in figure, in 
- 20230116 edit data_analysis.txyz_to_bond_plot,bond_plot_module
	   dynamic_points_analysis_2d.plot_bond_neighbor_change_oop,
	   static_points_analysis_2d.draw_bonds_conditional_bond_oop,
	   create class bond_plot_module in points_analysis_2D
- 20230107 simplify pin_seed_oop
	   edit data_analysis_cycle.data_analysis
	   edit workflow_part.workflow_simu_to_mysql_pin_hex_to_honeycomb_oop_klt_2m.
	   edit getDataAndScatter.
- 20230106 edit dynamic_points_analysis_2d.compute_nearest_neighbor_displacements & dynamic_coordination_bond_plot,
	   static_points_analysis_2d.get_first_minima_bond_length_distribution & draw_bonds_conditional_bond in points_analysis_2d_CUDA
	   fix error: 'ndarray' object cannot be interpreted as an integer
	   fix error: 'BlockManager' object has no attribute 'columns'
	   the main time cost is save numerous images.
- 20230105 fix bug of saving txyz as txyz_stable into 'txyz_stable.npy'
	   update get_displacement_field_xy.given txyz_stable, cut_edge is not necessary anymore.
	   move displacement_field_xy from proceed_gsd_file to dynamic_points_analysis_2d in points_analysis_2d
	   it seems that while-loop proceeding static_points_analysis_2d is too slow.
	   Implicit conversion to a NumPy array is not allowed. Please use `.get()` to construct a NumPy array explicitly.
	   https://cloud.tencent.com/developer/ask/sof/927734
- 20230104 reorganize structure of points_analysis_2D.dynamic_points_analysis_2d
	   edit proceed_gsd_file, dynamic_points_analysis_2d in points_analysis_2D
- 20230103 add monitor_neighbor_change_event(), bond_plot(), hist_neighbor_change_event() in points_analysis_2D.dynamic_points_analysis_2d
- 20230102 add sum_id_neighbors in static_points_analysis_2d.get_nearest_neighbor_dispalcement(), to check if neighbors of a center particle have changed over time.
- 20230101 set new standard of getDataAndScatter
- 20221219 invert the pixel coordination from image y-axis to coordination y-axis in particle_tracking.folder_frames_particle_tracking
	   test data_analysis_cycle.save_points_from_exp, right.
	   test and correct points_analysis_2d.dynamic_points_analysis_2d.compute_msd_normal, compute_atmsd_scan_t, compute_atmsd_t_chips
- 20221216 add data_analysis_cycle.save_points_from_exp, waiting test.
- 20221215 add Lindemann msd in dynamic_points_analysis_2d.
	   rename points_analysis_2D as dynamic_points_analysis_2d
- 20221213 add points_analysis_2D.msd.compute_nn_msd()
- 20221212 add introduction, mode='simu'or'exp', plot_msd,plot_trajectory, for points_analysis_2D.msd
- 20221209 edit points_analysis_2D.msd 
- 20221208 add kagome/honeycomb pin precise/low T
- 20221203 add seed in points_analysis_2D.proceed_gsd_file
	   add points_analysis_2D.msd ongoing
- 20221202 cut edge in workflow_part.workflow_simu_to_mysql_depin_from_kagome;
	   cut edge in workflow_part.workflow_simu_to_mysql_depin_from_honeycomb_part,
	   add RandomSeed in data_analysis_cycle.saveIndexPsi3Psi6;
	   cut edge in melt_hex_from_honeycomb;
	   add simple_simulation;
	   edit account in points_analysis_2D;
- 20221130 add psik_plot in data_analysis_cycle.save_from_gsd; 
	   add psi6 in workflow_mysql_to_data_pin_hex_to_honeycomb_part_klt_2m;
	   add condition="where kT = 1" in workflow_mysql_to_data_pin_hex_to_honeycomb_part_klt_2m(account='tplab')
- 20221129 point_analysis_2D mean repalced by median; symmetry_transformation/pin_seed_oop.py add mode=""
- 20221128 pin_seed_oop.py add self.mode
- 20221128 add file get_a_from_image
